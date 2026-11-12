#include "PDSolver.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <utility>
#include <vector>

#include "TetConstraints.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

// Small dense types for the Anderson normal equations (no heap allocation).
constexpr int kMaxWindow = PDSolver::kMaxAndersonWindow;
using SmallMatrix = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, 0, kMaxWindow, kMaxWindow>;
using SmallVector = Eigen::Matrix<ScalarType, Eigen::Dynamic, 1, 0, kMaxWindow, 1>;

// Frobenius inner product of two N x 3 matrices.
ScalarType Dot(const EigenMatrix3x& a, const EigenMatrix3x& b) { return a.cwiseProduct(b).sum(); }

}  // namespace

PDSolver::PDSolver()
    : m_h(1.0 / 60.0),
      m_gravity_constant(9.81),
      m_stiffness_drag(5e4),
      m_stiffness_contact(1e5),
      m_damping_coefficient(0.02),
      m_iterations_per_step(10),
      m_max_substeps_per_frame(4),
      m_track_energy(false),
      m_parallel_local_step(false),
      m_drop_backlog_when_slow(true),
      m_acceleration(PDAcceleration::kNone),
      m_chebyshev_rho(0.9),
      m_chebyshev_gamma(0.9),
      m_anderson_window(5),
      m_aa_count(0),
      m_aa_have_prev(false),
      m_system_ready(false),
      m_time_accumulator(0.0) {}

PDSolver::~PDSolver() = default;

void PDSolver::Clear() {
  m_constraints.clear();
  m_drag_constraint.reset();
  m_colliders.clear();
  m_contact_vertices.clear();
  m_contact_slots.clear();
  m_contact_stiffness = nullptr;
  m_positions.resize(0, 3);
  m_velocities.resize(0, 3);
  m_initial_positions.resize(0, 3);
  m_masses.resize(0);
  m_free_index.clear();
  m_free_vertices.clear();
  m_fixed_vertices.clear();
  m_batches.clear();
  m_j.resize(0, 0);
  m_gather.resize(0, 0);
  m_system_ready = false;
  m_time_accumulator = 0.0;
  m_stats = PDStepStats();
}

void PDSolver::SetVertices(const EigenMatrix3x& positions, const VectorX& masses) {
  m_positions = positions;
  m_velocities.setZero(positions.rows(), 3);
  m_masses = masses;
  m_free_index.assign(positions.rows(), 0);
  m_system_ready = false;
}

void PDSolver::SetFixed(unsigned int vertex, bool fixed) {
  if (vertex >= NumVertices()) return;
  m_free_index[vertex] = fixed ? -1 : 0;  // Real indices are assigned on rebuild.
  m_system_ready = false;
}

void PDSolver::AddConstraint(std::unique_ptr<Constraint> constraint) {
  m_constraints.push_back(std::move(constraint));
  m_system_ready = false;
}

void PDSolver::Finalize() {
  m_initial_positions = m_positions;
  m_time_accumulator = 0.0;
  RebuildSystemMatrix();
}

void PDSolver::SetFixedTarget(unsigned int vertex, const EigenVector3& target) {
  if (vertex >= NumVertices() || !IsFixed(vertex)) return;
  m_positions.row(vertex) = target.transpose();
}

int PDSolver::AddCollider(const Collider& collider) {
  m_colliders.push_back(collider);
  m_system_ready = false;  // The first collider adds the contact slots.
  return static_cast<int>(m_colliders.size()) - 1;
}

void PDSolver::SetContactVertices(const std::vector<unsigned int>& vertices) {
  m_contact_vertices = vertices;
  m_system_ready = false;
}

void PDSolver::RebuildSystemMatrix() {
  const unsigned int n = NumVertices();

  // Number the free and fixed vertices.
  m_free_vertices.clear();
  m_fixed_vertices.clear();
  std::vector<int> fixed_index(n, -1);
  for (unsigned int i = 0; i < n; ++i) {
    if (m_free_index[i] < 0) {
      fixed_index[i] = static_cast<int>(m_fixed_vertices.size());
      m_fixed_vertices.push_back(i);
    } else {
      m_free_index[i] = static_cast<int>(m_free_vertices.size());
      m_free_vertices.push_back(i);
    }
  }
  const int nf = static_cast<int>(m_free_vertices.size());
  const int nc = static_cast<int>(m_fixed_vertices.size());

  // Group the active constraints by kind and pair each volume constraint with a
  // strain constraint on the same tet (same G), then stack each batch's G into
  // J and copy its projection data into flat arrays. Rows are in batch order.
  std::vector<const Constraint*> springs, strains, volumes, fibres, generics;
  for (const auto& c : m_constraints) {
    switch (c->Kind()) {
      case ConstraintKind::kSpring: springs.push_back(c.get()); break;
      case ConstraintKind::kTetStrain: strains.push_back(c.get()); break;
      case ConstraintKind::kTetVolume: volumes.push_back(c.get()); break;
      case ConstraintKind::kTetFibre: fibres.push_back(c.get()); break;
      case ConstraintKind::kGeneric: generics.push_back(c.get()); break;
    }
  }
  if (m_drag_constraint) generics.push_back(m_drag_constraint.get());
  std::vector<std::pair<const TetStrainConstraint*, const TetVolumeConstraint*>> pairs;
  {
    std::map<std::array<unsigned int, 4>, std::vector<size_t>> strain_by_tet;
    for (size_t i = 0; i < strains.size(); ++i) {
      strain_by_tet[static_cast<const TetConstraint*>(strains[i])->Tet()].push_back(i);
    }
    std::vector<bool> paired(strains.size(), false);
    std::vector<const Constraint*> unpaired_volumes;
    for (const Constraint* c : volumes) {
      const auto* vol = static_cast<const TetVolumeConstraint*>(c);
      const auto it = strain_by_tet.find(vol->Tet());
      bool found = false;
      if (it != strain_by_tet.end()) {
        for (size_t i : it->second) {
          const auto* str = static_cast<const TetStrainConstraint*>(strains[i]);
          if (paired[i] || str->G() != vol->G()) continue;
          pairs.emplace_back(str, vol);
          paired[i] = true;
          found = true;
          break;
        }
      }
      if (!found) unpaired_volumes.push_back(c);
    }
    std::vector<const Constraint*> unpaired_strains;
    for (size_t i = 0; i < strains.size(); ++i) {
      if (!paired[i]) unpaired_strains.push_back(strains[i]);
    }
    strains.swap(unpaired_strains);
    volumes.swap(unpaired_volumes);
  }

  m_batches.clear();
  std::vector<SparseMatrixTriplet> j_triplets;
  j_triplets.reserve(m_constraints.size() * 12 + 2);
  std::vector<ScalarType> w_rows;  // W, one weight per row of J.
  int row = 0;
  // Appends c's rows to J and W.
  auto add_rows = [&](const Constraint* c) {
    c->AddG(row, &j_triplets);
    w_rows.insert(w_rows.end(), c->NumRows(), c->Weight());
    row += c->NumRows();
  };
  auto new_batch = [&](BatchKind kind, int rows_per_element) -> ConstraintBatch& {
    m_batches.emplace_back();
    m_batches.back().kind = kind;
    m_batches.back().first_row = row;
    m_batches.back().rows_per_element = rows_per_element;
    return m_batches.back();
  };
  if (!springs.empty()) {
    ConstraintBatch& b = new_batch(BatchKind::kSpring, 1);
    for (const Constraint* c : springs) {
      b.weight.push_back(c->Weight());
      b.rest_length.push_back(static_cast<const SpringConstraint*>(c)->GetRestLength());
      add_rows(c);
    }
  }
  if (!pairs.empty()) {
    ConstraintBatch& b = new_batch(BatchKind::kTetStrainVolume, 6);
    for (const auto& pr : pairs) {
      b.weight.push_back(pr.first->Weight());
      b.weight2.push_back(pr.second->Weight());
      b.param0.push_back(pr.first->MinStretch());
      b.param1.push_back(pr.first->MaxStretch());
      b.param2.push_back(pr.second->MinVolume());
      b.param3.push_back(pr.second->MaxVolume());
      add_rows(pr.first);
      add_rows(pr.second);
    }
  }
  if (!strains.empty()) {
    ConstraintBatch& b = new_batch(BatchKind::kTetStrain, 3);
    for (const Constraint* c : strains) {
      const auto* t = static_cast<const TetStrainConstraint*>(c);
      b.weight.push_back(c->Weight());
      b.param0.push_back(t->MinStretch());
      b.param1.push_back(t->MaxStretch());
      add_rows(c);
    }
  }
  if (!volumes.empty()) {
    ConstraintBatch& b = new_batch(BatchKind::kTetVolume, 3);
    for (const Constraint* c : volumes) {
      const auto* t = static_cast<const TetVolumeConstraint*>(c);
      b.weight.push_back(c->Weight());
      b.param0.push_back(t->MinVolume());
      b.param1.push_back(t->MaxVolume());
      add_rows(c);
    }
  }
  if (!fibres.empty()) {
    ConstraintBatch& b = new_batch(BatchKind::kTetFibre, 1);
    for (const Constraint* c : fibres) {
      const auto* t = static_cast<const TetFibreConstraint*>(c);
      b.weight.push_back(c->Weight());
      b.param0.push_back(t->Activation());
      b.param1.push_back(t->MaxContraction());
      b.fibre.push_back(t->Fibre());
      add_rows(c);
    }
  }
  // Contact slots: one row G = e_i per free contact vertex, whenever colliders exist.
  m_contact_slots.clear();
  if (!m_colliders.empty()) {
    if (m_contact_vertices.empty()) {
      m_contact_slots = m_free_vertices;
    } else {
      for (unsigned int i : m_contact_vertices) {
        if (i < n && m_free_index[i] >= 0) m_contact_slots.push_back(i);
      }
      std::sort(m_contact_slots.begin(), m_contact_slots.end());
      m_contact_slots.erase(std::unique(m_contact_slots.begin(), m_contact_slots.end()), m_contact_slots.end());
    }
  }
  if (!m_contact_slots.empty()) {
    ConstraintBatch& b = new_batch(BatchKind::kContact, 1);
    b.colliders = &m_colliders;
    m_contact_weight = ContactStiffness();
    const ScalarType w_contact = m_contact_weight;
    for (unsigned int i : m_contact_slots) {
      b.weight.push_back(w_contact);
      j_triplets.emplace_back(row, static_cast<int>(i), 1.0);
      w_rows.push_back(w_contact);
      ++row;
    }
  }
  if (!generics.empty()) {
    ConstraintBatch& b = new_batch(BatchKind::kGeneric, 0);
    for (const Constraint* c : generics) {
      b.weight.push_back(c->Weight());
      b.generic.push_back(c);
      b.generic_rows.push_back(row);
      add_rows(c);
    }
    b.generic_rows.push_back(row);
  }
  const int rows = row;
  m_j.resize(rows, n);
  m_j.setFromTriplets(j_triplets.begin(), j_triplets.end());

  // A = M/h^2 + J^T W J with W = per-row weights.
  const Eigen::Map<const VectorX> w(w_rows.data(), rows);
  const SparseMatrix j_col(m_j);
  const SparseMatrix jt_w = SparseMatrix(j_col.transpose()) * w.asDiagonal();
  SparseMatrix a = jt_w * j_col;
  const ScalarType inv_h2 = 1.0 / (m_h * m_h);
  {
    std::vector<SparseMatrixTriplet> mass;
    mass.reserve(n);
    for (unsigned int i = 0; i < n; ++i) mass.emplace_back(i, i, m_masses[i] * inv_h2);
    SparseMatrix m(n, n);
    m.setFromTriplets(mass.begin(), mass.end());
    a += m;
  }

  // Split into the free block A_ff and the coupling A_fc. Rows of fixed
  // vertices are dropped: their positions are known. The gather matrix keeps
  // the free rows of J^T W.
  std::vector<SparseMatrixTriplet> ff, fc, gather;
  ff.reserve(a.nonZeros());
  for (int col = 0; col < a.outerSize(); ++col) {
    for (SparseMatrix::InnerIterator it(a, col); it; ++it) {
      const int r = m_free_index[it.row()];
      if (r < 0) continue;
      const int c = m_free_index[it.col()];
      if (c >= 0) {
        ff.emplace_back(r, c, it.value());
      } else {
        fc.emplace_back(r, fixed_index[it.col()], it.value());
      }
    }
  }
  for (int col = 0; col < jt_w.outerSize(); ++col) {
    for (SparseMatrix::InnerIterator it(jt_w, col); it; ++it) {
      const int r = m_free_index[it.row()];
      if (r >= 0) gather.emplace_back(r, static_cast<int>(it.col()), it.value());
    }
  }
  m_gather.resize(nf, rows);
  m_gather.setFromTriplets(gather.begin(), gather.end());

  SparseMatrix a_ff(nf, nf);
  a_ff.setFromTriplets(ff.begin(), ff.end());
  m_coupling.resize(nf, nc);
  m_coupling.setFromTriplets(fc.begin(), fc.end());
  allocateBuffers();

  if (nf == 0) {
    m_system_ready = true;  // Everything is pinned; nothing to solve.
    return;
  }
  m_llt.compute(a_ff);
  m_system_ready = (m_llt.info() == Eigen::Success);
  if (!m_system_ready) std::cerr << "PDSolver: system matrix factorization failed" << std::endl;
}

void PDSolver::allocateBuffers() {
  const Eigen::Index n = NumVertices();
  const Eigen::Index nf = static_cast<Eigen::Index>(m_free_vertices.size());
  const Eigen::Index rows = m_j.rows();
  m_gx.setZero(rows, 3);
  m_projections.setZero(rows, 3);
  m_x.setZero(n, 3);
  m_y.setZero(n, 3);
  m_x_fixed.setZero(static_cast<Eigen::Index>(m_fixed_vertices.size()), 3);
  for (EigenMatrix3x* m : {&m_b0, &m_b, &m_solve_tmp, &m_u, &m_u_prev, &m_g, &m_g_accepted, &m_f, &m_f_prev,
                           &m_g_prev}) {
    m->setZero(nf, 3);
  }
  m_aa_df.clear();  // Sized on first use by andersonUpdate().
  m_aa_dg.clear();
  m_aa_count = 0;
  m_aa_have_prev = false;
}

void PDSolver::BeginDrag(unsigned int vertex, const EigenVector3& target) {
  if (vertex >= NumVertices() || IsFixed(vertex)) return;
  m_drag_constraint = std::make_unique<AttachmentConstraint>(&m_stiffness_drag, vertex, target);
  RebuildSystemMatrix();
}

void PDSolver::UpdateDrag(const EigenVector3& target) {
  if (m_drag_constraint) m_drag_constraint->SetFixedPoint(target);
}

void PDSolver::EndDrag() {
  if (!m_drag_constraint) return;
  m_drag_constraint.reset();
  RebuildSystemMatrix();
}

int PDSolver::DraggedVertex() const {
  return m_drag_constraint ? static_cast<int>(m_drag_constraint->GetConstrainedVertexIndex()) : -1;
}

void PDSolver::Reset() {
  EndDrag();
  m_positions = m_initial_positions;
  m_velocities.setZero(m_positions.rows(), 3);
  m_time_accumulator = 0.0;
}

void PDSolver::SetState(const EigenMatrix3x& positions, const EigenMatrix3x& velocities) {
  if (positions.rows() != m_positions.rows() || velocities.rows() != m_positions.rows()) return;
  m_positions = positions;
  m_velocities = velocities;
}

void PDSolver::Update(float dt) {
  if (!m_system_ready) return;
  // Clamp very long frames (e.g. window drag) so we don't spiral.
  m_time_accumulator += std::min(static_cast<double>(dt), m_h * m_max_substeps_per_frame);
  const auto start = std::chrono::steady_clock::now();
  int steps = 0;
  while (m_time_accumulator >= m_h) {
    Step();
    m_time_accumulator -= m_h;
    ++steps;
    // Steps cost more wall time than the time they simulate: more sub-steps
    // would only lengthen this frame and grow the next frame's backlog. Drop
    // it and run in slow motion instead.
    const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    if (m_drop_backlog_when_slow && elapsed > steps * m_h) {
      m_time_accumulator = 0.0;
      break;
    }
  }
}

int PDSolver::MaxThreads() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

// Called by every thread of a parallel region (or serially): the orphaned
// "omp for" splits the batch's elements over the team. Each element reads and
// writes only its own rows, so no two threads touch the same data and the
// projections do not depend on the thread count (only the energy sum is
// reordered). nowait: the next batch writes other rows, so no barrier is needed.
ScalarType PDSolver::projectBatch(const ConstraintBatch& b, const RowMatrix3x& gx, RowMatrix3x* p) {
  const int n = b.Size();
  // Row-major rows x 3: row r starts at 3 r.
  const ScalarType* g = gx.data() + 3 * static_cast<Eigen::Index>(b.first_row);
  ScalarType* q = p->data() + 3 * static_cast<Eigen::Index>(b.first_row);
  ScalarType energy = 0.0;
  switch (b.kind) {
    case BatchKind::kSpring:
#pragma omp for schedule(static) nowait
      for (int e = 0; e < n; ++e) {
        const Eigen::Map<const EigenVector3> d(g + 3 * e);
        Eigen::Map<EigenVector3> out(q + 3 * e);
        out = SpringConstraint::ProjectSpring(d, b.rest_length[e]);
        energy += 0.5 * b.weight[e] * (d - out).squaredNorm();
      }
      break;
    case BatchKind::kTetStrainVolume:
      // Rows 0-2 are F^T for the strain constraint, rows 3-5 the same F^T for
      // the volume constraint: one SVD, two projections.
#pragma omp for schedule(static) nowait
      for (int e = 0; e < n; ++e) {
        const Eigen::Map<const EigenMatrix3> f(g + 18 * e);
        Eigen::Map<EigenMatrix3> out_strain(q + 18 * e);
        Eigen::Map<EigenMatrix3> out_volume(q + 18 * e + 9);
        EigenMatrix3 u, v;
        EigenVector3 s;
        TetConstraint::SignedSVD(f, &u, &s, &v);
        out_strain = u * TetConstraint::ClampStretch(s, *b.param0[e], *b.param1[e]).asDiagonal() * v.transpose();
        out_volume =
            u * TetConstraint::ProjectVolumeSingularValues(s, *b.param2[e], *b.param3[e]).asDiagonal() * v.transpose();
        energy += 0.5 * (b.weight[e] * (f - out_strain).squaredNorm() + b.weight2[e] * (f - out_volume).squaredNorm());
      }
      break;
    case BatchKind::kTetStrain:
    case BatchKind::kTetVolume: {
      // The 3 rows F^T (row-major) are F column-major.
      const bool strain = b.kind == BatchKind::kTetStrain;
#pragma omp for schedule(static) nowait
      for (int e = 0; e < n; ++e) {
        const Eigen::Map<const EigenMatrix3> f(g + 9 * e);
        Eigen::Map<EigenMatrix3> out(q + 9 * e);
        out = strain ? TetConstraint::ProjectStrainF(f, *b.param0[e], *b.param1[e])
                     : TetConstraint::ProjectVolumeF(f, *b.param0[e], *b.param1[e]);
        energy += 0.5 * b.weight[e] * (f - out).squaredNorm();
      }
      break;
    }
    case BatchKind::kTetFibre:
#pragma omp for schedule(static) nowait
      for (int e = 0; e < n; ++e) {
        const Eigen::Map<const EigenVector3> fa(g + 3 * e);
        Eigen::Map<EigenVector3> out(q + 3 * e);
        const ScalarType target = TetFibreConstraint::TargetStretch(*b.param0[e], *b.param1[e]);
        out = TetFibreConstraint::ProjectFibre(fa, b.fibre[e], target);
        energy += 0.5 * b.weight[e] * (fa - out).squaredNorm();
      }
      break;
    case BatchKind::kContact:
      // Closest point outside every collider (the vertex itself if outside).
#pragma omp for schedule(static) nowait
      for (int e = 0; e < n; ++e) {
        const Eigen::Map<const EigenVector3> x(g + 3 * e);
        Eigen::Map<EigenVector3> out(q + 3 * e);
        EigenVector3 projected = x;
        ProjectOutOfColliders(*b.colliders, &projected);
        out = projected;
        energy += 0.5 * b.weight[e] * (x - out).squaredNorm();
      }
      break;
    case BatchKind::kGeneric:
#pragma omp for schedule(static) nowait
      for (int e = 0; e < n; ++e) {
        const int first = b.generic_rows[e];
        const int count = b.generic_rows[e + 1] - first;
        b.generic[e]->Project(gx.middleRows(first, count), p->middleRows(first, count));
        energy += 0.5 * b.weight[e] * (gx.middleRows(first, count) - p->middleRows(first, count)).squaredNorm();
      }
      break;
  }
  return energy;
}

void PDSolver::metricsBatch(const ConstraintBatch& b, const RowMatrix3x& gx, ScalarType* strain,
                            ScalarType* volume_error) {
  const int n = b.Size();
  const ScalarType* g = gx.data() + 3 * static_cast<Eigen::Index>(b.first_row);
  ScalarType s = *strain, v = *volume_error;
  switch (b.kind) {
    case BatchKind::kSpring:
#pragma omp for schedule(static) nowait
      for (int e = 0; e < n; ++e) {
        s = std::max(s, SpringConstraint::SpringStrain(Eigen::Map<const EigenVector3>(g + 3 * e), b.rest_length[e]));
      }
      break;
    case BatchKind::kTetStrainVolume:
    case BatchKind::kTetStrain:
    case BatchKind::kTetVolume: {
      // A pair's two F are the same: evaluate it once.
      const int stride = b.rows_per_element * 3;
#pragma omp for schedule(static) nowait
      for (int e = 0; e < n; ++e) {
        const Eigen::Map<const EigenMatrix3> f(g + stride * e);
        s = std::max(s, TetConstraint::StrainOfF(f));
        v = std::max(v, TetConstraint::VolumeErrorOfF(f));
      }
      break;
    }
    case BatchKind::kTetFibre:
#pragma omp for schedule(static) nowait
      for (int e = 0; e < n; ++e) {
        const ScalarType target = TetFibreConstraint::TargetStretch(*b.param0[e], *b.param1[e]);
        s = std::max(s, TetFibreConstraint::FibreStrain(Eigen::Map<const EigenVector3>(g + 3 * e), target));
      }
      break;
    case BatchKind::kContact:
      break;  // Penetration is reported by contactVelocityPass().
    case BatchKind::kGeneric:
#pragma omp for schedule(static) nowait
      for (int e = 0; e < n; ++e) {
        const int first = b.generic_rows[e];
        const ConstRows rows = gx.middleRows(first, b.generic_rows[e + 1] - first);
        s = std::max(s, b.generic[e]->Strain(rows));
        v = std::max(v, b.generic[e]->VolumeError(rows));
      }
      break;
  }
  *strain = s;
  *volume_error = v;
}

ScalarType PDSolver::localStep(const EigenMatrix3x& x, RowMatrix3x* gx, RowMatrix3x* p) const {
  if (!laidOutFor(x)) return 0.0;
  gx->noalias() = m_j * x;
  ScalarType energy = 0.0;
  const int num_batches = static_cast<int>(m_batches.size());
  // One parallel region for all batches (not one per batch).
#pragma omp parallel if (m_parallel_local_step) reduction(+ : energy)
  {
    for (int k = 0; k < num_batches; ++k) energy += projectBatch(m_batches[k], *gx, p);
  }
  return energy;
}

ScalarType PDSolver::inertiaEnergy(const EigenMatrix3x& x, const EigenMatrix3x& y) const {
  ScalarType e = 0.0;
  for (unsigned int i : m_free_vertices) e += m_masses[i] * (x.row(i) - y.row(i)).squaredNorm();
  return 0.5 * e / (m_h * m_h);
}

ScalarType PDSolver::ConstraintEnergy(const EigenMatrix3x& x) const {
  RowMatrix3x gx(m_j.rows(), 3), p(m_j.rows(), 3);
  return localStep(x, &gx, &p);
}

void PDSolver::maxMetrics(const EigenMatrix3x& x, RowMatrix3x* gx, ScalarType* strain,
                          ScalarType* volume_error) const {
  ScalarType s_max = 0.0, v_max = 0.0;
  if (laidOutFor(x)) {
    gx->noalias() = m_j * x;
    const int num_batches = static_cast<int>(m_batches.size());
    // Per-thread maxima, then combined (MSVC's OpenMP 2.0 has no max reduction).
#pragma omp parallel if (m_parallel_local_step)
    {
      ScalarType s = 0.0, v = 0.0;
      for (int k = 0; k < num_batches; ++k) metricsBatch(m_batches[k], *gx, &s, &v);
#pragma omp critical(pd_solver_metrics)
      {
        s_max = std::max(s_max, s);
        v_max = std::max(v_max, v);
      }
    }
  }
  if (strain) *strain = s_max;
  if (volume_error) *volume_error = v_max;
}

ScalarType PDSolver::MaxStrain(const EigenMatrix3x& x) const {
  RowMatrix3x gx(m_j.rows(), 3);
  ScalarType strain = 0.0;
  maxMetrics(x, &gx, &strain, nullptr);
  return strain;
}

ScalarType PDSolver::MaxVolumeError(const EigenMatrix3x& x) const {
  RowMatrix3x gx(m_j.rows(), 3);
  ScalarType volume_error = 0.0;
  maxMetrics(x, &gx, nullptr, &volume_error);
  return volume_error;
}

void PDSolver::globalSolve(const EigenMatrix3x& b, EigenMatrix3x* x) {
  // A_ff = P^T L L^T P. Same as m_llt.solve(b), but SimplicialLLT::solve
  // permutes its result in place, which allocates a temporary.
  if (m_llt.permutationP().size() > 0) {
    m_solve_tmp = m_llt.permutationP() * b;
  } else {
    m_solve_tmp = b;
  }
  m_llt.matrixL().solveInPlace(m_solve_tmp);
  m_llt.matrixU().solveInPlace(m_solve_tmp);
  if (m_llt.permutationPinv().size() > 0) {
    *x = m_llt.permutationPinv() * m_solve_tmp;
  } else {
    *x = m_solve_tmp;
  }
}

void PDSolver::scatterFree(const EigenMatrix3x& u, EigenMatrix3x* x) const {
  for (size_t k = 0; k < m_free_vertices.size(); ++k) x->row(m_free_vertices[k]) = u.row(k);
}

void PDSolver::andersonUpdate(const EigenMatrix3x& g, const EigenMatrix3x& f, EigenMatrix3x* u_next) {
  const int m = std::min(std::max(m_anderson_window, 1), kMaxAndersonWindow);
  if (static_cast<int>(m_aa_df.size()) != m || m_aa_df.front().rows() != g.rows()) {
    // First use or window changed: (re)size the history and restart.
    m_aa_df.assign(m, EigenMatrix3x::Zero(g.rows(), 3));
    m_aa_dg.assign(m, EigenMatrix3x::Zero(g.rows(), 3));
    m_aa_count = 0;
    m_aa_have_prev = false;
  }
  if (m_aa_have_prev) {
    const int slot = m_aa_count % m;
    m_aa_df[slot] = f - m_f_prev;
    m_aa_dg[slot] = g - m_g_prev;
    ++m_aa_count;
  }
  m_f_prev = f;
  m_g_prev = g;
  m_aa_have_prev = true;

  *u_next = g;
  const int cols = std::min(m_aa_count, m);
  if (cols == 0) return;
  // theta = argmin ||f - dF theta||^2 via the (tiny) normal equations, with a
  // relative Tikhonov term against nearly dependent columns.
  SmallMatrix normal(cols, cols);
  SmallVector rhs(cols);
  for (int i = 0; i < cols; ++i) {
    rhs(i) = Dot(m_aa_df[i], f);
    for (int j = 0; j <= i; ++j) normal(i, j) = normal(j, i) = Dot(m_aa_df[i], m_aa_df[j]);
  }
  const ScalarType trace = normal.trace();
  if (!(trace > 0.0)) return;
  normal.diagonal().array() += 1e-10 * trace;
  const SmallVector theta = normal.ldlt().solve(rhs);
  for (int i = 0; i < cols; ++i) *u_next -= theta(i) * m_aa_dg[i];
}

void PDSolver::contactVelocityPass() {
  m_stats.num_contacts = 0;
  m_stats.max_penetration = 0.0;
  const ScalarType w_contact = m_contact_weight;  // As factored.
  for (unsigned int i : m_contact_slots) {
    const EigenVector3 x = m_positions.row(i).transpose();
    ScalarType depth = 0.0;
    EigenVector3 n;
    const int k = DeepestPenetration(m_colliders, x, &depth, &n);
    if (k < 0) continue;
    ++m_stats.num_contacts;
    m_stats.max_penetration = std::max(m_stats.max_penetration, depth);

    // Velocity relative to the collider's surface, split into normal and
    // tangential parts. Approaching (vn < 0) is removed: inelastic contact.
    const Collider& c = m_colliders[k];
    const EigenVector3 vc = c.PointVelocity(x);
    const EigenVector3 rel = m_velocities.row(i).transpose() - vc;
    const ScalarType vn = std::max(rel.dot(n), 0.0);
    EigenVector3 vt = rel - rel.dot(n) * n;
    // Coulomb: the slot pushes with normal force f_n = w * depth over the
    // step, which can remove at most mu f_n h / m of tangential velocity
    // (static friction if that covers all of it).
    const ScalarType vt_norm = vt.norm();
    const ScalarType budget = c.friction * w_contact * depth * m_h / m_masses[i];
    vt *= vt_norm > budget ? 1.0 - budget / vt_norm : 0.0;
    m_velocities.row(i) = (vc + vn * n + vt).transpose();
  }
}

void PDSolver::Step() {
  if (!m_system_ready) return;
  const auto start = std::chrono::steady_clock::now();
  const ScalarType inv_h2 = 1.0 / (m_h * m_h);
  const size_t nf = m_free_vertices.size();
  const bool accelerate = m_acceleration != PDAcceleration::kNone;

  // Inertial prediction with gravity folded in (f = m g, so M^-1 f = g):
  // y = x + h v + h^2 g. Fixed vertices stay at their (target) positions.
  m_y = m_positions + m_h * m_velocities;
  m_y.col(1).array() -= m_h * m_h * m_gravity_constant;
  for (size_t k = 0; k < m_fixed_vertices.size(); ++k) {
    const unsigned int i = m_fixed_vertices[k];
    m_y.row(i) = m_positions.row(i);
    m_x_fixed.row(k) = m_positions.row(i);
  }

  // Everything on the RHS that does not change over the iterations, on the
  // free vertices: b0 = M/h^2 y - A_fc x_c (the known positions x_c of the
  // fixed vertices are moved to the RHS through A_fc).
  for (size_t k = 0; k < nf; ++k) {
    const unsigned int i = m_free_vertices[k];
    m_b0.row(k) = (m_masses[i] * inv_h2) * m_y.row(i);
  }
  if (!m_fixed_vertices.empty()) m_b0.noalias() -= m_coupling * m_x_fixed;

  // Initial guess x_0 = y.
  m_x = m_y;
  for (size_t k = 0; k < nf; ++k) m_u.row(k) = m_y.row(m_free_vertices[k]);
  m_aa_count = 0;
  m_aa_have_prev = false;
  int chebyshev_k = 0;  // Plain iterates since the last (re)start.
  ScalarType omega = 1.0;
  ScalarType e_accepted = std::numeric_limits<ScalarType>::infinity();

  m_stats.energy.clear();
  m_stats.energy.reserve(m_iterations_per_step + 1);
  m_stats.rejected_iterations = 0;
  using Clock = std::chrono::steady_clock;
  Clock::duration local_time{0}, global_time{0};
  for (int it = 0; it < m_iterations_per_step; ++it) {
    // Local step: each constraint projects into its own slot of P.
    const auto t0 = Clock::now();
    ScalarType e = localStep(m_x, &m_gx, &m_projections);
    if (accelerate || m_track_energy) e += inertiaEnergy(m_x, m_y);
    if (accelerate && e > e_accepted) {
      // The extrapolation overshot: fall back to the plain iterate from the
      // last accepted x (whose energy is <= e_accepted) and restart.
      m_u = m_g_accepted;
      scatterFree(m_u, &m_x);
      e = localStep(m_x, &m_gx, &m_projections) + inertiaEnergy(m_x, m_y);
      chebyshev_k = 0;
      m_aa_count = 0;
      m_aa_have_prev = false;
      ++m_stats.rejected_iterations;
    }
    const auto t1 = Clock::now();
    local_time += t1 - t0;
    e_accepted = e;
    if (m_track_energy) m_stats.energy.push_back(e);
    if (nf == 0) break;

    // Gather b = b0 + J^T W P, then the global step on the free vertices only;
    // one factor for x, y and z. g = G(x_k) is the plain PD iterate.
    m_b = m_b0;
    m_b.noalias() += m_gather * m_projections;
    globalSolve(m_b, &m_g);
    global_time += Clock::now() - t1;

    const bool last = it + 1 == m_iterations_per_step;
    if (!accelerate || last) {
      m_u = m_g;
    } else if (m_acceleration == PDAcceleration::kChebyshev) {
      m_g_accepted = m_g;
      ++chebyshev_k;
      const ScalarType rho2 = m_chebyshev_rho * m_chebyshev_rho;
      omega = chebyshev_k == 1 ? 1.0 : chebyshev_k == 2 ? 2.0 / (2.0 - rho2) : 4.0 / (4.0 - rho2 * omega);
      // m_f is scratch here. With omega = 1 the x_{k-1} terms cancel.
      m_f = omega * (m_chebyshev_gamma * (m_g - m_u) + m_u - m_u_prev) + m_u_prev;
      m_u_prev.swap(m_u);
      m_u.swap(m_f);
    } else {
      m_g_accepted = m_g;
      m_f = m_g - m_u;
      andersonUpdate(m_g, m_f, &m_u_prev);  // m_u_prev is scratch here.
      m_u.swap(m_u_prev);
    }
    scatterFree(m_u, &m_x);
  }
  if (m_track_energy) m_stats.energy.push_back(inertiaEnergy(m_x, m_y) + localStep(m_x, &m_gx, &m_projections));

  m_velocities = (m_x - m_positions) / m_h;
  m_positions = m_x;
  contactVelocityPass();
  m_velocities *= (1.0 - m_damping_coefficient);

  m_stats.step_ms = std::chrono::duration<double, std::milli>(Clock::now() - start).count();
  m_stats.local_ms = std::chrono::duration<double, std::milli>(local_time).count();
  m_stats.global_ms = std::chrono::duration<double, std::milli>(global_time).count();
  maxMetrics(m_positions, &m_gx, &m_stats.max_strain, &m_stats.max_volume_error);
}
