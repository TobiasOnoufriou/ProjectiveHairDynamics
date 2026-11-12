#include "PDSolver.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

PDSolver::PDSolver()
    : m_h(1.0 / 60.0),
      m_gravity_constant(9.81),
      m_stiffness_drag(5e4),
      m_damping_coefficient(0.02),
      m_iterations_per_step(10),
      m_max_substeps_per_frame(4),
      m_track_energy(false),
      m_parallel_local_step(false),
      m_drop_backlog_when_slow(true),
      m_row_offset(1, 0),
      m_system_ready(false),
      m_time_accumulator(0.0) {}

PDSolver::~PDSolver() = default;

void PDSolver::Clear() {
  m_constraints.clear();
  m_drag_constraint.reset();
  m_positions.resize(0, 3);
  m_velocities.resize(0, 3);
  m_initial_positions.resize(0, 3);
  m_masses.resize(0);
  m_free_index.clear();
  m_free_vertices.clear();
  m_fixed_vertices.clear();
  m_row_offset.assign(1, 0);
  m_weights.clear();
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

size_t PDSolver::numActiveConstraints() const {
  return m_constraints.size() + (m_drag_constraint ? 1 : 0);
}

size_t PDSolver::numLaidOut(const EigenMatrix3x& x) const {
  // Constraints stacked by the last RebuildSystemMatrix(), if it matches x.
  if (m_row_offset.empty() || m_j.cols() != x.rows()) return 0;
  return m_row_offset.size() - 1;
}

const Constraint& PDSolver::activeConstraint(size_t k) const {
  return k < m_constraints.size() ? *m_constraints[k] : *m_drag_constraint;
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

  // Stack every active constraint's G into J and record its rows and weight.
  const size_t num_active = numActiveConstraints();
  m_row_offset.assign(num_active + 1, 0);
  m_weights.resize(num_active);
  std::vector<SparseMatrixTriplet> j_triplets;
  j_triplets.reserve(num_active * 12);
  for (size_t k = 0; k < num_active; ++k) {
    const Constraint& c = activeConstraint(k);
    c.AddG(m_row_offset[k], &j_triplets);
    m_row_offset[k + 1] = m_row_offset[k] + c.NumRows();
    m_weights[k] = c.Weight();
  }
  const int rows = m_row_offset.back();
  m_j.resize(rows, n);
  m_j.setFromTriplets(j_triplets.begin(), j_triplets.end());
  m_gx.setZero(rows, 3);
  m_projections.setZero(rows, 3);

  // A = M/h^2 + J^T W J with W = per-row weights.
  VectorX w_rows(rows);
  for (size_t k = 0; k < num_active; ++k) {
    w_rows.segment(m_row_offset[k], m_row_offset[k + 1] - m_row_offset[k]).setConstant(m_weights[k]);
  }
  const SparseMatrix j_col(m_j);
  const SparseMatrix jt_w = SparseMatrix(j_col.transpose()) * w_rows.asDiagonal();
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

  if (nf == 0) {
    m_system_ready = true;  // Everything is pinned; nothing to solve.
    return;
  }
  m_llt.compute(a_ff);
  m_system_ready = (m_llt.info() == Eigen::Success);
  if (!m_system_ready) std::cerr << "PDSolver: system matrix factorization failed" << std::endl;
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

ScalarType PDSolver::localStep(const EigenMatrix3x& x, RowMatrix3x* gx, RowMatrix3x* p) const {
  const int num_active = static_cast<int>(numLaidOut(x));
  if (num_active == 0) return 0.0;
  gx->noalias() = m_j * x;
  // Each constraint reads and writes only its own rows: no races, and the same
  // projections whatever the thread count (only the energy sum is reordered).
  ScalarType energy = 0.0;
#pragma omp parallel for schedule(static) reduction(+ : energy) if (m_parallel_local_step)
  for (int k = 0; k < num_active; ++k) {
    const int first = m_row_offset[k];
    const int count = m_row_offset[k + 1] - first;
    activeConstraint(k).Project(gx->middleRows(first, count), p->middleRows(first, count));
    energy += 0.5 * m_weights[k] * (gx->middleRows(first, count) - p->middleRows(first, count)).squaredNorm();
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

void PDSolver::maxMetrics(const EigenMatrix3x& x, ScalarType* strain, ScalarType* volume_error) const {
  if (strain) *strain = 0.0;
  if (volume_error) *volume_error = 0.0;
  const int num = static_cast<int>(numLaidOut(x));
  if (num == 0) return;
  const RowMatrix3x gx = m_j * x;
  // Per-constraint values, then a serial max (MSVC's OpenMP 2.0 has no max reduction).
  std::vector<ScalarType> s(strain ? num : 0), v(volume_error ? num : 0);
#pragma omp parallel for schedule(static) if (m_parallel_local_step)
  for (int k = 0; k < num; ++k) {
    const int first = m_row_offset[k];
    const ConstRows rows = gx.middleRows(first, m_row_offset[k + 1] - first);
    if (strain) s[k] = activeConstraint(k).Strain(rows);
    if (volume_error) v[k] = activeConstraint(k).VolumeError(rows);
  }
  if (strain) *strain = *std::max_element(s.begin(), s.end());
  if (volume_error) *volume_error = *std::max_element(v.begin(), v.end());
}

ScalarType PDSolver::MaxStrain(const EigenMatrix3x& x) const {
  ScalarType strain = 0.0;
  maxMetrics(x, &strain, nullptr);
  return strain;
}

ScalarType PDSolver::MaxVolumeError(const EigenMatrix3x& x) const {
  ScalarType volume_error = 0.0;
  maxMetrics(x, nullptr, &volume_error);
  return volume_error;
}

void PDSolver::Step() {
  if (!m_system_ready) return;
  const auto start = std::chrono::steady_clock::now();
  const ScalarType inv_h2 = 1.0 / (m_h * m_h);

  // Inertial prediction with gravity folded in (f = m g, so M^-1 f = g):
  // y = x + h v + h^2 g. Fixed vertices stay at their (target) positions.
  EigenMatrix3x y = m_positions + m_h * m_velocities;
  y.col(1).array() -= m_h * m_h * m_gravity_constant;
  for (unsigned int i : m_fixed_vertices) y.row(i) = m_positions.row(i);

  // Constant part of the RHS on the free vertices: M/h^2 y - A_fc x_c. The known
  // positions x_c of the fixed vertices are moved to the RHS through A_fc.
  const size_t nf = m_free_vertices.size();
  EigenMatrix3x x_fixed(m_fixed_vertices.size(), 3);
  for (size_t k = 0; k < m_fixed_vertices.size(); ++k) x_fixed.row(k) = m_positions.row(m_fixed_vertices[k]);
  EigenMatrix3x b0_free(nf, 3);
  for (size_t k = 0; k < nf; ++k) {
    const unsigned int i = m_free_vertices[k];
    b0_free.row(k) = (m_masses[i] * inv_h2) * y.row(i);
  }
  if (!m_fixed_vertices.empty()) b0_free -= m_coupling * x_fixed;

  EigenMatrix3x x = y;
  EigenMatrix3x b_free(nf, 3);
  m_stats.energy.clear();
  for (int it = 0; it < m_iterations_per_step; ++it) {
    // Local step: each constraint projects into its own slot of P.
    const ScalarType constraint_energy = localStep(x, &m_gx, &m_projections);
    if (m_track_energy) m_stats.energy.push_back(inertiaEnergy(x, y) + constraint_energy);
    if (nf == 0) break;
    // Gather b = b0 + J^T W P, then the global step on the free vertices only;
    // one factor for x, y and z.
    b_free = b0_free;
    b_free.noalias() += m_gather * m_projections;
    const EigenMatrix3x x_free = m_llt.solve(b_free);
    for (size_t k = 0; k < nf; ++k) x.row(m_free_vertices[k]) = x_free.row(k);
  }
  if (m_track_energy) m_stats.energy.push_back(inertiaEnergy(x, y) + localStep(x, &m_gx, &m_projections));

  m_velocities = (x - m_positions) / m_h;
  m_velocities *= (1.0 - m_damping_coefficient);
  m_positions = x;

  m_stats.step_ms =
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
  maxMetrics(m_positions, &m_stats.max_strain, &m_stats.max_volume_error);
}
