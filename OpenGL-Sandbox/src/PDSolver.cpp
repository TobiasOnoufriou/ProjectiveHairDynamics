#include "PDSolver.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

PDSolver::PDSolver()
    : m_h(1.0 / 60.0),
      m_gravity_constant(9.81),
      m_stiffness_drag(5e4),
      m_damping_coefficient(0.02),
      m_iterations_per_step(10),
      m_max_substeps_per_frame(4),
      m_track_energy(false),
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

  // Full A = M/h^2 + sum_i w_i G_i^T G_i, as triplets.
  std::vector<SparseMatrixTriplet> triplets;
  triplets.reserve(n + m_constraints.size() * 4 + 1);
  const ScalarType inv_h2 = 1.0 / (m_h * m_h);
  for (unsigned int i = 0; i < n; ++i) triplets.emplace_back(i, i, m_masses[i] * inv_h2);
  for (const auto& c : m_constraints) c->AddToLHS(&triplets);
  if (m_drag_constraint) m_drag_constraint->AddToLHS(&triplets);

  // Split into the free block A_ff and the coupling A_fc. Rows of fixed
  // vertices are dropped: their positions are known.
  std::vector<SparseMatrixTriplet> ff, fc;
  ff.reserve(triplets.size());
  for (const SparseMatrixTriplet& t : triplets) {
    const int r = m_free_index[t.row()];
    if (r < 0) continue;
    const int c = m_free_index[t.col()];
    if (c >= 0) {
      ff.emplace_back(r, c, t.value());
    } else {
      fc.emplace_back(r, fixed_index[t.col()], t.value());
    }
  }

  SparseMatrix a_ff(nf, nf);
  a_ff.setFromTriplets(ff.begin(), ff.end());  // Duplicates are summed.
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
  while (m_time_accumulator >= m_h) {
    Step();
    m_time_accumulator -= m_h;
  }
}

ScalarType PDSolver::localStep(const EigenMatrix3x& x, EigenMatrix3x* b) const {
  ScalarType energy = 0.0;
  for (const auto& c : m_constraints) energy += c->AddToRHS(x, b);
  if (m_drag_constraint) energy += m_drag_constraint->AddToRHS(x, b);
  return energy;
}

ScalarType PDSolver::inertiaEnergy(const EigenMatrix3x& x, const EigenMatrix3x& y) const {
  ScalarType e = 0.0;
  for (unsigned int i : m_free_vertices) e += m_masses[i] * (x.row(i) - y.row(i)).squaredNorm();
  return 0.5 * e / (m_h * m_h);
}

ScalarType PDSolver::ConstraintEnergy(const EigenMatrix3x& x) const {
  EigenMatrix3x scratch = EigenMatrix3x::Zero(x.rows(), 3);
  return localStep(x, &scratch);
}

ScalarType PDSolver::MaxStrain(const EigenMatrix3x& x) const {
  ScalarType worst = 0.0;
  for (const auto& c : m_constraints) worst = std::max(worst, c->Strain(x));
  return worst;
}

void PDSolver::Step() {
  if (!m_system_ready) return;
  const auto start = std::chrono::steady_clock::now();
  const unsigned int n = NumVertices();
  const ScalarType inv_h2 = 1.0 / (m_h * m_h);

  // Inertial prediction with gravity folded in (f = m g, so M^-1 f = g):
  // y = x + h v + h^2 g. Fixed vertices stay at their (target) positions.
  EigenMatrix3x y = m_positions + m_h * m_velocities;
  y.col(1).array() -= m_h * m_h * m_gravity_constant;
  for (unsigned int i : m_fixed_vertices) y.row(i) = m_positions.row(i);

  // Constant part of the RHS: M/h^2 y.
  EigenMatrix3x b0 = (m_masses * inv_h2).asDiagonal() * y;

  // Known positions of the fixed vertices, moved to the RHS through A_fc.
  EigenMatrix3x x_fixed(m_fixed_vertices.size(), 3);
  for (size_t k = 0; k < m_fixed_vertices.size(); ++k) x_fixed.row(k) = m_positions.row(m_fixed_vertices[k]);
  const EigenMatrix3x coupling_rhs = m_coupling * x_fixed;

  const size_t nf = m_free_vertices.size();
  EigenMatrix3x x = y;
  EigenMatrix3x b(n, 3);
  EigenMatrix3x b_free(nf, 3);
  m_stats.energy.clear();
  for (int it = 0; it < m_iterations_per_step; ++it) {
    // Local step: each constraint projects independently.
    b = b0;
    const ScalarType constraint_energy = localStep(x, &b);
    if (m_track_energy) m_stats.energy.push_back(inertiaEnergy(x, y) + constraint_energy);
    if (nf == 0) break;
    // Global step on the free vertices only; one factor for x, y and z.
    for (size_t k = 0; k < nf; ++k) b_free.row(k) = b.row(m_free_vertices[k]);
    b_free -= coupling_rhs;
    const EigenMatrix3x x_free = m_llt.solve(b_free);
    for (size_t k = 0; k < nf; ++k) x.row(m_free_vertices[k]) = x_free.row(k);
  }
  if (m_track_energy) m_stats.energy.push_back(inertiaEnergy(x, y) + ConstraintEnergy(x));

  m_velocities = (x - m_positions) / m_h;
  m_velocities *= (1.0 - m_damping_coefficient);
  m_positions = x;

  m_stats.step_ms =
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
  m_stats.max_strain = MaxStrain(m_positions);
}
