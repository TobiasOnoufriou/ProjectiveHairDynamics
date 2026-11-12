#include "Simulation.h"

#include <algorithm>
#include <iostream>

Simulation::Simulation()
    : m_h(1.0 / 60.0),
      m_particle_mass(1.0),
      m_gravity_constant(9.81),
      m_stiffness_attachment(1e6),
      m_stiffness_stretch(1e6),
      m_stiffness_bending(5e2),
      m_stiffness_drag(5e4),
      m_damping_coefficient(0.02),
      m_iterations_per_step(10),
      m_max_substeps_per_frame(4),
      m_drag_constraint(nullptr),
      m_system_ready(false),
      m_time_accumulator(0.0) {}

Simulation::~Simulation() { clearConstraints(); }

void Simulation::clearConstraints() {
  for (Constraint* c : m_constraints) delete c;
  m_constraints.clear();
  delete m_drag_constraint;
  m_drag_constraint = nullptr;
}

void Simulation::InitHairStrands(const std::vector<EigenVector3>& roots,
                                 const EigenVector3& direction,
                                 unsigned int segments,
                                 ScalarType segment_length) {
  clearConstraints();
  m_segment_indices.clear();

  const unsigned int per_strand = segments + 1;
  const unsigned int n = static_cast<unsigned int>(roots.size()) * per_strand;
  const EigenVector3 dir = direction.normalized();

  m_positions.resize(n, 3);
  m_velocities.setZero(n, 3);

  for (unsigned int s = 0; s < roots.size(); ++s) {
    const unsigned int base = s * per_strand;
    for (unsigned int k = 0; k < per_strand; ++k) {
      m_positions.row(base + k) = (roots[s] + dir * (segment_length * k)).transpose();
    }

    m_constraints.push_back(new AttachmentConstraint(&m_stiffness_attachment, base, roots[s]));

    for (unsigned int k = 0; k + 1 < per_strand; ++k) {
      m_constraints.push_back(new SpringConstraint(&m_stiffness_stretch, base + k, base + k + 1,
                                                   segment_length));
      m_segment_indices.push_back(base + k);
      m_segment_indices.push_back(base + k + 1);
    }
    for (unsigned int k = 0; k + 2 < per_strand; ++k) {
      m_constraints.push_back(new SpringConstraint(&m_stiffness_bending, base + k, base + k + 2,
                                                   2.0 * segment_length));
    }
  }

  m_initial_positions = m_positions;
  m_time_accumulator = 0.0;
  RebuildSystemMatrix();
}

void Simulation::RebuildSystemMatrix() {
  const unsigned int n = NumParticles();
  std::vector<SparseMatrixTriplet> triplets;
  triplets.reserve(n + m_constraints.size() * 4);

  // M / h^2 (lumped, uniform mass).
  const ScalarType inertia = m_particle_mass / (m_h * m_h);
  for (unsigned int i = 0; i < n; ++i) triplets.emplace_back(i, i, inertia);

  // + sum_i w_i G_i^T G_i
  for (const Constraint* c : m_constraints) c->AddToLHS(&triplets);
  if (m_drag_constraint) m_drag_constraint->AddToLHS(&triplets);

  SparseMatrix lhs(n, n);
  lhs.setFromTriplets(triplets.begin(), triplets.end());  // Duplicates are summed.

  m_llt.compute(lhs);
  m_system_ready = (m_llt.info() == Eigen::Success);
  if (!m_system_ready) std::cerr << "Simulation: system matrix factorization failed" << std::endl;
}

void Simulation::BeginDrag(unsigned int particle, const EigenVector3& target) {
  if (particle >= NumParticles()) return;
  delete m_drag_constraint;
  m_drag_constraint = new AttachmentConstraint(&m_stiffness_drag, particle, target);
  RebuildSystemMatrix();
}

void Simulation::UpdateDrag(const EigenVector3& target) {
  if (m_drag_constraint) m_drag_constraint->SetFixedPoint(target);
}

void Simulation::EndDrag() {
  if (!m_drag_constraint) return;
  delete m_drag_constraint;
  m_drag_constraint = nullptr;
  RebuildSystemMatrix();
}

int Simulation::DraggedParticle() const {
  return m_drag_constraint ? static_cast<int>(m_drag_constraint->GetConstrainedVertexIndex()) : -1;
}

void Simulation::Reset() {
  EndDrag();
  m_positions = m_initial_positions;
  m_velocities.setZero(m_positions.rows(), 3);
  m_time_accumulator = 0.0;
}

void Simulation::Update(float dt) {
  if (!m_system_ready) return;
  // Clamp very long frames (e.g. window drag) so we don't spiral.
  m_time_accumulator += std::min(static_cast<double>(dt), m_h * m_max_substeps_per_frame);
  while (m_time_accumulator >= m_h) {
    Step();
    m_time_accumulator -= m_h;
  }
}

void Simulation::Step() {
  if (!m_system_ready) return;
  const unsigned int n = NumParticles();
  const ScalarType inertia = m_particle_mass / (m_h * m_h);

  // Inertial prediction and the constant part of the RHS: M/h^2 y + f_ext.
  const EigenMatrix3x y = m_positions + m_h * m_velocities;
  EigenMatrix3x b0 = inertia * y;
  b0.col(1).array() -= m_particle_mass * m_gravity_constant;

  EigenMatrix3x x = y;
  EigenMatrix3x b(n, 3);
  for (int it = 0; it < m_iterations_per_step; ++it) {
    // Local step: each constraint projects independently (serial for now; this
    // loop is what would move to a GPU kernel).
    b = b0;
    for (const Constraint* c : m_constraints) c->AddToRHS(x, &b);
    if (m_drag_constraint) m_drag_constraint->AddToRHS(x, &b);
    // Global step: same prefactored matrix for x, y and z.
    x = m_llt.solve(b);
  }

  m_velocities = (x - m_positions) / m_h;
  m_velocities *= (1.0 - m_damping_coefficient);
  m_positions = x;
}
