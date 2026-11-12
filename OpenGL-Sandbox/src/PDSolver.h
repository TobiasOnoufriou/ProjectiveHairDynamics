#pragma once
#include <memory>
#include <vector>

#include "Constraint.h"
#include "math_headers.h"

// Diagnostics of the most recent Step().
struct PDStepStats {
  double step_ms = 0.0;             // Wall-clock time of the step.
  // Total PD energy E(x_k) before iteration k and after the last one
  // (iterations + 1 values). Only filled when m_track_energy is set.
  std::vector<ScalarType> energy;
  ScalarType max_strain = 0.0;      // Max Constraint::Strain after the step.
};

// Generic Projective Dynamics solver (Bouaziz et al. 2014) on a set of
// vertices with per-vertex (lumped) masses and an arbitrary constraint list.
// Scene builders (e.g. HairScene) fill it via SetVertices / SetFixed /
// AddConstraint and then call Finalize().
//
// One step (h = fixed time step, M = diag(m_i)):
//   y  = x + h v + h^2 g                     (inertia + gravity, free vertices)
//   b0 = M/h^2 y
//   repeat iterations:
//     local:  b = b0 + sum_i w_i G_i^T p_i(x)
//     global: A_ff x_f = b_f - A_fc x_c      (prefactored)
//   v = (x_new - x_old) / h, then damped.
// A = M/h^2 + sum_i w_i G_i^T G_i. Fixed (Dirichlet) vertices x_c are eliminated
// from the system: only the free block A_ff is factored, and the coupling A_fc
// moves their known positions to the RHS. Moving a fixed vertex (e.g. a bone
// attachment) therefore only changes the RHS.
//
// Each iteration minimises g(x, p) = 1/(2h^2) ||x - y||_M^2 + sum_i w_i/2
// ||G_i x - p_i||^2 alternately over p and x, so the energy
// E(x) = min_p g(x, p) never increases across iterations.
class PDSolver {
 public:
  PDSolver();
  ~PDSolver();

  // --- Scene construction ----------------------------------------------------
  // Removes all vertices and constraints.
  void Clear();
  // positions: N x 3, masses: N lumped vertex masses. Marks every vertex free.
  void SetVertices(const EigenMatrix3x& positions, const VectorX& masses);
  // Pins a vertex exactly at its current position (Dirichlet). Changes the
  // system's structure, so RebuildSystemMatrix() (or Finalize) must follow.
  void SetFixed(unsigned int vertex, bool fixed = true);
  void AddConstraint(std::unique_ptr<Constraint> constraint);
  // Stores the current state as the reset pose and factors the system.
  void Finalize();

  // Moves a fixed vertex (kinematic target). RHS only, no refactorization.
  void SetFixedTarget(unsigned int vertex, const EigenVector3& target);

  // --- Time stepping ---------------------------------------------------------
  // Advances by frame time dt (seconds) using fixed sub-steps of m_h, so the
  // simulation speed does not depend on the frame rate.
  void Update(float dt);
  // A single fixed time step.
  void Step();
  // Puts every vertex back in the pose stored by Finalize().
  void Reset();
  // Must be called after changing any stiffness, mass, pinning or m_h.
  void RebuildSystemMatrix();

  // Mouse dragging: a temporary soft attachment pulling one free vertex towards
  // a target point. Begin/End refactor the system matrix; moving the target only
  // changes the right-hand side, so UpdateDrag is cheap.
  void BeginDrag(unsigned int vertex, const EigenVector3& target);
  void UpdateDrag(const EigenVector3& target);
  void EndDrag();
  bool IsDragging() const { return m_drag_constraint != nullptr; }
  int DraggedVertex() const;

  // --- Queries -----------------------------------------------------------------
  unsigned int NumVertices() const { return static_cast<unsigned int>(m_positions.rows()); }
  unsigned int NumFixed() const { return static_cast<unsigned int>(m_fixed_vertices.size()); }
  size_t NumConstraints() const { return m_constraints.size(); }
  bool IsFixed(unsigned int vertex) const { return m_free_index[vertex] < 0; }
  bool IsReady() const { return m_system_ready; }
  const EigenMatrix3x& GetPositions() const { return m_positions; }
  const EigenMatrix3x& GetVelocities() const { return m_velocities; }
  const VectorX& GetMasses() const { return m_masses; }
  const PDStepStats& GetStats() const { return m_stats; }

  // Sum of the constraint energies w_i/2 ||G_i x - p_i(x)||^2 (drag included).
  ScalarType ConstraintEnergy(const EigenMatrix3x& x) const;
  // Largest Constraint::Strain over all constraints at x.
  ScalarType MaxStrain(const EigenMatrix3x& x) const;

  // --- Tunable parameters ------------------------------------------------------
  // Changing m_h or m_stiffness_drag (while dragging) requires RebuildSystemMatrix().
  ScalarType m_h;                    // Fixed time step (s).
  ScalarType m_gravity_constant;     // Acceleration along -y (m/s^2).
  ScalarType m_stiffness_drag;
  ScalarType m_damping_coefficient;  // Fraction of velocity removed per step.
  int m_iterations_per_step;
  int m_max_substeps_per_frame;
  bool m_track_energy;               // Record PDStepStats::energy (costs one extra local pass).

 private:
  // Inertial part 1/(2h^2) sum_free m_i ||x_i - y_i||^2.
  ScalarType inertiaEnergy(const EigenMatrix3x& x, const EigenMatrix3x& y) const;
  // Local step over all constraints; returns the summed constraint energy.
  ScalarType localStep(const EigenMatrix3x& x, EigenMatrix3x* b) const;

  std::vector<std::unique_ptr<Constraint>> m_constraints;
  std::unique_ptr<AttachmentConstraint> m_drag_constraint;  // Null when nothing is dragged.

  EigenMatrix3x m_positions;          // N x 3, row i = vertex i.
  EigenMatrix3x m_velocities;         // N x 3.
  EigenMatrix3x m_initial_positions;
  VectorX m_masses;                   // N lumped masses.

  // Dirichlet bookkeeping: m_free_index[i] is i's row in the reduced system,
  // or -1 if vertex i is fixed.
  std::vector<int> m_free_index;
  std::vector<unsigned int> m_free_vertices;
  std::vector<unsigned int> m_fixed_vertices;

  Eigen::SimplicialLLT<SparseMatrix> m_llt;  // Factor of A_ff.
  SparseMatrix m_coupling;                   // A_fc (free x fixed).
  bool m_system_ready;
  double m_time_accumulator;
  PDStepStats m_stats;
};
