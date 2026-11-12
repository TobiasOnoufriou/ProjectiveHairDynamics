#pragma once
#include <vector>
#include "math_headers.h"
#include "Constraint.h"

// Projective Dynamics hair simulation (CPU).
//
// Each strand is a chain of particles. The root particle is held by an
// attachment constraint, neighbouring particles by stretch springs, and every
// second particle by weaker bending springs.
//
// One step (h = fixed time step, M = mass matrix):
//   y  = x + h v                       (inertia)
//   b0 = M/h^2 y + f_ext
//   repeat iterations:
//     local:  b = b0 + sum_i w_i G_i^T p_i(x)
//     global: x = (M/h^2 + sum_i w_i G_i^T G_i)^-1 b   (prefactored)
//   v = (x_new - x_old) / h, then damped.
class Simulation {
 public:
  Simulation();
  ~Simulation();

  // Builds a grid of strands. roots are the root positions, every strand starts
  // straight along direction with `segments` segments of `segment_length`.
  void InitHairStrands(const std::vector<EigenVector3>& roots,
                       const EigenVector3& direction,
                       unsigned int segments, ScalarType segment_length);

  // Advances by frame time dt (seconds) using fixed sub-steps of m_h, so the
  // simulation speed does not depend on the frame rate.
  void Update(float dt);
  // A single fixed time step.
  void Step();
  // Puts every strand back in its initial pose.
  void Reset();
  // Must be called after changing any stiffness or the time step.
  void RebuildSystemMatrix();

  unsigned int NumParticles() const { return static_cast<unsigned int>(m_positions.rows()); }
  const EigenMatrix3x& GetPositions() const { return m_positions; }
  const EigenMatrix3x& GetVelocities() const { return m_velocities; }
  // Pairs of particle indices to draw as line segments.
  const std::vector<unsigned int>& GetSegmentIndices() const { return m_segment_indices; }

  // Tunable parameters. Changing a stiffness or m_h requires RebuildSystemMatrix().
  ScalarType m_h;                   // Fixed time step (s).
  ScalarType m_particle_mass;
  ScalarType m_gravity_constant;
  ScalarType m_stiffness_attachment;
  ScalarType m_stiffness_stretch;
  ScalarType m_stiffness_bending;
  ScalarType m_damping_coefficient;  // Fraction of velocity removed per step.
  int m_iterations_per_step;
  int m_max_substeps_per_frame;

 private:
  void clearConstraints();

  std::vector<Constraint*> m_constraints;

  EigenMatrix3x m_positions;         // N x 3, row i = particle i.
  EigenMatrix3x m_velocities;        // N x 3.
  EigenMatrix3x m_initial_positions;
  std::vector<unsigned int> m_segment_indices;

  Eigen::SimplicialLLT<SparseMatrix> m_llt;
  bool m_system_ready;
  double m_time_accumulator;
};
