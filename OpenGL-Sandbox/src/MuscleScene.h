#pragma once
#include <vector>

#include "Scene.h"
#include "TetMesh.h"

// A kinematic bone: the rigid set of tendon vertices attached to it, all fixed
// (Dirichlet), plus a capsule collider owned by the solver. Posing the bone
// only moves their targets through PDSolver::SetFixedTarget (and the collider's
// pose), so animating bones changes the right-hand side and never the system
// matrix.
struct KinematicBone {
  EigenVector3 pivot = EigenVector3::Zero();  // Joint centre; the bone rotates about +z through it.
  std::vector<unsigned int> vertices;         // Solver vertices that move with the bone.
  EigenMatrix3x rest;                         // Their rest positions (row k = vertices[k]).

  // Position of a rest point when the bone is rotated by `angle` radians.
  EigenVector3 Transform(const EigenVector3& rest_point, ScalarType angle) const;
  // Moves every vertex of the bone to its posed target.
  void Apply(PDSolver* solver, ScalarType angle) const;
};

// Biceps-style muscle between two bones (Phase 4 demo).
//
// A static upper-arm bone runs along x up to the elbow; the forearm hinges about
// the elbow (z axis) and flexes upwards. The muscle is a fusiform tet mesh lying
// above the upper arm: its origin tendon is attached to the upper arm and its
// insertion tendon to the forearm just past the elbow, so flexing the elbow
// shortens the muscle.
//
// Material (Bouaziz et al. 2014; Saito, Zhou, Kavan 2015), all weighted by rest
// volume:
//   - every tet: corotational strain (passive tissue) + volume preservation,
//   - tendon tets (both ends): a much stiffer strain constraint and no fibres,
//   - belly tets: an active fibre constraint along a fibre field that follows
//     the spindle shape, with target stretch 1 - activation * max contraction.
// Contact (Phase 6): both bones are capsule colliders (drawn by the app), so
// the belly rests on the upper arm instead of sinking into it and the flexing
// forearm pushes the muscle instead of passing through it. The forearm
// collider is posed with the bone and carries its angular velocity for
// friction. The passive belly is soft, so under gravity it slumps onto the
// upper arm and drapes around it.
// Activation and the elbow angle only move projection targets (RHS). With
// m_animate the elbow flexes and extends periodically and the activation
// follows the flexion, so the belly shortens and bulges.
class MuscleScene : public Scene {
 public:
  MuscleScene();

  const char* Name() const override { return "Muscle"; }
  void Build(PDSolver* solver) override;
  const std::vector<unsigned int>& LineIndices() const override;
  const std::vector<unsigned int>& SurfaceIndices() const override { return m_surface_indices; }
  std::vector<SceneParam> Params() override;
  void Update(PDSolver* solver, double dt) override;

  // The muscle's tet mesh; its vertices are solver vertices 0..NumVertices()-1.
  const TetMesh& Mesh() const { return m_mesh; }
  // Whether tet t is tendon (no fibre, stiff) rather than muscle belly.
  bool IsTendon(size_t t) const { return m_is_tendon[t]; }
  const KinematicBone& UpperArm() const { return m_upper_arm; }
  const KinematicBone& Forearm() const { return m_forearm; }
  ScalarType ElbowAngleRadians() const;

  ScalarType m_stiffness_strain;  // Passive belly tissue.
  ScalarType m_stiffness_tendon;  // Strain stiffness of the tendon tets.
  ScalarType m_stiffness_volume;
  ScalarType m_stiffness_fibre;
  ScalarType m_activation;        // alpha in [0, 1] (RHS only).
  ScalarType m_max_contraction;   // Fibre stretch at full activation is 1 - this (RHS only).
  ScalarType m_elbow_angle_deg;   // Forearm flexion (RHS only).
  ScalarType m_cycle_seconds;     // Period of the flex/extend animation.
  ScalarType m_density;           // kg / m^3.
  ScalarType m_resolution;        // Cells along the muscle.
  ScalarType m_friction;          // Coulomb coefficient muscle-bone (RHS only).
  ScalarType m_stiffness_contact; // Contact slot weight.
  bool m_animate;
  bool m_bone_contact;            // Enables the bone colliders (RHS only).
  bool m_interior_contact;        // Contact slots on every free vertex, not only the boundary (rebuild).
  bool m_show_wireframe;

  static constexpr ScalarType kMaxAnimatedAngleDeg = 120.0;

 private:
  void MakeMuscle();
  // Poses the forearm bone and its collider at the current elbow angle; dt > 0
  // sets the collider's angular velocity from the change since the last pose.
  void PoseForearm(PDSolver* solver, double dt);

  TetMesh m_mesh;
  std::vector<char> m_is_tendon;  // Per tet.
  std::vector<unsigned int> m_origin_vertices;     // Muscle vertices on the upper arm.
  std::vector<unsigned int> m_insertion_vertices;  // Muscle vertices on the forearm.
  KinematicBone m_upper_arm;
  KinematicBone m_forearm;
  int m_upper_arm_collider = -1;  // Solver collider indices.
  int m_forearm_collider = -1;
  EigenVector3 m_forearm_rest_centre = EigenVector3::Zero();
  ScalarType m_posed_angle = 0.0;  // Elbow angle (rad) the forearm collider was last posed at.
  double m_phase;  // Animation phase in [0, 1).
  std::vector<unsigned int> m_edge_indices;
  std::vector<unsigned int> m_surface_indices;
};
