#pragma once
#include <vector>

#include "Scene.h"
#include "TetMesh.h"

// A kinematic bone: a rigid set of solver vertices (its render box plus the
// tendon vertices attached to it), all fixed (Dirichlet). Posing the bone only
// moves their targets through PDSolver::SetFixedTarget, so animating bones
// changes the right-hand side and never the system matrix.
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
  bool m_animate;
  bool m_show_wireframe;

  static constexpr ScalarType kMaxAnimatedAngleDeg = 120.0;

 private:
  void MakeMuscle();
  // Appends an axis-aligned box (4 vertices per face, for flat shading) to
  // *positions and its triangles to m_surface_indices; returns its vertex ids.
  std::vector<unsigned int> AddBox(const EigenVector3& lo, const EigenVector3& hi,
                                   std::vector<EigenVector3>* positions);

  TetMesh m_mesh;
  std::vector<char> m_is_tendon;  // Per tet.
  std::vector<unsigned int> m_origin_vertices;     // Muscle vertices on the upper arm.
  std::vector<unsigned int> m_insertion_vertices;  // Muscle vertices on the forearm.
  KinematicBone m_upper_arm;
  KinematicBone m_forearm;
  double m_phase;  // Animation phase in [0, 1).
  std::vector<unsigned int> m_edge_indices;
  std::vector<unsigned int> m_surface_indices;
};
