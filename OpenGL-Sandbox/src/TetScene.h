#pragma once
#include <vector>

#include "Scene.h"
#include "TetMesh.h"

// A tetrahedral volume, clamped (Dirichlet) at its minimum-x end, sagging under
// gravity like a cantilever. Masses are lumped from density * tet volume.
//
// The material is the PD continuum model of Bouaziz et al. 2014: one corotational
// strain constraint per tet (singular values of F clamped to
// [m_min_stretch, m_max_stretch]) plus one volume constraint per tet (det F
// towards 1). Both are weighted by rest volume, so the stiffness is a material
// property independent of the mesh resolution.
class TetScene : public Scene {
 public:
  // Default: procedural 1 x 0.2 x 0.2 m beam, m_resolution cells along its length.
  TetScene();
  // An arbitrary mesh (e.g. from LoadTetGen). Resolution has no effect then.
  explicit TetScene(TetMesh mesh);

  const char* Name() const override { return "Tet beam"; }
  void Build(PDSolver* solver) override;
  // Wireframe of the tet edges (toggle with m_show_wireframe).
  const std::vector<unsigned int>& LineIndices() const override;
  const std::vector<unsigned int>& SurfaceIndices() const override { return m_surface_indices; }
  std::vector<SceneParam> Params() override;

  const TetMesh& Mesh() const { return m_mesh; }

  ScalarType m_stiffness_strain;  // Weight per unit rest volume (Pa-like).
  ScalarType m_stiffness_volume;
  ScalarType m_min_stretch;       // Strain clamp range (RHS only).
  ScalarType m_max_stretch;
  ScalarType m_min_volume;        // Allowed det F range (RHS only).
  ScalarType m_max_volume;
  ScalarType m_density;     // kg / m^3.
  ScalarType m_resolution;  // Cells along the beam (procedural beam only).
  bool m_show_wireframe;

 private:
  void MakeBeam();

  bool m_procedural;
  TetMesh m_mesh;
  std::vector<unsigned int> m_edge_indices;
  std::vector<unsigned int> m_surface_indices;
};
