#pragma once
#include <vector>

#include "Scene.h"
#include "TetMesh.h"

// A tetrahedral volume, clamped (Dirichlet) at its minimum-x end, sagging under
// gravity like a cantilever. Masses are lumped from density * tet volume.
//
// Until the volumetric constraints of Phase 3 land, the material is a mass-spring
// network on the tet edges (Liu et al. 2013): every edge keeps its rest length.
// A tet mesh is a rigid framework, so this already holds its shape, but it is
// mesh-dependent and does not preserve volume; the corotational strain constraint
// (Bouaziz et al. 2014) replaces it.
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

  ScalarType m_stiffness_edge;
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
