#include "TetScene.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

TetScene::TetScene()
    : m_stiffness_edge(1e6),
      m_density(1000.0),
      m_resolution(10.0),
      m_show_wireframe(false),
      m_procedural(true) {
  MakeBeam();
}

TetScene::TetScene(TetMesh mesh)
    : m_stiffness_edge(1e6),
      m_density(1000.0),
      m_resolution(0.0),
      m_show_wireframe(false),
      m_procedural(false),
      m_mesh(std::move(mesh)) {}

void TetScene::MakeBeam() {
  // Square cross-section with cells of roughly equal size along all axes.
  const unsigned int nx = static_cast<unsigned int>(std::max(1.0, std::round(m_resolution)));
  const unsigned int ny = std::max(1u, nx / 5);
  const EigenVector3 size(1.0, 0.2, 0.2);
  m_mesh = MakeTetBeam(nx, ny, ny, size, EigenVector3(-0.5, 0.8, -0.1));
}

void TetScene::Build(PDSolver* solver) {
  if (m_procedural) MakeBeam();
  solver->Clear();
  solver->SetVertices(m_mesh.vertices, m_mesh.LumpedMasses(m_density));

  // Clamp the minimum-x end.
  const ScalarType min_x = m_mesh.vertices.col(0).minCoeff();
  const ScalarType extent = m_mesh.vertices.col(0).maxCoeff() - min_x;
  for (unsigned int i = 0; i < m_mesh.NumVertices(); ++i) {
    if (m_mesh.vertices(i, 0) <= min_x + 1e-6 * extent) solver->SetFixed(i);
  }

  m_edge_indices.clear();
  for (const auto& e : m_mesh.Edges()) {
    const ScalarType rest = (m_mesh.vertices.row(e[0]) - m_mesh.vertices.row(e[1])).norm();
    solver->AddConstraint(std::make_unique<SpringConstraint>(&m_stiffness_edge, e[0], e[1], rest));
    m_edge_indices.push_back(e[0]);
    m_edge_indices.push_back(e[1]);
  }

  m_surface_indices.clear();
  for (const TetMesh::Triangle& tri : m_mesh.BoundaryTriangles()) {
    m_surface_indices.insert(m_surface_indices.end(), tri.begin(), tri.end());
  }
  solver->Finalize();
}

const std::vector<unsigned int>& TetScene::LineIndices() const {
  static const std::vector<unsigned int> kNone;
  return m_show_wireframe ? m_edge_indices : kNone;
}

std::vector<SceneParam> TetScene::Params() {
  std::vector<SceneParam> params = {
      {"Edge stiffness", &m_stiffness_edge, 1e3f, 1e8f, true, SceneParam::kSystemMatrix},
      {"Density (kg/m^3)", &m_density, 10.0f, 1e4f, true, SceneParam::kRebuildScene},
  };
  if (m_procedural) params.push_back({"Cells along beam", &m_resolution, 1.0f, 40.0f, false, SceneParam::kRebuildScene});
  return params;
}
