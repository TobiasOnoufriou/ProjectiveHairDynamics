#include "TetScene.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

#include "TetConstraints.h"

TetScene::TetScene()
    : m_stiffness_strain(1e6),
      m_stiffness_volume(1e6),
      m_min_stretch(1.0),
      m_max_stretch(1.0),
      m_min_volume(1.0),
      m_max_volume(1.0),
      m_density(1000.0),
      m_resolution(10.0),
      m_ground_height(0.45),
      m_friction(0.5),
      m_stiffness_contact(1e5),
      m_show_wireframe(false),
      m_procedural(true) {
  MakeBeam();
}

TetScene::TetScene(TetMesh mesh)
    : m_stiffness_strain(1e6),
      m_stiffness_volume(1e6),
      m_min_stretch(1.0),
      m_max_stretch(1.0),
      m_min_volume(1.0),
      m_max_volume(1.0),
      m_density(1000.0),
      m_resolution(0.0),
      m_ground_height(0.45),
      m_friction(0.5),
      m_stiffness_contact(1e5),
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

  for (const TetMesh::Tet& tet : m_mesh.tets) {
    solver->AddConstraint(std::make_unique<TetStrainConstraint>(&m_stiffness_strain, tet, m_mesh.vertices,
                                                                &m_min_stretch, &m_max_stretch));
    solver->AddConstraint(std::make_unique<TetVolumeConstraint>(&m_stiffness_volume, tet, m_mesh.vertices,
                                                                &m_min_volume, &m_max_volume));
  }

  m_edge_indices.clear();
  for (const auto& e : m_mesh.Edges()) {
    m_edge_indices.push_back(e[0]);
    m_edge_indices.push_back(e[1]);
  }

  m_surface_indices.clear();
  for (const TetMesh::Triangle& tri : m_mesh.BoundaryTriangles()) {
    m_surface_indices.insert(m_surface_indices.end(), tri.begin(), tri.end());
  }

  m_ground = -1;
  if (m_procedural) {
    // The beam spans x in [-0.5, 0.5] at y in [0.8, 1.0]: the sphere's top is
    // 6 cm below its free half.
    m_ground = solver->AddCollider(Collider::Plane(EigenVector3(0.0, m_ground_height, 0.0), EigenVector3::UnitY()));
    solver->AddCollider(Collider::Sphere(EigenVector3(0.2, 0.62, 0.0), 0.12));
    // Interior vertices cannot reach a collider before the surface does.
    solver->SetContactVertices(m_surface_indices);
    solver->SetContactStiffness(&m_stiffness_contact);
    Update(solver, 0.0);
  }
  solver->Finalize();
}

void TetScene::Update(PDSolver* solver, double /*dt*/) {
  if (m_ground < 0) return;
  solver->GetCollider(m_ground).translation.y() = m_ground_height;
  for (int k = 0; k < static_cast<int>(solver->GetColliders().size()); ++k) {
    solver->GetCollider(k).friction = m_friction;
  }
}

const std::vector<unsigned int>& TetScene::LineIndices() const {
  static const std::vector<unsigned int> kNone;
  return m_show_wireframe ? m_edge_indices : kNone;
}

std::vector<SceneParam> TetScene::Params() {
  std::vector<SceneParam> params = {
      {"Strain stiffness", &m_stiffness_strain, 1e3f, 1e9f, true, SceneParam::kSystemMatrix},
      {"Volume stiffness", &m_stiffness_volume, 0.0f, 1e9f, true, SceneParam::kSystemMatrix},
      {"Min stretch", &m_min_stretch, 0.0f, 1.0f, false, SceneParam::kNone},
      {"Max stretch", &m_max_stretch, 1.0f, 3.0f, false, SceneParam::kNone},
      {"Density (kg/m^3)", &m_density, 10.0f, 1e4f, true, SceneParam::kRebuildScene},
  };
  if (m_procedural) {
    params.push_back({"Cells along beam", &m_resolution, 1.0f, 40.0f, false, SceneParam::kRebuildScene});
    params.push_back({"Ground height", &m_ground_height, -0.5f, 0.75f, false, SceneParam::kNone});
    params.push_back({"Friction", &m_friction, 0.0f, 2.0f, false, SceneParam::kNone});
    params.push_back({"Contact stiffness", &m_stiffness_contact, 1e3f, 1e7f, true, SceneParam::kSystemMatrix});
  }
  return params;
}
