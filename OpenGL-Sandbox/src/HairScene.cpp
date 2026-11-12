#include "HairScene.h"

#include <memory>
#include <vector>

namespace {

std::vector<EigenVector3> DefaultRoots() {
  std::vector<EigenVector3> roots;
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j) roots.emplace_back(-0.4 + 0.2 * i, 1.5, -0.4 + 0.2 * j);
  return roots;
}

}  // namespace

HairScene::HairScene() : HairScene(DefaultRoots(), EigenVector3(1.0, 0.0, 0.0), 20, 0.05) {
  m_obstacles = true;
}

HairScene::HairScene(const std::vector<EigenVector3>& roots, const EigenVector3& direction,
                     unsigned int segments, ScalarType segment_length)
    : m_stiffness_stretch(1e6),
      m_stiffness_bending(5e2),
      m_segment_mass(1.0),
      m_ground_height(0.7),
      m_friction(0.5),
      m_stiffness_contact(1e5),
      m_obstacles(false),
      m_roots(roots),
      m_direction(direction.normalized()),
      m_segments(segments),
      m_segment_length(segment_length) {}

void HairScene::Build(PDSolver* solver) {
  solver->Clear();
  m_segment_indices.clear();

  const unsigned int per_strand = m_segments + 1;
  const unsigned int n = static_cast<unsigned int>(m_roots.size()) * per_strand;

  EigenMatrix3x positions(n, 3);
  VectorX masses = VectorX::Zero(n);
  for (unsigned int s = 0; s < m_roots.size(); ++s) {
    const unsigned int base = s * per_strand;
    for (unsigned int k = 0; k < per_strand; ++k) {
      positions.row(base + k) = (m_roots[s] + m_direction * (m_segment_length * k)).transpose();
    }
    // Lumped mass: half of each segment to each of its two particles.
    for (unsigned int k = 0; k + 1 < per_strand; ++k) {
      masses[base + k] += 0.5 * m_segment_mass;
      masses[base + k + 1] += 0.5 * m_segment_mass;
    }
  }
  solver->SetVertices(positions, masses);

  for (unsigned int s = 0; s < m_roots.size(); ++s) {
    const unsigned int base = s * per_strand;
    solver->SetFixed(base);  // Root.
    for (unsigned int k = 0; k + 1 < per_strand; ++k) {
      solver->AddConstraint(std::make_unique<SpringConstraint>(&m_stiffness_stretch, base + k,
                                                               base + k + 1, m_segment_length));
      m_segment_indices.push_back(base + k);
      m_segment_indices.push_back(base + k + 1);
    }
    for (unsigned int k = 0; k + 2 < per_strand; ++k) {
      solver->AddConstraint(std::make_unique<SpringConstraint>(&m_stiffness_bending, base + k,
                                                               base + k + 2, 2.0 * m_segment_length));
    }
  }
  m_ground = -1;
  if (m_obstacles) {
    // Strands hang 1 m below their roots (y = 1.5), so the tips reach the
    // ground; the sphere sits under the strands rooted at x = 0.2 and 0.4.
    m_ground = solver->AddCollider(Collider::Plane(EigenVector3(0.0, m_ground_height, 0.0), EigenVector3::UnitY()));
    solver->AddCollider(Collider::Sphere(EigenVector3(0.3, 1.05, 0.0), 0.2));
    solver->SetContactStiffness(&m_stiffness_contact);
    Update(solver, 0.0);
  }
  solver->Finalize();
}

void HairScene::Update(PDSolver* solver, double /*dt*/) {
  if (m_ground < 0) return;
  solver->GetCollider(m_ground).translation.y() = m_ground_height;
  for (int k = 0; k < static_cast<int>(solver->GetColliders().size()); ++k) {
    solver->GetCollider(k).friction = m_friction;
  }
}

std::vector<SceneParam> HairScene::Params() {
  std::vector<SceneParam> params = {
      {"Stretch stiffness", &m_stiffness_stretch, 1e3f, 1e7f, true, SceneParam::kSystemMatrix},
      {"Bending stiffness", &m_stiffness_bending, 1.0f, 1e5f, true, SceneParam::kSystemMatrix},
      {"Segment mass", &m_segment_mass, 0.01f, 10.0f, true, SceneParam::kRebuildScene},
  };
  if (m_obstacles) {
    params.push_back({"Ground height", &m_ground_height, -0.5f, 1.4f, false, SceneParam::kNone});
    params.push_back({"Friction", &m_friction, 0.0f, 2.0f, false, SceneParam::kNone});
    params.push_back({"Contact stiffness", &m_stiffness_contact, 1e3f, 1e7f, true, SceneParam::kSystemMatrix});
  }
  return params;
}
