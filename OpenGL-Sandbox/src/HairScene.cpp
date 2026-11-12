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

HairScene::HairScene() : HairScene(DefaultRoots(), EigenVector3(1.0, 0.0, 0.0), 20, 0.05) {}

HairScene::HairScene(const std::vector<EigenVector3>& roots, const EigenVector3& direction,
                     unsigned int segments, ScalarType segment_length)
    : m_stiffness_stretch(1e6),
      m_stiffness_bending(5e2),
      m_segment_mass(1.0),
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
  solver->Finalize();
}

std::vector<SceneParam> HairScene::Params() {
  return {
      {"Stretch stiffness", &m_stiffness_stretch, 1e3f, 1e7f, true, SceneParam::kSystemMatrix},
      {"Bending stiffness", &m_stiffness_bending, 1.0f, 1e5f, true, SceneParam::kSystemMatrix},
      {"Segment mass", &m_segment_mass, 0.01f, 10.0f, true, SceneParam::kRebuildScene},
  };
}
