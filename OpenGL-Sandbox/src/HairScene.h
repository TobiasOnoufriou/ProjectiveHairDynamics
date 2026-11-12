#pragma once
#include <vector>

#include "Scene.h"

// Strand-based hair: each strand is a particle chain whose root is pinned
// exactly (Dirichlet), with stretch springs between neighbours (i, i+1) and
// weaker bending springs between every second particle (i, i+2).
//
// Masses are lumped per vertex: every segment carries m_segment_mass, split
// half to each end, so interior particles get one segment mass and the tip half.
class HairScene : public Scene {
 public:
  // Default: a 5x5 patch of horizontal strands that swing down under gravity.
  HairScene();
  // Strands rooted at `roots`, all starting straight along `direction`.
  HairScene(const std::vector<EigenVector3>& roots, const EigenVector3& direction,
            unsigned int segments, ScalarType segment_length);

  const char* Name() const override { return "Hair"; }
  void Build(PDSolver* solver) override;
  const std::vector<unsigned int>& LineIndices() const override { return m_segment_indices; }
  std::vector<SceneParam> Params() override;

  unsigned int ParticlesPerStrand() const { return m_segments + 1; }

  ScalarType m_stiffness_stretch;
  ScalarType m_stiffness_bending;
  ScalarType m_segment_mass;

 private:
  std::vector<EigenVector3> m_roots;
  EigenVector3 m_direction;
  unsigned int m_segments;
  ScalarType m_segment_length;
  std::vector<unsigned int> m_segment_indices;
};
