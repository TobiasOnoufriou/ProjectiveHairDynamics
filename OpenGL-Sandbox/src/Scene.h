#pragma once
#include <vector>

#include "PDSolver.h"
#include "math_headers.h"

// A tunable scene parameter, drawn generically by the UI (no ImGui here, so
// scenes stay testable headless).
struct SceneParam {
  enum Effect {
    kNone,          // Read every step (RHS only).
    kSystemMatrix,  // Changes the LHS: PDSolver::RebuildSystemMatrix().
    kRebuildScene,  // Changes geometry or masses: Scene::Build() again.
  };
  const char* name;
  ScalarType* value;
  float min, max;
  bool logarithmic;
  Effect effect;
};

// Builds a simulation into a PDSolver: vertices, masses, pinned vertices and
// constraints. Constraints point at stiffness values owned by the scene, so the
// scene must outlive the solver's current contents.
class Scene {
 public:
  virtual ~Scene() {}

  virtual const char* Name() const = 0;
  // Clears the solver and fills it with this scene, then finalizes it.
  virtual void Build(PDSolver* solver) = 0;
  // Pairs of vertex indices drawn as line segments.
  virtual const std::vector<unsigned int>& LineIndices() const = 0;
  // Boundary surface as index triples, counter-clockwise seen from outside, drawn
  // as a lit surface. Empty for line-only scenes such as hair.
  virtual const std::vector<unsigned int>& SurfaceIndices() const {
    static const std::vector<unsigned int> kNone;
    return kNone;
  }
  virtual std::vector<SceneParam> Params() { return {}; }
};
