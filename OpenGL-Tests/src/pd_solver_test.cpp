// Copyright 2022 Tobias Onoufriou
// Tests for the generic Projective Dynamics solver: energy behaviour of the
// local/global iteration, per-vertex masses and Dirichlet (fixed) vertices.
#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include "gtest/gtest.h"
#include "HairScene.h"
#include "PDSolver.h"

namespace {

// Relative slack for comparing energies that should be non-increasing.
bool NonIncreasing(const std::vector<ScalarType>& e, ScalarType* worst_rise) {
  *worst_rise = 0.0;
  for (size_t k = 1; k < e.size(); ++k) {
    *worst_rise = std::max(*worst_rise, (e[k] - e[k - 1]) / std::max(e[k - 1], 1e-12));
  }
  return *worst_rise <= 1e-9;
}

// A dumbbell: two free particles with different masses joined by one spring.
void BuildDumbbell(PDSolver* sim, const ScalarType* stiffness, ScalarType m0, ScalarType m1) {
  EigenMatrix3x x(2, 3);
  x << 0.0, 0.0, 0.0,
       1.0, 0.0, 0.0;
  VectorX m(2);
  m << m0, m1;
  sim->SetVertices(x, m);
  sim->AddConstraint(std::make_unique<SpringConstraint>(stiffness, 0, 1, 1.0));
  sim->Finalize();
}

}  // namespace

// The rest state has zero PD energy: every constraint is satisfied, and
// without gravity y = x so the inertial term vanishes too.
TEST(PDSolverTest, RestStateHasZeroEnergy) {
  HairScene scene({EigenVector3(0, 0, 0), EigenVector3(0.2, 0, 0)}, EigenVector3(0, 0, 1), 10, 0.05);
  PDSolver sim;
  sim.m_gravity_constant = 0.0;
  sim.m_track_energy = true;
  scene.Build(&sim);

  // Zero up to round-off in the rest lengths.
  EXPECT_NEAR(sim.ConstraintEnergy(sim.GetPositions()), 0.0, 1e-20);
  EXPECT_NEAR(sim.MaxStrain(sim.GetPositions()), 0.0, 1e-12);

  sim.Step();
  ASSERT_EQ(sim.GetStats().energy.size(), static_cast<size_t>(sim.m_iterations_per_step + 1));
  for (ScalarType e : sim.GetStats().energy) EXPECT_NEAR(e, 0.0, 1e-18);
}

// Local and global steps each minimise the same objective, so E(x_k) must
// never increase within a step. Checked on a swinging hair patch (exact roots,
// gravity) and while dragging, over many steps and iteration counts.
TEST(PDSolverTest, EnergyNeverIncreasesAcrossIterations) {
  HairScene scene;  // Default 5x5 patch.
  PDSolver sim;
  sim.m_track_energy = true;
  sim.m_iterations_per_step = 25;
  scene.Build(&sim);

  for (int step = 0; step < 120; ++step) {
    if (step == 40) sim.BeginDrag(30, EigenVector3(0.5, 2.0, 0.3));
    if (step == 80) sim.EndDrag();
    sim.Step();
    const std::vector<ScalarType>& e = sim.GetStats().energy;
    ASSERT_EQ(e.size(), 26u);
    ScalarType rise = 0.0;
    ASSERT_TRUE(NonIncreasing(e, &rise)) << "step " << step << ": relative rise " << rise;
  }
  // The first step from a perturbed state actually reduces energy.
  sim.Reset();
  sim.Step();
  EXPECT_LT(sim.GetStats().energy.back(), sim.GetStats().energy.front());
}

TEST(PDSolverTest, StatsAreRecorded) {
  HairScene scene;
  PDSolver sim;
  scene.Build(&sim);
  sim.Step();
  EXPECT_GE(sim.GetStats().step_ms, 0.0);
  EXPECT_TRUE(sim.GetStats().energy.empty());  // Tracking is off by default.
  EXPECT_GT(sim.GetStats().max_strain, 0.0);
}

// Spring forces are internal, so the mass-weighted centre of mass must move
// only with its initial momentum. This fails if the solver used a uniform mass.
TEST(PDSolverTest, PerVertexMassConservesMomentum) {
  ScalarType stiffness = 1e4;
  PDSolver sim;
  sim.m_gravity_constant = 0.0;
  sim.m_damping_coefficient = 0.0;
  BuildDumbbell(&sim, &stiffness, 1.0, 4.0);

  // Stretch the spring and let it oscillate.
  sim.SetVertices((EigenMatrix3x(2, 3) << 0.0, 0.0, 0.0, 1.5, 0.0, 0.0).finished(), sim.GetMasses());
  sim.RebuildSystemMatrix();
  const VectorX m = sim.GetMasses();
  const EigenVector3 com0 = (sim.GetPositions().transpose() * m) / m.sum();
  for (int i = 0; i < 100; ++i) {
    sim.Step();
    const EigenVector3 com = (sim.GetPositions().transpose() * m) / m.sum();
    ASSERT_LT((com - com0).norm(), 1e-9) << "step " << i;
  }
  // The light particle moves four times as far as the heavy one.
  const ScalarType d0 = std::abs(sim.GetPositions()(0, 0) - 0.0);
  const ScalarType d1 = std::abs(sim.GetPositions()(1, 0) - 1.5);
  EXPECT_NEAR(d0, 4.0 * d1, 1e-6);
}

// A fixed vertex stays exactly where it is put, and moving it (as an animated
// bone would) is an RHS-only change that the free vertices follow.
TEST(PDSolverTest, DirichletTargetMovesWithoutRefactor) {
  HairScene scene({EigenVector3(0, 0, 0)}, EigenVector3(0, -1, 0), 10, 0.05);
  PDSolver sim;
  scene.Build(&sim);
  ASSERT_EQ(sim.NumFixed(), 1u);

  // Slide the root over one second, then let the strand settle under gravity.
  const EigenVector3 target(0.3, 0.1, -0.2);
  for (int i = 1; i <= 60; ++i) {
    sim.SetFixedTarget(0, target * (i / 60.0));
    sim.Step();
    ASSERT_EQ((sim.GetPositions().row(0).transpose() - target * (i / 60.0)).norm(), 0.0);
  }
  for (int i = 0; i < 600; ++i) sim.Step();

  EXPECT_EQ((sim.GetPositions().row(0).transpose() - target).norm(), 0.0);
  // The strand followed its root and hangs straight below it again.
  EXPECT_LT((sim.GetPositions().row(10).transpose() - (target + EigenVector3(0, -0.5, 0))).norm(),
            0.02);
}

// SetFixedTarget on a free vertex is ignored (it would otherwise teleport it).
TEST(PDSolverTest, SetFixedTargetIgnoresFreeVertex) {
  ScalarType stiffness = 1e4;
  PDSolver sim;
  BuildDumbbell(&sim, &stiffness, 1.0, 1.0);
  sim.SetFixedTarget(1, EigenVector3(5, 5, 5));
  EXPECT_EQ(sim.GetPositions()(1, 0), 1.0);
}

// With every vertex pinned there is nothing to solve, but stepping is safe.
TEST(PDSolverTest, AllFixedIsStable) {
  ScalarType stiffness = 1e4;
  PDSolver sim;
  BuildDumbbell(&sim, &stiffness, 1.0, 1.0);
  sim.SetFixed(0);
  sim.SetFixed(1);
  sim.RebuildSystemMatrix();
  ASSERT_TRUE(sim.IsReady());
  const EigenMatrix3x x0 = sim.GetPositions();
  sim.Step();
  EXPECT_EQ((sim.GetPositions() - x0).norm(), 0.0);
}
