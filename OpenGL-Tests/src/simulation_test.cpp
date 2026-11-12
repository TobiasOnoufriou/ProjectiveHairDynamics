// Copyright 2022 Tobias Onoufriou
// Behavioural tests for the Projective Dynamics hair solver (Simulation).
#include <vector>

#include "gtest/gtest.h"
#include "Simulation.h"

namespace {

constexpr unsigned int kSegments = 20;
constexpr ScalarType kSegmentLength = 0.05;
constexpr unsigned int kPerStrand = kSegments + 1;

// One strand rooted at the origin, pointing along `direction`.
void InitSingleStrand(Simulation* sim, const EigenVector3& direction) {
  std::vector<EigenVector3> roots = {EigenVector3(0.0, 0.0, 0.0)};
  sim->InitHairStrands(roots, direction, kSegments, kSegmentLength);
}

void StepN(Simulation* sim, int n) {
  for (int i = 0; i < n; ++i) sim->Step();
}

// Largest relative deviation of a stretch segment from its rest length.
ScalarType MaxStretchError(const Simulation& sim) {
  const EigenMatrix3x& x = sim.GetPositions();
  ScalarType worst = 0.0;
  for (unsigned int k = 0; k + 1 < kPerStrand; ++k) {
    const ScalarType len = (x.row(k + 1) - x.row(k)).norm();
    worst = std::max(worst, std::abs(len - kSegmentLength) / kSegmentLength);
  }
  return worst;
}

}  // namespace

TEST(SimulationTest, InitBuildsStraightStrand) {
  Simulation sim;
  InitSingleStrand(&sim, EigenVector3(1.0, 0.0, 0.0));

  ASSERT_EQ(sim.NumParticles(), kPerStrand);
  EXPECT_EQ(sim.GetSegmentIndices().size(), 2u * kSegments);
  EXPECT_NEAR(MaxStretchError(sim), 0.0, 1e-12);
  EXPECT_NEAR(sim.GetPositions()(kSegments, 0), kSegments * kSegmentLength, 1e-12);
}

// With no external force, the rest pose is already the energy minimum, so a
// step must not move anything.
TEST(SimulationTest, RestPoseIsStableWithoutGravity) {
  Simulation sim;
  sim.m_gravity_constant = 0.0;
  InitSingleStrand(&sim, EigenVector3(1.0, 0.0, 0.0));
  const EigenMatrix3x initial = sim.GetPositions();

  StepN(&sim, 30);

  EXPECT_LT((sim.GetPositions() - initial).cwiseAbs().maxCoeff(), 1e-9);
  EXPECT_LT(sim.GetVelocities().cwiseAbs().maxCoeff(), 1e-9);
}

// A horizontal strand swings down under gravity and settles hanging below its
// root, staying pinned at the root and close to its rest length.
TEST(SimulationTest, StrandHangsUnderGravity) {
  Simulation sim;
  InitSingleStrand(&sim, EigenVector3(1.0, 0.0, 0.0));

  StepN(&sim, 600);  // 10 s at 60 Hz.

  const EigenMatrix3x& x = sim.GetPositions();
  const ScalarType strand_length = kSegments * kSegmentLength;
  EXPECT_LT(x.row(0).norm(), 1e-2) << "root drifted from its attachment";
  EXPECT_LT(x(kSegments, 1), -0.9 * strand_length) << "tip is not hanging below the root";
  EXPECT_LT(std::abs(x(kSegments, 0)), 0.1 * strand_length) << "strand has not settled";
  EXPECT_LT(MaxStretchError(sim), 0.05);
  EXPECT_LT(sim.GetVelocities().cwiseAbs().maxCoeff(), 0.05) << "strand is still moving";
}

TEST(SimulationTest, DragPullsParticleTowardsTarget) {
  Simulation sim;
  sim.m_gravity_constant = 0.0;
  InitSingleStrand(&sim, EigenVector3(1.0, 0.0, 0.0));
  const unsigned int tip = kSegments;
  const EigenVector3 target = sim.GetPositions().row(tip).transpose() + EigenVector3(0.0, 0.0, 0.3);
  const ScalarType start_dist = (sim.GetPositions().row(tip).transpose() - target).norm();

  sim.BeginDrag(tip, target);
  ASSERT_TRUE(sim.IsDragging());
  EXPECT_EQ(sim.DraggedParticle(), static_cast<int>(tip));
  StepN(&sim, 120);
  const ScalarType end_dist = (sim.GetPositions().row(tip).transpose() - target).norm();
  sim.EndDrag();

  EXPECT_LT(end_dist, 0.5 * start_dist);
  EXPECT_FALSE(sim.IsDragging());
  EXPECT_EQ(sim.DraggedParticle(), -1);
}

TEST(SimulationTest, ResetRestoresInitialPose) {
  Simulation sim;
  InitSingleStrand(&sim, EigenVector3(1.0, 0.0, 0.0));
  const EigenMatrix3x initial = sim.GetPositions();

  StepN(&sim, 60);
  ASSERT_GT((sim.GetPositions() - initial).cwiseAbs().maxCoeff(), 1e-3);
  sim.Reset();

  EXPECT_EQ((sim.GetPositions() - initial).cwiseAbs().maxCoeff(), 0.0);
  EXPECT_EQ(sim.GetVelocities().cwiseAbs().maxCoeff(), 0.0);
}

// Update() advances in fixed sub-steps, so the result must not depend on how
// the frame time is split up.
TEST(SimulationTest, UpdateIsFrameRateIndependent) {
  Simulation a, b;
  InitSingleStrand(&a, EigenVector3(1.0, 0.0, 0.0));
  InitSingleStrand(&b, EigenVector3(1.0, 0.0, 0.0));

  for (int i = 0; i < 60; ++i) a.Update(1.0f / 60.0f);   // 60 fps for 1 s.
  for (int i = 0; i < 30; ++i) b.Update(1.0f / 30.0f);   // 30 fps for 1 s.

  EXPECT_LT((a.GetPositions() - b.GetPositions()).cwiseAbs().maxCoeff(), 1e-9);
}
