// Copyright 2022 Tobias Onoufriou
// Behavioural tests for the hair scene running on the PD solver.
#include <algorithm>
#include <vector>

#include "gtest/gtest.h"
#include "HairScene.h"
#include "PDSolver.h"

namespace {

constexpr unsigned int kSegments = 20;
constexpr ScalarType kSegmentLength = 0.05;
constexpr unsigned int kPerStrand = kSegments + 1;

// One strand rooted at the origin, pointing along `direction`.
HairScene SingleStrand(const EigenVector3& direction) {
  return HairScene({EigenVector3(0.0, 0.0, 0.0)}, direction, kSegments, kSegmentLength);
}

void StepN(PDSolver* sim, int n) {
  for (int i = 0; i < n; ++i) sim->Step();
}

// Largest relative deviation of a stretch segment from its rest length.
ScalarType MaxStretchError(const PDSolver& sim) {
  const EigenMatrix3x& x = sim.GetPositions();
  ScalarType worst = 0.0;
  for (unsigned int k = 0; k + 1 < kPerStrand; ++k) {
    const ScalarType len = (x.row(k + 1) - x.row(k)).norm();
    worst = std::max(worst, std::abs(len - kSegmentLength) / kSegmentLength);
  }
  return worst;
}

}  // namespace

TEST(HairSceneTest, BuildsStraightStrand) {
  HairScene scene = SingleStrand(EigenVector3(1.0, 0.0, 0.0));
  PDSolver sim;
  scene.Build(&sim);

  ASSERT_EQ(sim.NumVertices(), kPerStrand);
  EXPECT_TRUE(sim.IsReady());
  EXPECT_EQ(scene.LineIndices().size(), 2u * kSegments);
  EXPECT_NEAR(MaxStretchError(sim), 0.0, 1e-12);
  EXPECT_NEAR(sim.GetPositions()(kSegments, 0), kSegments * kSegmentLength, 1e-12);
}

// Lumped mass: interior particles get one segment mass, the ends half.
TEST(HairSceneTest, LumpedMasses) {
  HairScene scene = SingleStrand(EigenVector3(1.0, 0.0, 0.0));
  scene.m_segment_mass = 2.0;
  PDSolver sim;
  scene.Build(&sim);

  const VectorX& m = sim.GetMasses();
  EXPECT_DOUBLE_EQ(m[0], 1.0);
  EXPECT_DOUBLE_EQ(m[kSegments / 2], 2.0);
  EXPECT_DOUBLE_EQ(m[kSegments], 1.0);
  EXPECT_NEAR(m.sum(), 2.0 * kSegments, 1e-12);
}

TEST(HairSceneTest, RootIsPinnedExactly) {
  HairScene scene = SingleStrand(EigenVector3(1.0, 0.0, 0.0));
  PDSolver sim;
  scene.Build(&sim);
  EXPECT_TRUE(sim.IsFixed(0));
  EXPECT_FALSE(sim.IsFixed(1));

  StepN(&sim, 120);
  EXPECT_EQ(sim.GetPositions().row(0).norm(), 0.0);
}

// With no external force, the rest pose is already the energy minimum, so a
// step must not move anything.
TEST(HairSceneTest, RestPoseIsStableWithoutGravity) {
  HairScene scene = SingleStrand(EigenVector3(1.0, 0.0, 0.0));
  PDSolver sim;
  sim.m_gravity_constant = 0.0;
  scene.Build(&sim);
  const EigenMatrix3x initial = sim.GetPositions();

  StepN(&sim, 30);

  EXPECT_LT((sim.GetPositions() - initial).cwiseAbs().maxCoeff(), 1e-9);
  EXPECT_LT(sim.GetVelocities().cwiseAbs().maxCoeff(), 1e-9);
}

// A horizontal strand swings down under gravity and settles hanging below its
// root, close to its rest length.
TEST(HairSceneTest, StrandHangsUnderGravity) {
  HairScene scene = SingleStrand(EigenVector3(1.0, 0.0, 0.0));
  PDSolver sim;
  scene.Build(&sim);

  StepN(&sim, 600);  // 10 s at 60 Hz.

  const EigenMatrix3x& x = sim.GetPositions();
  const ScalarType strand_length = kSegments * kSegmentLength;
  EXPECT_LT(x(kSegments, 1), -0.9 * strand_length) << "tip is not hanging below the root";
  EXPECT_LT(std::abs(x(kSegments, 0)), 0.1 * strand_length) << "strand has not settled";
  EXPECT_LT(MaxStretchError(sim), 0.05);
  EXPECT_LT(sim.GetVelocities().cwiseAbs().maxCoeff(), 0.05) << "strand is still moving";
}

TEST(HairSceneTest, DragPullsParticleTowardsTarget) {
  HairScene scene = SingleStrand(EigenVector3(1.0, 0.0, 0.0));
  PDSolver sim;
  sim.m_gravity_constant = 0.0;
  scene.Build(&sim);
  const unsigned int tip = kSegments;
  const EigenVector3 target = sim.GetPositions().row(tip).transpose() + EigenVector3(0.0, 0.0, 0.3);
  const ScalarType start_dist = (sim.GetPositions().row(tip).transpose() - target).norm();

  sim.BeginDrag(tip, target);
  ASSERT_TRUE(sim.IsDragging());
  EXPECT_EQ(sim.DraggedVertex(), static_cast<int>(tip));
  StepN(&sim, 120);
  const ScalarType end_dist = (sim.GetPositions().row(tip).transpose() - target).norm();
  sim.EndDrag();

  EXPECT_LT(end_dist, 0.5 * start_dist);
  EXPECT_FALSE(sim.IsDragging());
  EXPECT_EQ(sim.DraggedVertex(), -1);
}

// The pinned root cannot be dragged.
TEST(HairSceneTest, DragIgnoresFixedVertex) {
  HairScene scene = SingleStrand(EigenVector3(1.0, 0.0, 0.0));
  PDSolver sim;
  scene.Build(&sim);
  sim.BeginDrag(0, EigenVector3(0.0, 1.0, 0.0));
  EXPECT_FALSE(sim.IsDragging());
}

TEST(HairSceneTest, ResetRestoresInitialPose) {
  HairScene scene = SingleStrand(EigenVector3(1.0, 0.0, 0.0));
  PDSolver sim;
  scene.Build(&sim);
  const EigenMatrix3x initial = sim.GetPositions();

  StepN(&sim, 60);
  ASSERT_GT((sim.GetPositions() - initial).cwiseAbs().maxCoeff(), 1e-3);
  sim.Reset();

  EXPECT_EQ((sim.GetPositions() - initial).cwiseAbs().maxCoeff(), 0.0);
  EXPECT_EQ(sim.GetVelocities().cwiseAbs().maxCoeff(), 0.0);
}

// Update() advances in fixed sub-steps, so the result must not depend on how
// the frame time is split up.
TEST(HairSceneTest, UpdateIsFrameRateIndependent) {
  HairScene scene = SingleStrand(EigenVector3(1.0, 0.0, 0.0));
  PDSolver a, b;
  scene.Build(&a);
  HairScene scene_b = SingleStrand(EigenVector3(1.0, 0.0, 0.0));
  scene_b.Build(&b);
  // Wall-clock dependent; must not kick in (e.g. in slow Debug builds).
  a.m_drop_backlog_when_slow = b.m_drop_backlog_when_slow = false;

  for (int i = 0; i < 60; ++i) a.Update(1.0f / 60.0f);   // 60 fps for 1 s.
  for (int i = 0; i < 30; ++i) b.Update(1.0f / 30.0f);   // 30 fps for 1 s.

  EXPECT_LT((a.GetPositions() - b.GetPositions()).cwiseAbs().maxCoeff(), 1e-9);
}
