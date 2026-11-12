// Copyright 2022 Tobias Onoufriou
// Tests for Phase 6 contact: collider geometry, the fixed-weight contact slots
// of PDSolver (projection out of colliders, RHS-only posing), Coulomb friction
// and contact in the beam and muscle scenes.
#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "Collider.h"
#include "HairScene.h"
#include "MuscleScene.h"
#include "PDSolver.h"
#include "TetScene.h"

namespace {

constexpr double kPi = 3.14159265358979323846;

// A few points in and around the unit cube, for checking distance functions.
std::vector<EigenVector3> SamplePoints() {
  std::vector<EigenVector3> points;
  for (int i = 0; i < 7; ++i) {
    for (int j = 0; j < 7; ++j) {
      for (int k = 0; k < 7; ++k) points.emplace_back(-0.9 + 0.31 * i, -0.8 + 0.27 * j, -0.85 + 0.29 * k);
    }
  }
  return points;
}

std::vector<Collider> AllShapes() {
  const EigenMatrix3 tilt =
      Eigen::AngleAxis<ScalarType>(0.4, EigenVector3(1.0, 2.0, -0.5).normalized()).toRotationMatrix();
  return {
      Collider::Plane(EigenVector3(0.1, -0.2, 0.0), EigenVector3(0.3, 1.0, 0.2)),
      Collider::Sphere(EigenVector3(0.1, 0.0, -0.1), 0.5),
      Collider::Capsule(EigenVector3(-0.4, -0.1, 0.0), EigenVector3(0.3, 0.4, 0.2), 0.25),
      Collider::Box(EigenVector3(0.05, 0.1, 0.0), EigenVector3(0.4, 0.2, 0.3), tilt),
  };
}

// One free particle of mass m with no constraints.
void BuildParticle(PDSolver* sim, const EigenVector3& x, ScalarType m) {
  sim->Clear();
  sim->SetVertices(EigenMatrix3x(x.transpose()), VectorX::Constant(1, m));
}

}  // namespace

// --- Collider geometry -----------------------------------------------------------

TEST(ColliderTest, SignedDistanceOfSimpleShapes) {
  EigenVector3 n;
  const Collider ground = Collider::Plane(EigenVector3(0, 1, 0), EigenVector3::UnitY());
  EXPECT_NEAR(ground.SignedDistance(EigenVector3(3, 0.25, -2), &n), -0.75, 1e-12);
  EXPECT_LT((n - EigenVector3::UnitY()).norm(), 1e-12);

  const Collider ball = Collider::Sphere(EigenVector3(1, 0, 0), 0.5);
  EXPECT_NEAR(ball.SignedDistance(EigenVector3(1, 2, 0), &n), 1.5, 1e-12);
  EXPECT_NEAR(ball.SignedDistance(EigenVector3(1.2, 0, 0), &n), -0.3, 1e-12);
  EXPECT_LT((n - EigenVector3::UnitX()).norm(), 1e-12);

  // Capsule along y from (0,0,0) to (0,2,0): the cylinder and the end caps.
  const Collider capsule = Collider::Capsule(EigenVector3::Zero(), EigenVector3(0, 2, 0), 0.5);
  EXPECT_NEAR(capsule.SignedDistance(EigenVector3(0.2, 1.5, 0), &n), -0.3, 1e-12);
  EXPECT_LT((n - EigenVector3::UnitX()).norm(), 1e-12);
  EXPECT_NEAR(capsule.SignedDistance(EigenVector3(0, 3, 0), &n), 0.5, 1e-12);
  EXPECT_LT((n - EigenVector3::UnitY()).norm(), 1e-12);

  // Box: inside, the nearest face wins; outside, distance to the clamped point.
  const Collider box = Collider::Box(EigenVector3::Zero(), EigenVector3(1, 0.5, 2));
  EXPECT_NEAR(box.SignedDistance(EigenVector3(0.2, -0.4, 0.0), &n), -0.1, 1e-12);
  EXPECT_LT((n + EigenVector3::UnitY()).norm(), 1e-12);
  EXPECT_NEAR(box.SignedDistance(EigenVector3(2, 1.5, 0), &n), std::sqrt(2.0), 1e-12);
}

// For every shape and a cloud of points, p - d n lies on the surface and a
// point inside projects onto it, while a point outside stays where it is.
TEST(ColliderTest, ProjectionLandsOnSurface) {
  for (const Collider& c : AllShapes()) {
    const std::vector<Collider> list = {c};
    int inside = 0;
    for (const EigenVector3& p : SamplePoints()) {
      EigenVector3 n;
      const ScalarType d = c.SignedDistance(p, &n);
      EXPECT_NEAR(n.norm(), 1.0, 1e-12);
      EXPECT_NEAR(c.SignedDistance(p - d * n, nullptr), 0.0, 1e-9);
      EigenVector3 q = p;
      const ScalarType depth = ProjectOutOfColliders(list, &q);
      if (d < 0.0) {
        ++inside;
        EXPECT_NEAR(depth, -d, 1e-12);
        EXPECT_LT((q - (p - d * n)).norm(), 1e-12);
      } else {
        EXPECT_EQ(depth, 0.0);
        EXPECT_EQ(q, p);
      }
    }
    EXPECT_GT(inside, 0) << "shape " << static_cast<int>(c.shape);
  }
}

// The projection is the closest outside point: no sampled surface point is
// closer to an interior point than its projection.
TEST(ColliderTest, ProjectionIsClosestSurfacePoint) {
  for (const Collider& c : AllShapes()) {
    if (c.shape == Collider::Shape::kPlane) continue;
    std::vector<EigenVector3> surface;
    std::vector<unsigned int> tris;
    TriangulateCollider(c, 2.0, &surface, &tris);
    for (const EigenVector3& p : SamplePoints()) {
      if (c.SignedDistance(p, nullptr) >= 0.0) continue;
      EigenVector3 q = p;
      ProjectOutOfColliders({c}, &q);
      for (const EigenVector3& s : surface) EXPECT_LE((q - p).norm(), (s - p).norm() + 1e-9);
    }
  }
}

TEST(ColliderTest, OverlappingCollidersResolveBoth) {
  // A point in the corner where a box sits on the ground ends up outside both.
  const std::vector<Collider> list = {Collider::Plane(EigenVector3::Zero(), EigenVector3::UnitY()),
                                      Collider::Box(EigenVector3(0, 0.25, 0), EigenVector3(0.5, 0.5, 0.5))};
  EigenVector3 p(0.45, -0.05, 0.0);
  ProjectOutOfColliders(list, &p);
  for (const Collider& c : list) EXPECT_GE(c.SignedDistance(p, nullptr), -1e-12);
}

TEST(ColliderTest, TriangulationIsOnSurfaceAndFacesOutwards) {
  for (const Collider& c : AllShapes()) {
    std::vector<EigenVector3> points;
    std::vector<unsigned int> tris;
    TriangulateCollider(c, 2.0, &points, &tris);
    ASSERT_FALSE(tris.empty());
    ASSERT_EQ(tris.size() % 3, 0u);
    for (const EigenVector3& p : points) EXPECT_NEAR(c.SignedDistance(p, nullptr), 0.0, 1e-9);
    for (size_t t = 0; t < tris.size(); t += 3) {
      const EigenVector3 a = points[tris[t]], b = points[tris[t + 1]], d = points[tris[t + 2]];
      const EigenVector3 area = (b - a).cross(d - a);
      if (area.norm() < 1e-12) continue;  // Pole triangles.
      EigenVector3 n;
      c.SignedDistance((a + b + d) / 3.0, &n);
      EXPECT_GT(area.dot(n), 0.0) << "shape " << static_cast<int>(c.shape) << " triangle " << t / 3 << " area " << area.norm();
    }
  }
}

TEST(ColliderTest, PointVelocityOfRotatingCollider) {
  Collider c = Collider::Box(EigenVector3(1, 0, 0), EigenVector3(0.5, 0.1, 0.1));
  c.linear_velocity = EigenVector3(0, 2, 0);
  c.angular_velocity = EigenVector3(0, 0, 3);
  // v + w x (p - centre) = (0, 2, 0) + 3 z x (0.5, 0, 0) = (0, 3.5, 0).
  EXPECT_LT((c.PointVelocity(EigenVector3(1.5, 0, 0)) - EigenVector3(0, 3.5, 0)).norm(), 1e-12);
}

// --- Contact slots in the solver --------------------------------------------------

TEST(ContactTest, SlotsOnlyForFreeContactVertices) {
  PDSolver sim;
  EigenMatrix3x x(4, 3);
  x << 0, 1, 0,  1, 1, 0,  2, 1, 0,  3, 1, 0;
  sim.SetVertices(x, VectorX::Ones(4));
  sim.SetFixed(0);
  sim.Finalize();
  EXPECT_EQ(sim.NumContactSlots(), 0u);  // No colliders, no slots.

  sim.AddCollider(Collider::Plane(EigenVector3::Zero(), EigenVector3::UnitY()));
  EXPECT_FALSE(sim.IsReady());
  sim.RebuildSystemMatrix();
  EXPECT_EQ(sim.NumContactSlots(), 3u);  // Every free vertex by default.

  sim.SetContactVertices({0, 2, 2, 3, 17});  // Fixed, duplicate and out-of-range ones are skipped.
  sim.RebuildSystemMatrix();
  EXPECT_EQ(sim.NumContactSlots(), 2u);
  // Outside every collider the slots cost nothing.
  EXPECT_EQ(sim.ConstraintEnergy(sim.GetPositions()), 0.0);
}

// A dropped particle comes to rest on the ground: shallow penetration, no
// bounce, and the energy never increases within a step.
TEST(ContactTest, ParticleRestsOnGround) {
  PDSolver sim;
  sim.m_track_energy = true;
  BuildParticle(&sim, EigenVector3(0.0, 0.5, 0.0), 0.1);
  sim.AddCollider(Collider::Plane(EigenVector3::Zero(), EigenVector3::UnitY()));
  sim.Finalize();
  ScalarType max_height_after_impact = -1.0;
  bool hit = false;
  for (int step = 0; step < 180; ++step) {
    sim.Step();
    const std::vector<ScalarType>& e = sim.GetStats().energy;
    for (size_t k = 1; k < e.size(); ++k) ASSERT_LE(e[k], e[k - 1] * (1.0 + 1e-9) + 1e-12) << "step " << step;
    hit |= sim.GetStats().num_contacts > 0;
    if (hit) max_height_after_impact = std::max(max_height_after_impact, sim.GetPositions()(0, 1));
  }
  ASSERT_TRUE(hit);
  EXPECT_LT(max_height_after_impact, 1e-3) << "particle bounced";
  // At rest the slot holds the weight: w * depth = m g.
  const ScalarType depth = -sim.GetPositions()(0, 1);
  EXPECT_NEAR(depth, 0.1 * sim.m_gravity_constant / sim.m_stiffness_contact, 1e-6);
  EXPECT_LT(sim.GetVelocities().norm(), 1e-6);
}

// Coulomb friction: a particle sliding on the ground at v0 stops after about
// v0^2 / (2 mu g); without friction it keeps sliding.
TEST(ContactTest, CoulombFrictionStopsSliding) {
  const ScalarType v0 = 2.0;
  ScalarType distance[2];
  const ScalarType mus[2] = {0.0, 0.5};
  for (int k = 0; k < 2; ++k) {
    PDSolver sim;
    sim.m_damping_coefficient = 0.0;
    BuildParticle(&sim, EigenVector3::Zero(), 1.0);
    const int ground = sim.AddCollider(Collider::Plane(EigenVector3::Zero(), EigenVector3::UnitY()));
    sim.GetCollider(ground).friction = mus[k];
    sim.Finalize();
    // Settle into the contact, then launch along x.
    for (int i = 0; i < 30; ++i) sim.Step();
    sim.SetState(sim.GetPositions(), (EigenMatrix3x(1, 3) << v0, 0.0, 0.0).finished());
    const ScalarType x0 = sim.GetPositions()(0, 0);
    for (int i = 0; i < 120; ++i) sim.Step();
    distance[k] = sim.GetPositions()(0, 0) - x0;
    if (k == 1) EXPECT_LT(std::abs(sim.GetVelocities()(0, 0)), 1e-9) << "did not stop";
  }
  EXPECT_NEAR(distance[0], v0 * 2.0, 0.05);  // 2 s at constant speed.
  const ScalarType expected = v0 * v0 / (2.0 * mus[1] * 9.81);
  EXPECT_NEAR(distance[1], expected, 0.15 * expected);
}

// Moving a collider only changes the projection targets: a sphere swept
// through a resting particle pushes it out of the way, with no rebuild.
TEST(ContactTest, MovingColliderPushesParticle) {
  PDSolver sim;
  sim.m_gravity_constant = 0.0;
  BuildParticle(&sim, EigenVector3::Zero(), 1.0);
  const int ball = sim.AddCollider(Collider::Sphere(EigenVector3(-1.0, 0.05, 0.0), 0.3));
  sim.Finalize();
  for (int step = 0; step <= 60; ++step) {
    sim.GetCollider(ball).translation.x() = -1.0 + step / 30.0;  // Ends at x = 1.
    sim.Step();
    ASSERT_TRUE(sim.IsReady());
  }
  const EigenVector3 x = sim.GetPositions().row(0).transpose();
  EXPECT_GT(sim.GetColliders()[ball].SignedDistance(x, nullptr), -0.02);  // Not left inside.
  EXPECT_GT(std::abs(x.y()) + x.x(), 0.1);                                 // Pushed along.
}

// Serial and parallel steps must agree with contact slots too.
TEST(ContactTest, ParallelMatchesSerialWithContact) {
  EigenMatrix3x result[2];
  for (int parallel = 0; parallel < 2; ++parallel) {
    TetScene scene;
    scene.m_stiffness_strain = 1e4;
    PDSolver sim;
    sim.m_parallel_local_step = parallel == 1;
    scene.Build(&sim);
    for (int i = 0; i < 60; ++i) sim.Step();
    result[parallel] = sim.GetPositions();
  }
  EXPECT_LT((result[0] - result[1]).cwiseAbs().maxCoeff(), 1e-12);
}

// --- Scenes -------------------------------------------------------------------------

// The default beam is stiff enough to clear the sphere. A soft one drapes over
// it onto the ground without its surface sinking into either, the energy never
// rising within a step and no tet inverting.
TEST(ContactSceneTest, SoftBeamDrapesOverSphere) {
  TetScene scene;
  scene.m_stiffness_strain = 1e4;
  PDSolver sim;
  sim.m_track_energy = true;
  sim.m_acceleration = PDAcceleration::kAnderson;
  scene.Build(&sim);
  EXPECT_GT(sim.NumContactSlots(), 0u);
  EXPECT_LT(sim.NumContactSlots(), sim.NumVertices());  // Boundary vertices only.
  for (int step = 0; step < 240; ++step) {
    scene.Update(&sim, sim.m_h);
    sim.Step();
    const std::vector<ScalarType>& e = sim.GetStats().energy;
    for (size_t k = 1; k < e.size(); ++k) {
      ASSERT_LE(e[k], e[k - 1] * (1.0 + 1e-8) + 1e-12) << "step " << step << " iteration " << k;
    }
  }
  const EigenMatrix3x& x = sim.GetPositions();
  EXPECT_GT(sim.GetStats().num_contacts, 0);
  EXPECT_LT(sim.GetStats().max_penetration, 5e-3);
  ScalarType lowest = x.col(1).minCoeff();
  EXPECT_LT(lowest, scene.m_ground_height + 0.01) << "beam did not reach the ground";
  EXPECT_GT(lowest, scene.m_ground_height - 0.01) << "beam sank through the ground";
  for (size_t t = 0; t < scene.Mesh().NumTets(); ++t) {
    EXPECT_GT(scene.Mesh().SignedVolume(t, x), 0.0) << "tet " << t << " inverted";
  }
}

// The default hair patch settles on the sphere and the ground.
TEST(ContactSceneTest, HairSettlesOnObstacles) {
  HairScene scene;
  PDSolver sim;
  scene.Build(&sim);
  ASSERT_EQ(sim.GetColliders().size(), 2u);
  for (int step = 0; step < 300; ++step) sim.Step();
  EXPECT_GT(sim.GetStats().num_contacts, 10);
  EXPECT_LT(sim.GetStats().max_penetration, 2e-3);
  EXPECT_GT(sim.GetPositions().col(1).minCoeff(), scene.m_ground_height - 2e-3);

  // Lowering the ground is RHS only: the strands follow it down.
  scene.m_ground_height = 0.0;
  for (int step = 0; step < 120; ++step) {
    scene.Update(&sim, sim.m_h);
    sim.Step();
  }
  EXPECT_LT(sim.GetPositions().col(1).minCoeff(), 0.6);
}

// Muscle volume over rest volume at solver positions x.
ScalarType VolumeRatio(const MuscleScene& scene, const EigenMatrix3x& x) {
  ScalarType v = 0.0;
  for (size_t t = 0; t < scene.Mesh().NumTets(); ++t) v += scene.Mesh().SignedVolume(t, x);
  return v / scene.Mesh().TotalVolume();
}

// Deepest penetration of any free muscle vertex into either bone, measured
// whether or not the bone colliders are enabled.
ScalarType MaxBonePenetration(const MuscleScene& scene, const PDSolver& sim) {
  std::vector<Collider> bones = sim.GetColliders();
  for (Collider& c : bones) c.enabled = true;
  ScalarType worst = 0.0;
  for (unsigned int i = 0; i < scene.Mesh().NumVertices(); ++i) {
    if (sim.IsFixed(i)) continue;
    ScalarType depth = 0.0;
    DeepestPenetration(bones, sim.GetPositions().row(i).transpose(), &depth, nullptr);
    worst = std::max(worst, depth);
  }
  return worst;
}

// Activated flexion to 120 degrees with gravity: the contracting belly is
// pulled down onto the upper arm. Without contact it sinks deep into the bone;
// with contact it stays (nearly) outside both bones and keeps its volume.
TEST(ContactSceneTest, MuscleStaysOutsideBones) {
  ScalarType penetration[2];
  for (int contact = 0; contact < 2; ++contact) {
    MuscleScene scene;
    scene.m_animate = false;
    scene.m_bone_contact = contact == 1;
    scene.m_activation = 1.0;
    PDSolver sim;
    sim.m_acceleration = PDAcceleration::kAnderson;  // As in the app.
    scene.Build(&sim);
    ScalarType worst = 0.0, worst_volume = 0.0;
    for (int step = 0; step < 240; ++step) {
      scene.m_elbow_angle_deg = std::min(1.0, step / 120.0) * MuscleScene::kMaxAnimatedAngleDeg;
      scene.Update(&sim, sim.m_h);
      sim.Step();
      worst = std::max(worst, MaxBonePenetration(scene, sim));
      worst_volume = std::max(worst_volume, std::abs(VolumeRatio(scene, sim.GetPositions()) - 1.0));
    }
    penetration[contact] = worst;
    RecordProperty(contact ? "volume_error_with_contact" : "volume_error_without_contact",
                   std::to_string(worst_volume));
    EXPECT_LT(worst_volume, 0.03) << "contact " << contact;
  }
  RecordProperty("penetration_without_contact", std::to_string(penetration[0]));
  RecordProperty("penetration_with_contact", std::to_string(penetration[1]));
  EXPECT_GT(penetration[0], 0.008);
  EXPECT_LT(penetration[1], 0.7 * penetration[0]);
  EXPECT_LT(penetration[1], 0.008);
}

// The full flex/extend animation with bone contact, plain PD: the energy never
// rises within a step (the colliders only move between steps).
TEST(ContactSceneTest, MuscleAnimationEnergyMonotone) {
  MuscleScene scene;
  PDSolver sim;
  sim.m_track_energy = true;
  scene.Build(&sim);
  EXPECT_EQ(sim.NumContactSlots(), sim.NumVertices() - sim.NumFixed());  // Interior contact by default.
  int contacts = 0;
  for (int step = 0; step < 240; ++step) {
    scene.Update(&sim, sim.m_h);
    sim.Step();
    contacts = std::max(contacts, sim.GetStats().num_contacts);
    const std::vector<ScalarType>& e = sim.GetStats().energy;
    for (size_t k = 1; k < e.size(); ++k) {
      ASSERT_LE(e[k], e[k - 1] * (1.0 + 1e-8) + 1e-12) << "step " << step << " iteration " << k;
    }
  }
  EXPECT_GT(contacts, 0);
}

// The forearm collider follows the bone (the insertion tendon keeps its
// distance to it) and carries its angular velocity.
TEST(ContactSceneTest, ForearmColliderFollowsBone) {
  MuscleScene scene;
  scene.m_animate = false;
  PDSolver sim;
  scene.Build(&sim);
  ASSERT_EQ(sim.GetColliders().size(), 2u);
  for (const Collider& c : sim.GetColliders()) {
    EXPECT_EQ(c.shape, Collider::Shape::kCapsule);
    EXPECT_TRUE(c.visible);  // The app draws the bones from their colliders.
  }
  const KinematicBone& bone = scene.Forearm();
  std::vector<ScalarType> rest_distance;
  for (unsigned int v : bone.vertices) {
    rest_distance.push_back(sim.GetColliders()[1].SignedDistance(sim.GetPositions().row(v).transpose(), nullptr));
  }
  const ScalarType dt = 0.1;
  scene.m_elbow_angle_deg = 30.0;
  scene.Update(&sim, dt);
  const Collider& fore = sim.GetColliders()[1];
  for (size_t k = 0; k < bone.vertices.size(); ++k) {
    const EigenVector3 x = sim.GetPositions().row(bone.vertices[k]).transpose();
    EXPECT_NEAR(fore.SignedDistance(x, nullptr), rest_distance[k], 1e-12);
  }
  EXPECT_NEAR(fore.angular_velocity.z(), (30.0 * kPi / 180.0) / dt, 1e-9);
  // The bones do not overlap up to the animated maximum flexion.
  scene.m_elbow_angle_deg = MuscleScene::kMaxAnimatedAngleDeg;
  scene.Update(&sim, dt);
  const Collider& upper = sim.GetColliders()[0];
  std::vector<EigenVector3> points;
  std::vector<unsigned int> tris;
  TriangulateCollider(sim.GetColliders()[1], 1.0, &points, &tris);
  for (const EigenVector3& p : points) EXPECT_GT(upper.SignedDistance(p, nullptr), 0.0);
}
