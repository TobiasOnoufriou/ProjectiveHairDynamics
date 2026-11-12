// Copyright 2022 Tobias Onoufriou
// Tests for the Phase 4 muscle model: fibre directions, the activated fibre
// constraint (Saito, Zhou, Kavan 2015) and the biceps scene with kinematic bones.
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "gtest/gtest.h"
#include "MuscleScene.h"
#include "PDSolver.h"
#include "TetConstraints.h"
#include "TetMesh.h"

namespace {

const std::array<unsigned int, 4> kTet = {0, 1, 2, 3};
constexpr double kPi = 3.14159265358979323846;

EigenMatrix3x RestTet() {
  EigenMatrix3x x(4, 3);
  x << 0.0, 0.0, 0.0,
       1.0, 0.1, 0.0,
       0.2, 0.9, 0.1,
       0.1, 0.3, 1.2;
  return x;
}

EigenMatrix3x Transform(const EigenMatrix3x& x, const EigenMatrix3& a, const EigenVector3& t) {
  EigenMatrix3x out = x * a.transpose();
  out.rowwise() += t.transpose();
  return out;
}

SparseMatrix GOf(const Constraint& c, int cols) {
  std::vector<SparseMatrixTriplet> triplets;
  c.AddG(0, &triplets);
  SparseMatrix g(c.NumRows(), cols);
  g.setFromTriplets(triplets.begin(), triplets.end());
  return g;
}

// Muscle volume (signed, so inversions count against it) at solver positions x.
ScalarType MuscleVolume(const TetMesh& mesh, const EigenMatrix3x& x) {
  ScalarType v = 0.0;
  for (size_t t = 0; t < mesh.NumTets(); ++t) v += mesh.SignedVolume(t, x);
  return v;
}

// Half extent and centre of the muscle along z. The muscle bends in the x-y
// plane, so its z extent is its thickness across the bending plane, and a
// shifted centre means it buckled sideways.
void MuscleZExtent(const TetMesh& mesh, const EigenMatrix3x& x, ScalarType* half_width, ScalarType* centre) {
  const ScalarType lo = x.topRows(mesh.NumVertices()).col(2).minCoeff();
  const ScalarType hi = x.topRows(mesh.NumVertices()).col(2).maxCoeff();
  *half_width = 0.5 * (hi - lo);
  *centre = 0.5 * (hi + lo);
}

// Mean fibre stretch |F a| over the belly (non-tendon) tets.
ScalarType MeanBellyFibreStretch(const MuscleScene& scene, const EigenMatrix3x& x) {
  const TetMesh& mesh = scene.Mesh();
  ScalarType sum = 0.0;
  int count = 0;
  for (size_t t = 0; t < mesh.NumTets(); ++t) {
    if (scene.IsTendon(t)) continue;
    ScalarType volume = 0.0;
    const Eigen::Matrix<ScalarType, 3, 4> g = TetConstraint::RestG(mesh.tets[t], mesh.vertices, &volume);
    EigenMatrix3x xt(4, 3);
    for (int c = 0; c < 4; ++c) xt.row(c) = x.row(mesh.tets[t][c]);
    const EigenMatrix3 f = (g * xt).transpose();
    sum += (f * mesh.fibres[t]).norm();
    ++count;
  }
  return sum / count;
}

}  // namespace

TEST(FibreTest, BeamFibresFollowLongestAxis) {
  TetMesh mesh = MakeTetBeam(2, 1, 1, EigenVector3(0.2, 1.0, 0.3));
  ASSERT_EQ(mesh.fibres.size(), mesh.NumTets());
  for (const EigenVector3& a : mesh.fibres) EXPECT_EQ(a, EigenVector3::UnitY());

  mesh.SetFibresFromField([](const EigenVector3& c) { return EigenVector3(3.0, 4.0 * c.y(), 0.0); });
  for (size_t t = 0; t < mesh.NumTets(); ++t) {
    EXPECT_NEAR(mesh.fibres[t].norm(), 1.0, 1e-14);
    EXPECT_NEAR(mesh.fibres[t].y() / mesh.fibres[t].x(), 4.0 * mesh.Centroid(t).y() / 3.0, 1e-12);
  }
}

// G x is the deformed fibre F a for any affine map x -> A x + t.
TEST(FibreTest, GMapsPositionsToDeformedFibre) {
  ScalarType k = 1.0, alpha = 0.0, c_max = 0.3;
  const EigenVector3 a = EigenVector3(1.0, 2.0, -0.5).normalized();
  TetFibreConstraint c(&k, kTet, RestTet(), 3.0 * a, &alpha, &c_max);
  EXPECT_NEAR(c.Fibre().norm(), 1.0, 1e-15);

  EigenMatrix3 f;
  f << 1.2, 0.3, -0.1,
       0.0, 0.9, 0.2,
       0.4, -0.2, 1.1;
  const EigenMatrix3x x = Transform(RestTet(), f, EigenVector3(3, -1, 2));
  const EigenVector3 fa = (GOf(c, 4) * x).transpose();
  EXPECT_LT((fa - f * a).cwiseAbs().maxCoeff(), 1e-12);
}

// Projection rescales F a to the target stretch 1 - alpha * c_max, keeping its
// direction. Activation changes only the projection, never G or the weight.
TEST(FibreTest, ProjectionAndActivation) {
  ScalarType k = 7.0, alpha = 0.0, c_max = 0.4;
  TetFibreConstraint c(&k, kTet, RestTet(), EigenVector3::UnitX(), &alpha, &c_max);
  const SparseMatrix g0 = GOf(c, 4);
  const ScalarType w0 = c.Weight();

  RowMatrix3x gx(1, 3), p(1, 3);
  gx << 0.3, -1.2, 0.4;
  for (ScalarType a : {0.0, 0.5, 1.0, 1.5, -0.2}) {
    alpha = a;
    const ScalarType target = 1.0 - std::min(std::max(a, 0.0), 1.0) * c_max;
    EXPECT_NEAR(c.TargetStretch(), target, 1e-15);
    c.Project(gx, p);
    EXPECT_NEAR(p.row(0).norm(), target, 1e-12);
    EXPECT_NEAR(p.row(0).dot(gx.row(0)), target * gx.row(0).norm(), 1e-12);  // Same direction.
    EXPECT_EQ((GOf(c, 4) - g0).norm(), 0.0);
    EXPECT_EQ(c.Weight(), w0);
  }
  // Degenerate fibre: any vector of the target length is a closest point.
  gx.setZero();
  c.Project(gx, p);
  EXPECT_NEAR(p.row(0).norm(), c.TargetStretch(), 1e-15);
}

// No energy at rest without activation, also after rigid motions; activation
// makes the rest shape cost w/2 (alpha c_max)^2.
TEST(FibreTest, EnergyAtRestAndUnderActivation) {
  ScalarType k = 1e3, alpha = 0.0, c_max = 0.3;
  const EigenMatrix3x rest = RestTet();
  PDSolver sim;
  sim.SetVertices(rest, VectorX::Ones(4));
  auto fibre = std::make_unique<TetFibreConstraint>(&k, kTet, rest, EigenVector3(1, 1, 0), &alpha, &c_max);
  const ScalarType w = fibre->Weight();
  sim.AddConstraint(std::move(fibre));
  sim.Finalize();

  const EigenMatrix3 r = Eigen::AngleAxis<ScalarType>(1.1, EigenVector3(1, -2, 0.5).normalized()).toRotationMatrix();
  EXPECT_NEAR(sim.ConstraintEnergy(rest), 0.0, 1e-24);
  EXPECT_NEAR(sim.ConstraintEnergy(Transform(rest, r, EigenVector3(1, 2, 3))), 0.0, 1e-20);
  alpha = 1.0;
  EXPECT_NEAR(sim.ConstraintEnergy(rest), 0.5 * w * c_max * c_max, 1e-12 * w);
  EXPECT_NEAR(sim.ConstraintEnergy(Transform(rest, r, EigenVector3(1, 2, 3))), 0.5 * w * c_max * c_max, 1e-10 * w);
}

// A free, weightless block with fibres along x and volume preservation: full
// activation shortens it along the fibres and, keeping its volume, makes it
// thicker across them (the muscle bulge).
TEST(FibreTest, ActivatedBlockShortensAndBulges) {
  TetMesh mesh = MakeTetBeam(6, 2, 2, EigenVector3(1.0, 0.3, 0.3));
  ScalarType k_strain = 1e3, k_volume = 1e6, k_fibre = 1e5, alpha = 1.0, c_max = 0.3, one = 1.0;
  PDSolver sim;
  sim.m_gravity_constant = 0.0;
  sim.m_iterations_per_step = 20;
  sim.SetVertices(mesh.vertices, mesh.LumpedMasses(1000.0));
  for (size_t t = 0; t < mesh.NumTets(); ++t) {
    const TetMesh::Tet& tet = mesh.tets[t];
    sim.AddConstraint(std::make_unique<TetStrainConstraint>(&k_strain, tet, mesh.vertices, &one, &one));
    sim.AddConstraint(std::make_unique<TetVolumeConstraint>(&k_volume, tet, mesh.vertices, &one, &one));
    sim.AddConstraint(std::make_unique<TetFibreConstraint>(&k_fibre, tet, mesh.vertices, mesh.fibres[t], &alpha,
                                                           &c_max));
  }
  sim.Finalize();
  for (int step = 0; step < 300; ++step) sim.Step();

  const EigenMatrix3x& x = sim.GetPositions();
  const EigenVector3 extent = x.colwise().maxCoeff() - x.colwise().minCoeff();
  EXPECT_LT(extent.x(), 0.85);  // Rest length 1, fibre target 0.7.
  EXPECT_GT(extent.x(), 0.65);
  EXPECT_GT(extent.y(), 0.32);  // Rest thickness 0.3.
  EXPECT_GT(extent.z(), 0.32);
  EXPECT_NEAR(MuscleVolume(mesh, x) / mesh.TotalVolume(), 1.0, 0.02);
}

TEST(MuscleSceneTest, BuildsValidMuscleAndBones) {
  MuscleScene scene;
  PDSolver sim;
  scene.Build(&sim);
  ASSERT_TRUE(sim.IsReady());
  const TetMesh& mesh = scene.Mesh();

  size_t tendon = 0;
  for (size_t t = 0; t < mesh.NumTets(); ++t) {
    EXPECT_GT(mesh.SignedVolume(t), 1e-9) << "tet " << t;
    EXPECT_NEAR(mesh.fibres[t].norm(), 1.0, 1e-12);
    EXPECT_GT(mesh.fibres[t].x(), 0.5);  // Fibres run along the muscle.
    tendon += scene.IsTendon(t) ? 1 : 0;
  }
  EXPECT_GT(tendon, 0u);
  EXPECT_LT(tendon, mesh.NumTets() / 2);

  // Muscle vertices come first; the rest are bone box corners, all fixed.
  EXPECT_EQ(sim.NumVertices(), mesh.NumVertices() + 48);
  for (unsigned int i = mesh.NumVertices(); i < sim.NumVertices(); ++i) EXPECT_TRUE(sim.IsFixed(i));
  for (const KinematicBone* bone : {&scene.UpperArm(), &scene.Forearm()}) {
    EXPECT_GT(bone->vertices.size(), 24u);  // Box plus attached tendon vertices.
    for (unsigned int v : bone->vertices) EXPECT_TRUE(sim.IsFixed(v));
  }
  // The muscle surface is closed and the boxes are 12 triangles each.
  EXPECT_EQ(scene.SurfaceIndices().size(), 3 * (mesh.BoundaryTriangles().size() + 24));
  EXPECT_NEAR(sim.ConstraintEnergy(sim.GetPositions()), 0.0, 1e-16);
}

// Posing the elbow moves the forearm and the insertion tendon rigidly about the
// elbow and leaves the upper arm where it is.
TEST(MuscleSceneTest, ElbowRotatesForearmOnly) {
  MuscleScene scene;
  scene.m_animate = false;
  PDSolver sim;
  scene.Build(&sim);
  const EigenMatrix3x before = sim.GetPositions();
  scene.m_elbow_angle_deg = 90.0;
  scene.Update(&sim, 0.0);
  const EigenMatrix3x& x = sim.GetPositions();
  const KinematicBone& fore = scene.Forearm();
  for (size_t k = 0; k < fore.vertices.size(); ++k) {
    const EigenVector3 r = fore.rest.row(k).transpose() - fore.pivot;
    const EigenVector3 expected = fore.pivot + EigenVector3(-r.y(), r.x(), r.z());  // 90 deg about +z.
    EXPECT_LT((x.row(fore.vertices[k]).transpose() - expected).norm(), 1e-12);
  }
  for (unsigned int v : scene.UpperArm().vertices) EXPECT_EQ(x.row(v), before.row(v));
}

// Animated flex/extend with activation: moving targets between steps, but within
// each step the PD energy must still never increase, and the muscle keeps its volume.
TEST(MuscleSceneTest, AnimationKeepsEnergyMonotoneAndVolume) {
  MuscleScene scene;
  PDSolver sim;
  sim.m_track_energy = true;
  scene.Build(&sim);
  const ScalarType rest_volume = scene.Mesh().TotalVolume();
  for (int step = 0; step < 240; ++step) {  // One full 4 s cycle at 60 Hz.
    scene.Update(&sim, sim.m_h);
    sim.Step();
    const std::vector<ScalarType>& e = sim.GetStats().energy;
    for (size_t k = 1; k < e.size(); ++k) {
      ASSERT_LE(e[k], e[k - 1] * (1.0 + 1e-8) + 1e-12) << "step " << step << " iteration " << k;
    }
    ASSERT_NEAR(MuscleVolume(scene.Mesh(), sim.GetPositions()) / rest_volume, 1.0, 0.05) << "step " << step;
  }
}

// Flexing the elbow to 120 degrees shortens the muscle path. Activated, the
// belly fibres contract to take up the slack and the belly thickens (bulges)
// symmetrically. Passive, the fibres keep their length and the slack muscle
// buckles sideways instead.
TEST(MuscleSceneTest, ActiveFlexionContractsAndBulges) {
  ScalarType width[2], offset[2], stretch[2];
  for (int active = 0; active < 2; ++active) {
    MuscleScene scene;
    scene.m_animate = false;
    PDSolver sim;
    sim.m_gravity_constant = 0.0;
    scene.Build(&sim);
    ScalarType rest_width = 0.0, rest_centre = 0.0;
    MuscleZExtent(scene.Mesh(), sim.GetPositions(), &rest_width, &rest_centre);
    scene.m_activation = active;
    for (int step = 0; step < 240; ++step) {
      scene.m_elbow_angle_deg = std::min(1.0, step / 120.0) * MuscleScene::kMaxAnimatedAngleDeg;
      scene.Update(&sim, sim.m_h);
      sim.Step();
    }
    MuscleZExtent(scene.Mesh(), sim.GetPositions(), &width[active], &offset[active]);
    width[active] /= rest_width;
    stretch[active] = MeanBellyFibreStretch(scene, sim.GetPositions());
    EXPECT_NEAR(MuscleVolume(scene.Mesh(), sim.GetPositions()) / scene.Mesh().TotalVolume(), 1.0, 0.03);
  }
  RecordProperty("passive_width", std::to_string(width[0]));
  RecordProperty("active_width", std::to_string(width[1]));
  RecordProperty("passive_fibre_stretch", std::to_string(stretch[0]));
  RecordProperty("active_fibre_stretch", std::to_string(stretch[1]));

  EXPECT_LT(stretch[1], 0.93);               // Active fibres contract...
  EXPECT_GT(stretch[0], 0.97);               // ...passive ones do not.
  EXPECT_GT(width[1], 1.1);                  // Active belly bulges...
  EXPECT_GT(width[1], width[0]);             // ...more than the passive one...
  EXPECT_LT(std::abs(offset[1]), 0.005);     // ...and stays centred,
  EXPECT_GT(std::abs(offset[0]), std::abs(offset[1]));  // while the passive one buckles sideways.
}
