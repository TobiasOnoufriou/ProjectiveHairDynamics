// Copyright 2022 Tobias Onoufriou
// Tests for the volumetric PD constraints (corotational strain and volume
// preservation, Bouaziz et al. 2014) and the race-free parallel local step.
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <vector>

#include "gtest/gtest.h"
#include "PDSolver.h"
#include "TetConstraints.h"
#include "TetMesh.h"
#include "TetScene.h"

namespace {

const std::array<unsigned int, 4> kTet = {0, 1, 2, 3};

// A positively oriented, deliberately irregular rest tet.
EigenMatrix3x RestTet() {
  EigenMatrix3x x(4, 3);
  x << 0.0, 0.0, 0.0,
       1.0, 0.1, 0.0,
       0.2, 0.9, 0.1,
       0.1, 0.3, 1.2;
  return x;
}

EigenMatrix3 Rotation(ScalarType angle, EigenVector3 axis) {
  return Eigen::AngleAxis<ScalarType>(angle, axis.normalized()).toRotationMatrix();
}

// Applies x -> R x + t to every row.
EigenMatrix3x Transform(const EigenMatrix3x& x, const EigenMatrix3& r, const EigenVector3& t) {
  EigenMatrix3x out = x * r.transpose();
  out.rowwise() += t.transpose();
  return out;
}

// Projection of a deformation gradient F through a constraint's Project().
EigenMatrix3 ProjectF(const Constraint& c, const EigenMatrix3& f) {
  RowMatrix3x gx = f.transpose();
  RowMatrix3x p(3, 3);
  c.Project(gx, p);
  return p.transpose();
}

// One free tet with both constraints; energy is evaluated through the solver.
struct SingleTet {
  ScalarType stiffness = 1e4, one = 1.0;
  PDSolver sim;
  explicit SingleTet(bool strain, bool volume) {
    const EigenMatrix3x rest = RestTet();
    sim.SetVertices(rest, VectorX::Ones(4));
    if (strain) sim.AddConstraint(std::make_unique<TetStrainConstraint>(&stiffness, kTet, rest, &one, &one));
    if (volume) sim.AddConstraint(std::make_unique<TetVolumeConstraint>(&stiffness, kTet, rest, &one, &one));
    sim.Finalize();
  }
};

bool NonIncreasing(const std::vector<ScalarType>& e, ScalarType tol, ScalarType* worst_rise) {
  *worst_rise = 0.0;
  for (size_t k = 1; k < e.size(); ++k) {
    *worst_rise = std::max(*worst_rise, (e[k] - e[k - 1]) / std::max(e[k - 1], 1e-12));
  }
  return *worst_rise <= tol;
}

}  // namespace

// G x reproduces the deformation gradient: an affine map x -> A x + t has F = A.
TEST(TetConstraintTest, GMapsPositionsToDeformationGradient) {
  ScalarType k = 1.0, one = 1.0;
  const EigenMatrix3x rest = RestTet();
  TetStrainConstraint c(&k, kTet, rest, &one, &one);
  std::vector<SparseMatrixTriplet> triplets;
  c.AddG(0, &triplets);
  SparseMatrix g(3, 4);
  g.setFromTriplets(triplets.begin(), triplets.end());

  EigenMatrix3 a;
  a << 1.2, 0.3, -0.1,
       0.0, 0.9, 0.2,
       0.4, -0.2, 1.1;
  const EigenMatrix3x x = Transform(rest, a, EigenVector3(3, -1, 2));
  const EigenMatrix3 f = (g * x).transpose();
  EXPECT_LT((f - a).cwiseAbs().maxCoeff(), 1e-12);
  EXPECT_NEAR(c.RestVolume(), std::abs(TetMesh{rest, {kTet}}.SignedVolume(0)), 1e-15);
}

TEST(TetConstraintTest, RestShapeHasZeroEnergy) {
  SingleTet tet(true, true);
  EXPECT_NEAR(tet.sim.ConstraintEnergy(RestTet()), 0.0, 1e-24);
  EXPECT_NEAR(tet.sim.MaxStrain(RestTet()), 0.0, 1e-12);
  EXPECT_NEAR(tet.sim.MaxVolumeError(RestTet()), 0.0, 1e-12);
}

// The energy only measures deformation: rigid motions of the rest shape cost
// nothing, and rotating a deformed shape does not change its energy.
TEST(TetConstraintTest, EnergyIsRotationInvariant) {
  for (bool strain : {true, false}) {
    SingleTet tet(strain, !strain);
    const EigenMatrix3x rest = RestTet();
    for (int i = 0; i < 5; ++i) {
      const EigenMatrix3 r = Rotation(0.7 + 1.3 * i, EigenVector3(1.0, -0.5 + i, 0.3 * i));
      const EigenVector3 t(0.5 * i, -1.0, 2.0);
      EXPECT_NEAR(tet.sim.ConstraintEnergy(Transform(rest, r, t)), 0.0, 1e-20) << "strain=" << strain;

      EigenMatrix3x deformed = rest;
      deformed.row(3) += EigenVector3(0.2, -0.1, 0.4).transpose();
      deformed.row(1) *= 0.8;
      const ScalarType e = tet.sim.ConstraintEnergy(deformed);
      ASSERT_GT(e, 1e-3);
      EXPECT_NEAR(tet.sim.ConstraintEnergy(Transform(deformed, r, t)), e, 1e-10 * e) << "strain=" << strain;
    }
  }
}

// With the range [1, 1] the strain projection is the closest proper rotation,
// also for inverted tets (no reflection). A wider range clamps singular values.
TEST(TetConstraintTest, StrainProjection) {
  ScalarType k = 1.0, lo = 1.0, hi = 1.0;
  TetStrainConstraint c(&k, kTet, RestTet(), &lo, &hi);

  const EigenMatrix3 rot = Rotation(0.4, EigenVector3(1, 2, 3));
  EigenMatrix3 f = rot * EigenVector3(1.3, 0.8, 1.1).asDiagonal() * Rotation(-1.1, EigenVector3(0, 1, 1));
  EigenMatrix3 p = ProjectF(c, f);
  EXPECT_LT((p.transpose() * p - EigenMatrix3::Identity()).cwiseAbs().maxCoeff(), 1e-12);
  EXPECT_NEAR(p.determinant(), 1.0, 1e-12);
  EXPECT_LT((p - rot * Rotation(-1.1, EigenVector3(0, 1, 1))).cwiseAbs().maxCoeff(), 1e-12);

  // Inverted: F = diag(1, 1, -0.5). Closest rotation is the identity.
  f = EigenVector3(1.0, 1.0, -0.5).asDiagonal();
  p = ProjectF(c, f);
  EXPECT_NEAR(p.determinant(), 1.0, 1e-12);
  EXPECT_LT((p - EigenMatrix3::Identity()).cwiseAbs().maxCoeff(), 1e-12);

  lo = 0.9;
  hi = 1.1;
  f = rot * EigenVector3(1.5, 1.0, 0.5).asDiagonal();
  p = ProjectF(c, f);
  EXPECT_LT((p - rot * EigenVector3(1.1, 1.0, 0.9).asDiagonal().toDenseMatrix()).cwiseAbs().maxCoeff(), 1e-12);
  // Inside the range nothing changes.
  f = rot * EigenVector3(1.05, 1.0, 0.95).asDiagonal();
  EXPECT_LT((ProjectF(c, f) - f).cwiseAbs().maxCoeff(), 1e-12);
}

// The volume projection lands on det = 1, keeps the singular vectors, and is
// closer to F than the naive uniform rescaling F / det(F)^(1/3).
TEST(TetConstraintTest, VolumeProjection) {
  ScalarType k = 1.0, one = 1.0;
  TetVolumeConstraint c(&k, kTet, RestTet(), &one, &one);
  const EigenMatrix3 rot = Rotation(0.9, EigenVector3(-1, 2, 0.5));
  for (const EigenVector3& s : {EigenVector3(2.0, 1.0, 1.0), EigenVector3(0.7, 0.9, 0.8),
                                EigenVector3(1.4, 1.2, 0.3), EigenVector3(1.0, 1.0, -0.2)}) {
    const EigenMatrix3 f = rot * s.asDiagonal();
    const EigenMatrix3 p = ProjectF(c, f);
    EXPECT_NEAR(p.determinant(), 1.0, 1e-9) << s.transpose();
    if (f.determinant() > 0) {
      const EigenMatrix3 scaled = f / std::cbrt(f.determinant());
      EXPECT_LE((f - p).norm(), (f - scaled).norm() + 1e-12) << s.transpose();
    }
    // p = rot * diag(d): same singular vectors as F.
    const EigenMatrix3 d = rot.transpose() * p;
    EXPECT_LT((d - EigenMatrix3(d.diagonal().asDiagonal())).cwiseAbs().maxCoeff(), 1e-9);
  }
  // Already volume preserving: unchanged.
  const EigenMatrix3 f = rot * EigenVector3(2.0, 0.5, 1.0).asDiagonal();
  EXPECT_LT((ProjectF(c, f) - f).cwiseAbs().maxCoeff(), 1e-12);
}

TEST(TetConstraintTest, UniformScaleMetrics) {
  SingleTet tet(true, true);
  const EigenMatrix3x scaled = 1.1 * RestTet();
  EXPECT_NEAR(tet.sim.MaxStrain(scaled), 0.1, 1e-12);
  EXPECT_NEAR(tet.sim.MaxVolumeError(scaled), 1.1 * 1.1 * 1.1 - 1.0, 1e-12);
}

// Weights are stiffness * rest volume, so refining the mesh does not change the
// total energy of the same affine deformation.
TEST(TetConstraintTest, EnergyIndependentOfResolution) {
  EigenMatrix3 a;
  a << 1.1, 0.1, 0.0,
       0.0, 0.95, 0.05,
       0.0, 0.0, 1.02;
  std::vector<ScalarType> energies;
  for (unsigned int n : {1u, 3u}) {
    TetMesh mesh = MakeTetBeam(n, n, n, EigenVector3(1, 1, 1));
    ScalarType k = 1.0, one = 1.0;
    PDSolver sim;
    sim.SetVertices(mesh.vertices, mesh.LumpedMasses(1.0));
    for (const TetMesh::Tet& t : mesh.tets) {
      sim.AddConstraint(std::make_unique<TetStrainConstraint>(&k, t, mesh.vertices, &one, &one));
    }
    sim.Finalize();
    energies.push_back(sim.ConstraintEnergy(Transform(mesh.vertices, a, EigenVector3::Zero())));
  }
  EXPECT_NEAR(energies[0], energies[1], 1e-10 * energies[0]);
}

// Exact projections keep the PD energy non-increasing over local/global
// iterations on the volumetric cantilever, while dragging too.
TEST(TetConstraintTest, CantileverEnergyNeverIncreases) {
  TetScene scene;
  PDSolver sim;
  sim.m_track_energy = true;
  sim.m_iterations_per_step = 20;
  scene.Build(&sim);
  for (int step = 0; step < 120; ++step) {
    if (step == 40) sim.BeginDrag(sim.NumVertices() - 1, EigenVector3(0.8, 1.2, 0.3));
    if (step == 80) sim.EndDrag();
    sim.Step();
    ScalarType rise = 0.0;
    ASSERT_TRUE(NonIncreasing(sim.GetStats().energy, 1e-8, &rise)) << "step " << step << ": rise " << rise;
  }
}

// Each constraint projects into its own slot, so the OpenMP local step must give
// bit-identical positions to the serial one.
TEST(TetConstraintTest, ParallelLocalStepMatchesSerial) {
  EigenMatrix3x result[2];
  for (int parallel = 0; parallel < 2; ++parallel) {
    TetScene scene;
    PDSolver sim;
    sim.m_parallel_local_step = parallel == 1;
    scene.Build(&sim);
    for (int step = 0; step < 60; ++step) {
      if (step == 20) sim.BeginDrag(sim.NumVertices() - 1, EigenVector3(0.8, 1.2, 0.3));
      sim.Step();
    }
    result[parallel] = sim.GetPositions();
  }
  EXPECT_EQ((result[0] - result[1]).cwiseAbs().maxCoeff(), 0.0);
}
