// Copyright 2022 Tobias Onoufriou
// Tests for the Phase 5 solver work: the signed SVD, the flattened per-type
// constraint batches (including strain + volume pairs sharing one SVD), Chebyshev / Anderson acceleration and the
// allocation-free step.
#include <algorithm>
#include <cmath>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include "gtest/gtest.h"
#include "HairScene.h"
#include "MuscleScene.h"
#include "PDSolver.h"
#include "TetConstraints.h"
#include "TetMesh.h"
#include "TetScene.h"

namespace {

// Forwards everything to a wrapped constraint but reports ConstraintKind::kGeneric,
// so the solver projects it through the virtual interface instead of its flat
// kernels.
class VirtualOnly : public Constraint {
 public:
  explicit VirtualOnly(std::unique_ptr<Constraint> c) : Constraint(nullptr), m_c(std::move(c)) {}
  ScalarType Weight() const override { return m_c->Weight(); }
  int NumRows() const override { return m_c->NumRows(); }
  void AddG(int first_row, std::vector<SparseMatrixTriplet>* t) const override { m_c->AddG(first_row, t); }
  void Project(const ConstRows& gx, Rows p) const override { m_c->Project(gx, p); }
  ScalarType Strain(const ConstRows& gx) const override { return m_c->Strain(gx); }
  ScalarType VolumeError(const ConstRows& gx) const override { return m_c->VolumeError(gx); }

 private:
  std::unique_ptr<Constraint> m_c;
};

// A block with every built-in constraint kind (strain, volume, fibre, springs),
// pinned at x = 0 and sagging under gravity; activation contracts the fibres.
struct MixedBlock {
  TetMesh mesh = MakeTetBeam(6, 2, 2, EigenVector3(1.0, 0.25, 0.25));
  ScalarType k_strain = 1e4, k_volume = 1e5, k_fibre = 1e4, k_spring = 1e3;
  ScalarType alpha = 0.6, c_max = 0.3, lo = 0.95, hi = 1.05, one = 1.0;

  void Build(PDSolver* sim, bool virtual_only) {
    auto add = [&](std::unique_ptr<Constraint> c) {
      if (virtual_only) c = std::make_unique<VirtualOnly>(std::move(c));
      sim->AddConstraint(std::move(c));
    };
    sim->SetVertices(mesh.vertices, mesh.LumpedMasses(1000.0));
    for (size_t t = 0; t < mesh.NumTets(); ++t) {
      const TetMesh::Tet& tet = mesh.tets[t];
      add(std::make_unique<TetStrainConstraint>(&k_strain, tet, mesh.vertices, &lo, &hi));
      add(std::make_unique<TetVolumeConstraint>(&k_volume, tet, mesh.vertices, &one, &one));
      add(std::make_unique<TetFibreConstraint>(&k_fibre, tet, mesh.vertices, mesh.fibres[t], &alpha, &c_max));
    }
    for (const auto& e : mesh.Edges()) {
      const ScalarType len = (mesh.vertices.row(e[0]) - mesh.vertices.row(e[1])).norm();
      add(std::make_unique<SpringConstraint>(&k_spring, e[0], e[1], len));
    }
    for (unsigned int i = 0; i < mesh.NumVertices(); ++i) {
      if (mesh.vertices(i, 0) < 1e-9) sim->SetFixed(i);
    }
    sim->Finalize();
  }
};

// Runs a scene for `warmup` plain steps, then one step with the given mode and
// iteration count; returns the energy gap E_final - E* relative to E_0 - E*,
// where E* comes from a 2000-iteration plain step from the same state.
template <typename SceneT>
ScalarType RelativeEnergyGap(PDAcceleration mode, int iterations, int warmup) {
  ScalarType final_energy[2], first_energy = 0.0;
  for (int run = 0; run < 2; ++run) {
    SceneT scene;
    PDSolver sim;
    scene.Build(&sim);
    for (int s = 0; s < warmup; ++s) {
      scene.Update(&sim, sim.m_h);
      sim.Step();
    }
    scene.Update(&sim, sim.m_h);
    sim.m_track_energy = true;
    sim.m_acceleration = run == 0 ? PDAcceleration::kNone : mode;
    sim.m_iterations_per_step = run == 0 ? 2000 : iterations;
    sim.Step();
    final_energy[run] = sim.GetStats().energy.back();
    first_energy = sim.GetStats().energy.front();
  }
  return (final_energy[1] - final_energy[0]) / (first_energy - final_energy[0]);
}

bool NonIncreasing(const std::vector<ScalarType>& e, ScalarType* worst_rise) {
  *worst_rise = 0.0;
  for (size_t k = 1; k < e.size(); ++k) {
    *worst_rise = std::max(*worst_rise, (e[k] - e[k - 1]) / std::max(std::abs(e[k - 1]), 1e-12));
  }
  return *worst_rise <= 1e-9;
}

}  // namespace

// The signed SVD reproduces F with proper rotations U, V and |s| equal to
// the singular values, s(2) < 0 exactly when F is inverted, on random,
// inverted, nearly repeated and rank-deficient matrices.
TEST(SignedSVDTest, MatchesJacobiSVD) {
  std::mt19937 rng(7);
  std::uniform_real_distribution<ScalarType> uni(-1.0, 1.0);
  auto random_rotation = [&]() {
    return Eigen::Quaternion<ScalarType>(uni(rng), uni(rng), uni(rng), uni(rng)).normalized().toRotationMatrix();
  };
  std::vector<EigenMatrix3> cases = {EigenMatrix3::Identity(), EigenMatrix3::Zero(),
                                     EigenMatrix3(EigenVector3(1, 1, -1).asDiagonal())};
  for (int i = 0; i < 300; ++i) {
    EigenVector3 s(uni(rng) + 1.5, uni(rng) + 1.2, uni(rng));
    if (i % 5 == 1) s(1) = s(0) * (1.0 + 1e-9);   // Nearly repeated.
    if (i % 5 == 2) s(2) = 0.0;                   // Flat tet.
    if (i % 5 == 3) s(1) = s(2) = 0.0;            // Collapsed to a line.
    if (i % 5 == 4) s *= 1e-3;                    // Tiny scale.
    cases.push_back(random_rotation() * s.asDiagonal() * random_rotation());
  }
  for (size_t c = 0; c < cases.size(); ++c) {
    const EigenMatrix3& f = cases[c];
    EigenMatrix3 u, v;
    EigenVector3 s;
    TetConstraint::SignedSVD(f, &u, &s, &v);
    const ScalarType scale = std::max(f.norm(), 1e-300);
    EXPECT_LT((u * s.asDiagonal() * v.transpose() - f).norm() / scale, 1e-13) << "case " << c << "\n" << f;
    EXPECT_LT((u.transpose() * u - EigenMatrix3::Identity()).norm(), 1e-13);
    EXPECT_LT((v.transpose() * v - EigenMatrix3::Identity()).norm(), 1e-13);
    EXPECT_NEAR(u.determinant(), 1.0, 1e-13);
    EXPECT_NEAR(v.determinant(), 1.0, 1e-13);

    Eigen::JacobiSVD<EigenMatrix3> svd(f);
    EigenVector3 ref = svd.singularValues();
    if (f.determinant() < 0.0) ref(2) = -ref(2);
    EXPECT_LT((s - ref).norm() / scale, 1e-13) << "case " << c << "\n" << f;
    EXPECT_LT((TetConstraint::SignedSingularValues(f) - s).norm() / scale, 1e-13);
  }
}

// The solver's flat per-type kernels give the same energies, metrics and
// trajectory as projecting every constraint through its virtual Project().
TEST(FlatBatchTest, MatchesVirtualProjection) {
  MixedBlock scenes[2];
  PDSolver sims[2];
  for (int v = 0; v < 2; ++v) {
    sims[v].m_iterations_per_step = 15;
    scenes[v].Build(&sims[v], v == 1);
  }
  ASSERT_EQ(sims[0].NumConstraints(), sims[1].NumConstraints());
  for (int step = 0; step < 40; ++step) {
    if (step == 10) {
      for (PDSolver& sim : sims) sim.BeginDrag(sim.NumVertices() - 1, EigenVector3(1.2, 0.4, 0.2));
    }
    for (PDSolver& sim : sims) sim.Step();
  }
  const EigenMatrix3x& x = sims[0].GetPositions();
  EXPECT_LT((x - sims[1].GetPositions()).cwiseAbs().maxCoeff(), 1e-10);
  EXPECT_NEAR(sims[0].ConstraintEnergy(x), sims[1].ConstraintEnergy(x), 1e-9 * sims[1].ConstraintEnergy(x));
  EXPECT_NEAR(sims[0].MaxStrain(x), sims[1].MaxStrain(x), 1e-7);
  EXPECT_NEAR(sims[0].MaxVolumeError(x), sims[1].MaxVolumeError(x), 1e-12);
  EXPECT_GT(sims[0].MaxStrain(x), 0.0);
}

// Changing a value behind a constraint's pointer (here the activation) is seen
// by the flat batches without a rebuild.
TEST(FlatBatchTest, ActivationIsReadLive) {
  MixedBlock block;
  PDSolver sim;
  sim.m_gravity_constant = 0.0;
  block.Build(&sim, false);
  const EigenMatrix3x& x = sim.GetPositions();
  const ScalarType e_active = sim.ConstraintEnergy(x);
  block.alpha = 0.0;
  EXPECT_LT(sim.ConstraintEnergy(x), 1e-3 * e_active);
}

// Every acceleration is safeguarded, so the accepted energies never rise, also
// with a drag (generic batch) and moving bones.
TEST(AccelerationTest, EnergyNeverIncreases) {
  for (PDAcceleration mode : {PDAcceleration::kChebyshev, PDAcceleration::kAnderson}) {
    MuscleScene scene;
    scene.m_resolution = 8;
    PDSolver sim;
    sim.m_acceleration = mode;
    sim.m_track_energy = true;
    sim.m_chebyshev_rho = 0.99;  // Aggressive, so the safeguard is exercised.
    scene.Build(&sim);
    int rejected = 0;
    for (int step = 0; step < 60; ++step) {
      if (step == 30) sim.BeginDrag(0, EigenVector3(0.0, 0.1, 0.1));
      scene.Update(&sim, sim.m_h);
      sim.Step();
      ScalarType rise = 0.0;
      ASSERT_TRUE(NonIncreasing(sim.GetStats().energy, &rise))
          << "mode " << static_cast<int>(mode) << " step " << step << ": rise " << rise;
      rejected += sim.GetStats().rejected_iterations;
    }
    if (mode == PDAcceleration::kChebyshev) EXPECT_GT(rejected, 0);
  }
}

// At the same iteration count, both accelerations end closer to the converged
// step than plain PD (cantilever under gravity and the animated muscle).
TEST(AccelerationTest, ConvergesFasterThanPlain) {
  for (int iterations : {10, 20}) {
    const ScalarType plain = RelativeEnergyGap<TetScene>(PDAcceleration::kNone, iterations, 20);
    const ScalarType cheb = RelativeEnergyGap<TetScene>(PDAcceleration::kChebyshev, iterations, 20);
    const ScalarType anderson = RelativeEnergyGap<TetScene>(PDAcceleration::kAnderson, iterations, 20);
    EXPECT_LT(cheb, plain) << iterations << " iterations";
    EXPECT_LT(anderson, plain) << iterations << " iterations";
  }
  const ScalarType plain = RelativeEnergyGap<MuscleScene>(PDAcceleration::kNone, 10, 30);
  EXPECT_LT(RelativeEnergyGap<MuscleScene>(PDAcceleration::kChebyshev, 10, 30), plain);
  EXPECT_LT(RelativeEnergyGap<MuscleScene>(PDAcceleration::kAnderson, 10, 30), plain);
}

// Once its buffers are sized, Step() makes no Eigen heap allocation, in every
// acceleration mode, serial and parallel, with energy tracking and a drag.
TEST(PDSolverTest, StepDoesNotAllocate) {
  MuscleScene scene;
  scene.m_resolution = 8;
  PDSolver sim;
  sim.m_track_energy = true;
  scene.Build(&sim);
  sim.BeginDrag(0, EigenVector3(0.0, 0.1, 0.1));
  for (PDAcceleration mode : {PDAcceleration::kNone, PDAcceleration::kChebyshev, PDAcceleration::kAnderson}) {
    for (bool parallel : {false, true}) {
      sim.m_acceleration = mode;
      sim.m_parallel_local_step = parallel;
      sim.Step();  // Warm-up: sizes the Anderson history.
      const bool asserts = pd::g_eigen_asserts;
      pd::g_eigen_asserts = true;
      Eigen::internal::set_is_malloc_allowed(false);
      for (int i = 0; i < 3; ++i) sim.Step();  // Aborts on an allocation.
      Eigen::internal::set_is_malloc_allowed(true);
      pd::g_eigen_asserts = asserts;
    }
  }
  SUCCEED();
}
