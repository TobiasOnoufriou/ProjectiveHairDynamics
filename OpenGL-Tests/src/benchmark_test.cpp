// Copyright 2022 Tobias Onoufriou
// Step-cost benchmark (disabled by default). Run with
//   OpenGL-Tests.exe --gtest_also_run_disabled_tests --gtest_filter=*Benchmark*
#include <chrono>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <limits>
#include <vector>

#include "gtest/gtest.h"
#include "HairScene.h"
#include "MuscleScene.h"
#include "PDSolver.h"
#include "TetConstraints.h"
#include "TetScene.h"

TEST(PDSolverBenchmark, DISABLED_StepCostVsIterations) {
  for (bool track : {false, true}) {
    for (int iters : {1, 2, 10, 50}) {
      HairScene scene;  // Default 5x5 patch, 525 vertices.
      PDSolver sim;
      sim.m_track_energy = track;
      sim.m_iterations_per_step = iters;
      scene.Build(&sim);
      const int steps = 200;
      const auto t0 = std::chrono::steady_clock::now();
      for (int i = 0; i < steps; ++i) sim.Step();
      const double ms =
          std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
      std::printf("track_energy=%d iterations=%2d: %.3f ms/step\n", track, iters, ms / steps);
    }
  }
}

// Serial versus OpenMP local step on the volumetric cantilever (strain + volume
// constraint per tet). Results are identical (see ParallelLocalStepMatchesSerial),
// so only the time differs.
TEST(PDSolverBenchmark, DISABLED_ParallelLocalStep) {
  std::printf("max threads: %d\n", PDSolver::MaxThreads());
  for (double cells : {10.0, 20.0, 40.0}) {
    for (bool parallel : {false, true}) {
      TetScene scene;
      scene.m_resolution = cells;
      PDSolver sim;
      sim.m_parallel_local_step = parallel;
      scene.Build(&sim);
      const int steps = 50;
      double local = 0.0, global = 0.0;
      const auto t0 = std::chrono::steady_clock::now();
      for (int i = 0; i < steps; ++i) {
        sim.Step();
        local += sim.GetStats().local_ms;
        global += sim.GetStats().global_ms;
      }
      const double ms =
          std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
      std::printf("%5zu tets, parallel=%d: %.3f ms/step (local %.3f, global %.3f)\n", scene.Mesh().NumTets(),
                  parallel, ms / steps, local / steps, global / steps);
    }
  }
}

// Muscle scene (strain + volume + fibre per belly tet) while the elbow animates.
TEST(PDSolverBenchmark, DISABLED_MuscleScene) {
  for (double cells : {8.0, 12.0, 16.0, 24.0}) {
    for (bool parallel : {false, true}) {
      MuscleScene scene;
      scene.m_resolution = cells;
      PDSolver sim;
      sim.m_parallel_local_step = parallel;
      scene.Build(&sim);
      const int steps = 60;
      const auto t0 = std::chrono::steady_clock::now();
      for (int i = 0; i < steps; ++i) {
        scene.Update(&sim, sim.m_h);
        sim.Step();
      }
      const double ms =
          std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
      std::printf("%5zu tets, parallel=%d: %.3f ms/step\n", scene.Mesh().NumTets(), parallel, ms / steps);
    }
  }
}

// Cost per tet of the 3x3 signed SVD and the projections, on the deformation
// gradients of the sagging cantilever (1920 tets). A strain + volume pair shares
// one SVD (PDSolver's kTetStrainVolume batch).
TEST(PDSolverBenchmark, DISABLED_TetKernels) {
  TetScene scene;
  scene.m_resolution = 20;
  PDSolver sim;
  scene.Build(&sim);
  for (int i = 0; i < 30; ++i) sim.Step();
  const TetMesh& mesh = scene.Mesh();
  std::vector<EigenMatrix3> fs;
  for (size_t t = 0; t < mesh.NumTets(); ++t) {
    ScalarType volume = 0.0;
    const Eigen::Matrix<ScalarType, 3, 4> g = TetConstraint::RestG(mesh.tets[t], mesh.vertices, &volume);
    EigenMatrix3x xt(4, 3);
    for (int c = 0; c < 4; ++c) xt.row(c) = sim.GetPositions().row(mesh.tets[t][c]);
    fs.push_back((g * xt).transpose());
  }
  const int reps = 200;
  const double n = static_cast<double>(reps * fs.size());
  ScalarType sink = 0.0;
  auto time = [&](const char* name, auto fn) {
    const auto t0 = std::chrono::steady_clock::now();
    for (int r = 0; r < reps; ++r) {
      for (const EigenMatrix3& f : fs) sink += fn(f);
    }
    const double ns = std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count();
    std::printf("%-26s %7.1f ns\n", name, ns / n);
  };
  time("SignedSVD", [](const EigenMatrix3& f) {
    EigenMatrix3 u, v;
    EigenVector3 s;
    TetConstraint::SignedSVD(f, &u, &s, &v);
    return s(2) + u(0, 0);
  });
  time("SignedSingularValues", [](const EigenMatrix3& f) { return TetConstraint::SignedSingularValues(f)(2); });
  time("ProjectStrainF", [](const EigenMatrix3& f) { return TetConstraint::ProjectStrainF(f, 1.0, 1.0)(0, 0); });
  time("ProjectVolumeF", [](const EigenMatrix3& f) { return TetConstraint::ProjectVolumeF(f, 1.0, 1.0)(0, 0); });
  time("strain + volume, one SVD", [](const EigenMatrix3& f) {
    EigenMatrix3 u, v;
    EigenVector3 s;
    TetConstraint::SignedSVD(f, &u, &s, &v);
    const EigenMatrix3 ps = u * TetConstraint::ClampStretch(s, 1.0, 1.0).asDiagonal() * v.transpose();
    const EigenMatrix3 pv = u * TetConstraint::ProjectVolumeSingularValues(s, 1.0, 1.0).asDiagonal() * v.transpose();
    return ps(0, 0) + pv(0, 0);
  });
  std::printf("(checksum %g)\n", sink);
}

namespace {

// Error versus wall-clock time of the accelerations. From one saved
// mid-simulation state every configuration takes one step; the error is the
// relative energy gap (E_K - E*) / (E_0 - E*), with E* from a 2000-iteration
// plain step, and the time is the best of 5 untracked repetitions of the step.
template <typename SceneT>
void ConvergenceVsTime(double resolution, int warmup, bool parallel) {
  SceneT scene;
  scene.m_resolution = resolution;
  PDSolver sim;
  sim.m_parallel_local_step = parallel;
  scene.Build(&sim);
  for (int i = 0; i < warmup; ++i) {
    scene.Update(&sim, sim.m_h);
    sim.Step();
  }
  scene.Update(&sim, sim.m_h);
  const EigenMatrix3x x0 = sim.GetPositions(), v0 = sim.GetVelocities();
  auto run = [&](PDAcceleration mode, int iterations, bool track) {
    sim.SetState(x0, v0);
    sim.m_acceleration = mode;
    sim.m_iterations_per_step = iterations;
    sim.m_track_energy = track;
    sim.Step();
  };
  run(PDAcceleration::kNone, 2000, true);
  const ScalarType e_star = sim.GetStats().energy.back();
  const ScalarType e0 = sim.GetStats().energy.front();
  std::printf("\n%s, %zu tets, parallel=%d\n", scene.Name(), scene.Mesh().NumTets(), parallel);

  struct Config {
    const char* name;
    PDAcceleration mode;
    ScalarType rho;
    int window;
  };
  const Config configs[] = {
      {"plain", PDAcceleration::kNone, 0.0, 0},
      {"chebyshev rho=0.90", PDAcceleration::kChebyshev, 0.90, 0},
      {"chebyshev rho=0.95", PDAcceleration::kChebyshev, 0.95, 0},
      {"chebyshev rho=0.99", PDAcceleration::kChebyshev, 0.99, 0},
      {"anderson m=3", PDAcceleration::kAnderson, 0.0, 3},
      {"anderson m=5", PDAcceleration::kAnderson, 0.0, 5},
      {"anderson m=8", PDAcceleration::kAnderson, 0.0, 8},
  };
  for (int iterations : {5, 10, 20, 40}) {
    for (const Config& c : configs) {
      sim.m_chebyshev_rho = c.rho > 0.0 ? c.rho : sim.m_chebyshev_rho;
      sim.m_anderson_window = c.window > 0 ? c.window : sim.m_anderson_window;
      run(c.mode, iterations, true);
      const ScalarType gap = (sim.GetStats().energy.back() - e_star) / (e0 - e_star);
      const int rejected = sim.GetStats().rejected_iterations;
      double best = std::numeric_limits<double>::infinity();
      for (int r = 0; r < 5; ++r) {
        run(c.mode, iterations, false);
        best = std::min(best, sim.GetStats().step_ms);
      }
      std::printf("  %2d it  %-20s gap %.2e  %7.3f ms  rejected %d\n", iterations, c.name, gap, best, rejected);
    }
  }
}

}  // namespace

TEST(PDSolverBenchmark, DISABLED_ConvergenceVsTime) {
  ConvergenceVsTime<TetScene>(20, 30, false);
  ConvergenceVsTime<MuscleScene>(12, 30, false);
  ConvergenceVsTime<MuscleScene>(12, 30, true);
}
