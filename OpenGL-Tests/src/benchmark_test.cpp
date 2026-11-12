// Copyright 2022 Tobias Onoufriou
// Step-cost benchmark (disabled by default). Run with
//   OpenGL-Tests.exe --gtest_also_run_disabled_tests --gtest_filter=*Benchmark*
#include <chrono>
#include <cstdio>

#include "gtest/gtest.h"
#include "HairScene.h"
#include "PDSolver.h"

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
