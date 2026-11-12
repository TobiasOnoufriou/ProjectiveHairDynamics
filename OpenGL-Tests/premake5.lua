-- Unit tests for the simulation code. Links the GL-free PD-Core library and
-- compiles the sandbox's scene builders (also GL-free) directly: no OpenGL
-- context, window or assimp is needed and the tests can run headless in CI.
project "OpenGL-Tests"
  kind "ConsoleApp"
  language "C++"
  cppdialect "C++17"
  staticruntime "on"

  targetdir ("../bin/" .. outputdir .. "/%{prj.name}")
  objdir ("../bin-int/" .. outputdir .. "/%{prj.name}")

  files
  {
    "src/**.h",
    "src/**.cpp",
    "vendor/googletest/src/gtest-all.cc",
    -- GL-free scene builders from the sandbox (the solver itself is PD-Core).
    "../OpenGL-Sandbox/src/Scene.h",
    "../OpenGL-Sandbox/src/HairScene.h",
    "../OpenGL-Sandbox/src/HairScene.cpp",
    "../OpenGL-Sandbox/src/TetScene.h",
    "../OpenGL-Sandbox/src/TetScene.cpp",
    "../OpenGL-Sandbox/src/MuscleScene.h",
    "../OpenGL-Sandbox/src/MuscleScene.cpp"
  }

  includedirs
  {
    "src",
    "vendor/googletest/include",
    "vendor/googletest",
    "../PD-Core/src",
    "../OpenGL-Sandbox/src",
    "../OpenGL-Core/vendor"
  }

  links { "PD-Core" }

  -- Same Eigen configuration as PD-Core; its runtime no-malloc check is
  -- switchable per test (see the header).
  forceincludes { "pd_eigen_config.h" }

  -- Only our own OpenMP loops (PDSolver::m_parallel_local_step) may use
  -- threads, so serial-vs-parallel timings are not muddied by Eigen's own.
  defines { "EIGEN_DONT_PARALLELIZE" }

  filter "system:windows"
    systemversion "latest"
    buildoptions { "/openmp" }

  filter "configurations:Debug"
    runtime "Debug"
    symbols "on"

  filter "configurations:Release"
    runtime "Release"
    optimize "on"
