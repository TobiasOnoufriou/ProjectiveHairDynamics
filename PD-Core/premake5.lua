-- PD-Core: the GL-free simulation library (solver, constraints, tet meshes,
-- colliders, asset format). Linked by the sandbox, the tests and the tools.
-- Depends only on Eigen (header-only); no OpenGL, glm, ImGui or assimp.
project "PD-Core"
  kind "StaticLib"
  language "C++"
  cppdialect "C++17"
  staticruntime "on"

  targetdir ("../bin/" .. outputdir .. "/%{prj.name}")
  objdir ("../bin-int/" .. outputdir .. "/%{prj.name}")

  files
  {
    "src/**.h",
    "src/**.cpp"
  }

  includedirs
  {
    "src",
    "../OpenGL-Core/vendor"
  }

  -- Must match every project that links PD-Core (see the header).
  forceincludes { "pd_eigen_config.h" }

  -- Only PDSolver's own OpenMP loop (m_parallel_local_step) uses threads.
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
