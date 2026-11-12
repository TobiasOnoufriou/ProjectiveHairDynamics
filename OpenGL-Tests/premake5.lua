-- Unit tests for the simulation code. The solver only depends on Eigen and glm,
-- so its sources are compiled straight into this executable: no OpenGL context,
-- window or assimp is needed and the tests can run headless in CI.
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
    -- Simulation code under test.
    "../OpenGL-Sandbox/src/PDSolver.h",
    "../OpenGL-Sandbox/src/PDSolver.cpp",
    "../OpenGL-Sandbox/src/Scene.h",
    "../OpenGL-Sandbox/src/HairScene.h",
    "../OpenGL-Sandbox/src/HairScene.cpp",
    "../OpenGL-Sandbox/src/Constraint.h",
    "../OpenGL-Sandbox/src/math_headers.h",
    "../OpenGL-Sandbox/src/math_headers.cpp"
  }

  includedirs
  {
    "vendor/googletest/include",
    "vendor/googletest",
    "../OpenGL-Sandbox/src",
    "../OpenGL-Core/vendor",
    "../OpenGL-Core/%{IncludeDir.glm}"
  }

  filter "system:windows"
    systemversion "latest"

  filter "configurations:Debug"
    runtime "Debug"
    symbols "on"

  filter "configurations:Release"
    runtime "Release"
    optimize "on"
