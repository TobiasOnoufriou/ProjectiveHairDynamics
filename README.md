# Projective Hair Dynamics

A real-time **Projective Dynamics** simulator written in C++17 with OpenGL, Eigen and ImGui.

It currently simulates strand-based hair: each strand is a particle chain held by a root attachment, stretch springs
and bending springs, solved with the local/global method of Liu et al. 2013 and Bouaziz et al. 2014 using a prefactored
sparse Cholesky system. You can drag strands with the mouse and tune stiffness, damping, gravity and iteration count live.

The project is moving towards volumetric **muscle / soft-tissue** simulation (tetrahedral meshes, bone attachments,
fibre activation). See the roadmap in [CLAUDE.md](CLAUDE.md).

## Building (Windows)

```
git clone --recursive https://github.com/TobiasOnoufriou/ProjectiveHairDynamics
```

1. Run `scripts/Win-Premake.bat` to generate `OpenGL-Sandbox.sln`.
2. Open it in Visual Studio. The bundled premake (alpha14) targets VS2019 (`v142`); on VS2022, retarget the
   projects to `v143` when prompted, or build from the command line with `-p:PlatformToolset=v143`.
3. Build and run `OpenGL-Sandbox`. Its working directory must be `OpenGL-Sandbox/` (the Visual Studio default),
   because shaders and assets are loaded by relative path.

Use the Release configuration when judging performance.

## Controls

- **WASD** move the camera; hold **Left Alt** and move the mouse to look around (configurable in the Camera tab).
- **Left-click and drag** a strand to pull it.
- The **Simulation** tab has run/step/reset buttons and live parameters.

## Tests

`OpenGL-Tests` (googletest) contains headless tests for the solver. Build the `OpenGL-Tests` project and run
`bin/<config>-windows-x86_64/OpenGL-Tests/OpenGL-Tests.exe`. CI builds and runs them on every push.
