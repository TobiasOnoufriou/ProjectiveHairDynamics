# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What this project is

A real-time **Projective Dynamics (PD)** simulator in C++17 / OpenGL / Eigen / ImGui.
It started as strand-based hair simulation and is **moving to volumetric muscle / soft-tissue
simulation** (tetrahedral meshes, bone attachments, fibre activation). See the roadmap below.

The hair solver must keep working as a regression scene while the muscle work lands.

## Layout

- `OpenGL-Core/` - engine framework (window, layers, events, ImGui, camera, shaders). Rarely changed.
  Vendored deps live in `OpenGL-Core/vendor/` (GLFW, Glad, imgui, glm, eigen, assimp, spdlog). Do not edit vendor code.
- `OpenGL-Sandbox/src/` - the actual project:
  - `Simulation.h/.cpp` - PD solver (local/global loop, prefactored `SimplicialLLT`, fixed time step with accumulator).
  - `Constraint.h` - PD constraints (`AttachmentConstraint`, `SpringConstraint`). Each adds `w*G^T G` to the LHS and `w*G^T p` to the RHS.
  - `SandboxLayer.*` - scene setup, mouse picking/dragging, ImGui panels.
  - `StrandRenderer.*` - draws strands as GL lines. `Mesh.*`, `ModelLoader.*` - assimp model loading/rendering.
  - `math_headers.h` - Eigen typedefs (`ScalarType = double`, `EigenMatrix3x` = N x 3, row i = particle i).
- `OpenGL-Tests/` - googletest project. Compiles `Simulation.cpp` and `math_headers.cpp` directly (no GL context),
  so solver tests run headless. Add new solver sources to its `premake5.lua` `files` list.
- `OpenGL-Examples/` - unrelated framework examples.

## Build and run (Windows, VS2019 toolset)

```bash
scripts/Win-Premake.bat
```

Generates `OpenGL-Sandbox.sln` and the `.vcxproj` files (all gitignored). The bundled premake is alpha14, which only
knows `vs2019` (toolset `v142`). Only VS2022 is installed, so build with a toolset override:

```bash
"C:/Program Files/Microsoft Visual Studio/2022/Community/MSBuild/Current/Bin/MSBuild.exe" OpenGL-Sandbox.sln -t:OpenGL-Tests -p:Configuration=Release -p:Platform=x64 -p:PlatformToolset=v143 -m -v:minimal
```

```bash
bin/Release-windows-x86_64/OpenGL-Tests/OpenGL-Tests.exe
```

Use `-t:OpenGL-Sandbox` for the app. The user's local `.vcxproj` files may contain manual tweaks (for example a
`v143` retarget). Back them up before re-running premake. **Run the app with the working directory set to
`OpenGL-Sandbox/`**, because assets and shaders are loaded by relative paths such as `src/shaders/...`.
Use Release for any performance measurement; Eigen in Debug is 10-50x slower.

## Solver conventions (keep these)

- Positions/velocities are `EigenMatrix3x` (N x 3). Constraints must be **separable per coordinate** (`G (x) I3`)
  so the system matrix stays N x N and one factorization solves x, y and z. Tetrahedral strain constraints also fit
  this form (G maps the 4 tet vertices to the 3x3 deformation gradient).
- The system matrix `M/h^2 + sum w G^T G` depends only on masses, `h`, weights and topology. Anything that changes
  those needs `RebuildSystemMatrix()`. Anything that only moves a target (drag, bone attachments, muscle activation
  via rest shape) must change the RHS only. Never refactor per frame.
- Stiffness values are referenced by pointer (`const ScalarType*`) so ImGui sliders can change them.
- Time stepping is fixed-`h` with an accumulator in `Simulation::Update`. Do not make `h` depend on frame time.

## Code style

- Google-ish C++ style, enforced by the `cpplint` GitHub workflow. 2-space indent in new code.
- Members: new code uses `m_name`. Older files mix `name_`; do not mass-rename unrelated code.
- Comments explain the math (which energy, which projection, which paper). Cite papers by author + year.
- Prefer flat arrays/structs over `new`-allocated polymorphic objects in hot loops.

## Known issues

1. Simulation uses a single uniform particle mass. Volumetric meshes need lumped per-vertex mass (from tet volume x density).
2. Damping (`v *= 1 - c`) also damps rigid motion. Replace with damping relative to the mean velocity, or a Rayleigh-style term.
3. Root/bone pinning uses very stiff soft springs (1e6). Prefer eliminating fixed vertices from the system (Dirichlet).
4. Existing code predates the cpplint workflow and may still have style violations.

## Roadmap: hair to muscle

Work through phases in order. Each phase should leave the app building, running and testable.

### Phase 0: Repo hygiene (done)
- [x] Untrack Glad build outputs; fix Mesh slab test, include case, camera by `const&`, remove dead PD fields; README.
- [x] `OpenGL-Tests` compiles the solver headlessly with passing tests; CI builds and runs them; cpplint scoped to our sources.

### Phase 1: Generic solver
- [ ] Split `Simulation` into a reusable `PDSolver` (positions, per-vertex masses, constraint list, factorization,
      step) and scene builders (`HairScene`, later `MuscleScene`). Keep the hair demo as a selectable scene.
- [ ] Per-vertex mass (issue 8); exact Dirichlet pinning (issue 10).
- [ ] ImGui metrics: ms/step, total PD energy per iteration, max strain/volume error.
- [ ] Unit tests: rest state has zero energy; energy never increases across local/global iterations.

### Phase 2: Tetrahedral meshes
- [ ] Procedural tet beam (box split into tets). No external dependency, used for tests.
- [ ] Loader for TetGen `.node/.ele` (or Medit `.mesh`) produced offline from a surface mesh.
- [ ] Extract boundary triangles and render the surface with normals. Optionally embed a high-res render mesh with barycentric weights.

### Phase 3: Volumetric constraints
- [ ] Corotational strain constraint: F = Ds * Dm^-1, project to nearest rotation (SVD, fix reflections),
      optional singular-value clamping range (Bouaziz et al. 2014).
- [ ] Volume-preservation constraint (clamp det F towards 1).
- [ ] Demo: cantilever beam fixed at one end sagging under gravity. Tests: rotation invariance, rest-shape zero energy.

### Phase 4: Muscle model
- [ ] Per-tet fibre direction (long axis for a beam; later a vector field).
- [ ] Anisotropic fibre constraint whose target stretch along the fibre depends on an activation `a in [0,1]`
      (rest shape changes, matrix does not). ImGui activation slider.
- [ ] Tendons and bones: attach end vertices to kinematic bone boxes; animate bones (e.g. elbow flex) by moving
      attachment targets only. Demo: biceps-style muscle between two bones, contracting and bulging.

### Phase 5: Performance
- [ ] Flatten constraints into per-type arrays; fold constant RHS terms into `b0`; no per-step allocations.
- [ ] Chebyshev acceleration (Wang 2015); evaluate Anderson acceleration (Peng et al. 2018).
- [ ] Parallel local step (OpenMP). Compare error versus wall-clock time, not versus iteration count.

### Phase 6: Contact
- [ ] Ground plane and bone collisions as fixed-weight slots or post-projection (no refactor).
- [ ] Muscle-muscle / self contact (later; see Ly et al. 2020 for contact and friction in PD).

### Stretch goals
- Quasi-Newton hyperelastic materials (Liu, Bouaziz, Kavan 2017). GPU local step / Jacobi global solve.
  Real anatomical meshes. Multiple muscles with shared bones.

## References

- Liu et al. 2013, *Fast Simulation of Mass-Spring Systems*. Bouaziz et al. 2014, *Projective Dynamics*.
- Wang 2015, *A Chebyshev Semi-Iterative Approach for Accelerating Projective and Position-based Dynamics*.
- Liu, Bouaziz, Kavan 2017, *Quasi-Newton Methods for Real-time Simulation of Hyperelastic Materials*.
- Peng et al. 2018, *Anderson Acceleration for Geometry Optimization and Physics Simulation*.
- Saito, Zhou, Kavan 2015, *Computational Bodybuilding*. Ly et al. 2020, *Projective Dynamics with Dry Frictional Contact*.
