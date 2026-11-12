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
  - `PDSolver.h/.cpp` - generic PD solver: per-vertex lumped masses, exact Dirichlet (fixed) vertices eliminated
    from the system, constraint list stacked into one sparse `J` (A = M/h^2 + J^T W J), prefactored `SimplicialLLT`
    of the free block, race-free local step (each constraint projects into its own rows of a projection buffer,
    optional OpenMP via `m_parallel_local_step`, gather `J^T W P` is one sparse product), fixed time step with
    accumulator, per-step stats (ms/step, energy per iteration, max strain, max volume error).
  - `Constraint.h` - PD constraint interface plus `AttachmentConstraint` (soft drag) and `SpringConstraint`. A constraint
    declares `NumRows()`, its `G` entries (`AddG`), a `Weight()`, and `Project(Gx) -> p`, which may only write its
    own slice. The solver computes energies `w/2 ||Gx - p||^2` and assembles LHS/RHS.
  - `TetConstraints.*` - per-tet `G = Dm^-T E` (so `Gx = F^T`), corotational strain (clamped signed SVD) and volume
    (det F range, Gauss-Newton on singular values) constraints, weighted by stiffness * rest volume. Active fibre
    constraint: one row `G = a^T Dm^-T E` (so `Gx = (F a)^T`), projected to stretch `1 - activation * max_contraction`.
  - `TetMesh.*` - tet mesh data (positive-orientation convention), optional per-tet fibre directions, procedural Kuhn-subdivided beam, TetGen
    `.node/.ele` load/save, boundary triangles, edges, lumped masses, vertex normals. No GL, tested headless.
  - `TetScene.*` - clamped cantilever tet beam with a strain + volume constraint per tet.
  - `MuscleScene.*` - biceps-style fusiform muscle between a static upper-arm bone and a hinged forearm.
    `KinematicBone` = rigid set of fixed vertices (render box + attached tendon vertices) posed via
    `SetFixedTarget` only. Tendon tets are stiff, belly tets get a `TetFibreConstraint`. `Scene::Update` drives
    the elbow/activation animation each frame before `PDSolver::Update`.
  - `SurfaceRenderer.*` - lit triangle surface (`Scene::SurfaceIndices`), normals recomputed per upload.
  - `Scene.h`, `HairScene.*` - scene builders that fill a `PDSolver` and expose tunable `SceneParam`s (no ImGui, so
    they are testable headless). New scenes are registered in `kScenes` in `SandboxLayer.cpp`.
  - `SandboxLayer.*` - scene selection, mouse picking/dragging, ImGui panels (Simulation, Metrics, Camera).
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

1. Damping (`v *= 1 - c`) also damps rigid motion. Replace with damping relative to the mean velocity, or a Rayleigh-style term.
2. Existing code predates the cpplint workflow and may still have style violations.

## Roadmap: hair to muscle

Work through phases in order. Each phase should leave the app building, running and testable.

### Phase 0: Repo hygiene (done)
- [x] Untrack Glad build outputs; fix Mesh slab test, include case, camera by `const&`, remove dead PD fields; README.
- [x] `OpenGL-Tests` compiles the solver headlessly with passing tests; CI builds and runs them; cpplint scoped to our sources.

### Phase 1: Generic solver (done)
- [x] Split `Simulation` into a reusable `PDSolver` (positions, per-vertex masses, constraint list, factorization,
      step) and scene builders (`HairScene`, later `MuscleScene`). Keep the hair demo as a selectable scene.
- [x] Per-vertex mass; exact Dirichlet pinning (`PDSolver::SetFixed`, moved via `SetFixedTarget`, RHS only).
- [x] ImGui metrics: ms/step, total PD energy per iteration, max strain/volume error.
- [x] Unit tests: rest state has zero energy; energy never increases across local/global iterations.

### Phase 2: Tetrahedral meshes (done)
- [x] Procedural tet beam (box split into tets). No external dependency, used for tests.
- [x] Loader for TetGen `.node/.ele` produced offline from a surface mesh (Medit `.mesh` not supported yet).
- [x] Extract boundary triangles and render the surface with normals.
- [ ] Optional: embed a high-res render mesh with barycentric weights.

### Phase 3: Volumetric constraints (done)
- [x] Corotational strain constraint: F = Ds * Dm^-1, project to nearest rotation (SVD, fix reflections),
      optional singular-value clamping range (Bouaziz et al. 2014).
- [x] Volume-preservation constraint (clamp det F towards 1).
- [x] Demo: cantilever beam fixed at one end sagging under gravity. Tests: rotation invariance, rest-shape zero energy.
- [x] Race-free local step: each constraint projects into its own slot of a projection buffer, then a gather adds
      `w G^T p` into `b`. ImGui toggle for an OpenMP-parallel projection loop (`m_parallel_local_step`); compare
      ms/step on and off. Test: serial and parallel steps give identical positions.

### Phase 4: Muscle model (done)
- [x] Per-tet fibre direction (`TetMesh::fibres`; long axis for a beam, `SetFibresFromField` for a vector field;
      the muscle's fibres follow its spindle shape).
- [x] Anisotropic fibre constraint whose target stretch along the fibre depends on an activation `a in [0,1]`
      (projection changes, matrix does not). ImGui activation slider.
- [x] Tendons and bones: attach end vertices to kinematic bone boxes; animate bones (elbow flex) by moving
      attachment targets only. Demo: biceps-style muscle between two bones, contracting and bulging.
- Notes: at 120 deg flexion the muscle passes through the elbow region (no contact until Phase 6), so tets near
  the insertion are strongly stretched. Default 12 cells (1152 tets) costs ~21 ms/step serial, ~8 ms with OpenMP;
  the per-tet SVDs dominate (Phase 5).

### Phase 5: Performance
- [ ] Flatten constraints into per-type arrays; fold constant RHS terms into `b0`; no per-step allocations.
- [ ] Chebyshev acceleration (Wang 2015); evaluate Anderson acceleration (Peng et al. 2018).
- [ ] Evaluate the OpenMP toggle from Phase 3 at scale. Compare error versus wall-clock time, not versus iteration count.
- [ ] Optional GPU backend behind the same toggle: local step and a Jacobi/Chebyshev global solve both on the GPU
      (moving only the local step would be dominated by per-iteration transfers). Needs the flattened arrays above.

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
