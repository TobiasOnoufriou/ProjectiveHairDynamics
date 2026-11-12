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
    accumulator, per-step stats (ms/step split into local/global, energy per iteration, max strain, max volume error).
    `RebuildSystemMatrix` flattens the built-in constraint kinds (`Constraint::Kind()`) into per-type `ConstraintBatch`
    arrays with inlined kernels; strain + volume constraints on the same tet are paired and share one SVD; other
    constraints (drag, user types) go through the virtual `Project`. Optional safeguarded Chebyshev / Anderson
    acceleration (`m_acceleration`). `Step()` makes no heap allocations once the buffers are sized.
    Contact: kinematic `Collider`s (`AddCollider`, posed via `GetCollider`, RHS only) and fixed-weight contact slots
    (`SetContactVertices`, `SetContactStiffness`), plus a post-step Coulomb friction / inelastic normal velocity pass.
  - `Collider.*` - plane / sphere / capsule / box colliders: signed distance, closest-point projection out of a set of
    colliders, triangulation for drawing. No GL, tested headless.
  - `Constraint.h` - PD constraint interface plus `AttachmentConstraint` (soft drag) and `SpringConstraint`. A constraint
    declares `NumRows()`, its `G` entries (`AddG`), a `Weight()`, and `Project(Gx) -> p`, which may only write its
    own slice. The solver computes energies `w/2 ||Gx - p||^2` and assembles LHS/RHS.
  - `TetConstraints.*` - per-tet `G = Dm^-T E` (so `Gx = F^T`), corotational strain (clamped signed SVD) and volume
    (det F range, Gauss-Newton on singular values) constraints, weighted by stiffness * rest volume. Active fibre
    constraint: one row `G = a^T Dm^-T E` (so `Gx = (F a)^T`), projected to stretch `1 - activation * max_contraction`.
  - `TetMesh.*` - tet mesh data (positive-orientation convention), optional per-tet fibre directions, procedural Kuhn-subdivided beam, TetGen
    `.node/.ele` load/save, boundary triangles, edges, lumped masses, vertex normals. No GL, tested headless.
  - `TetScene.*` - clamped cantilever tet beam with a strain + volume constraint per tet, a sphere and a ground plane.
  - `MuscleScene.*` - biceps-style fusiform muscle between a static upper-arm bone and a hinged forearm.
    `KinematicBone` = rigid set of fixed tendon vertices posed via `SetFixedTarget` only; the bones themselves are
    capsule colliders (drawn by the app). Tendon tets are stiff, belly tets get a `TetFibreConstraint`.
    `Scene::Update` drives the elbow/activation animation and the forearm collider each frame before
    `PDSolver::Update`.
  - `SurfaceRenderer.*` - lit triangle surface (`Scene::SurfaceIndices`, and the visible colliders), normals
    recomputed per upload.
  - `Scene.h`, `HairScene.*` - scene builders that fill a `PDSolver` and expose tunable `SceneParam`s (no ImGui, so
    they are testable headless). New scenes are registered in `kScenes` in `SandboxLayer.cpp`.
  - `SandboxLayer.*` - scene selection, mouse picking/dragging, ImGui panels (Simulation, Metrics, Camera).
  - `StrandRenderer.*` - draws strands as GL lines. `Mesh.*`, `ModelLoader.*` - assimp model loading/rendering.
  - `math_headers.h` - Eigen typedefs (`ScalarType = double`, `EigenMatrix3x` = N x 3, row i = particle i).
- `OpenGL-Tests/` - googletest project. Compiles the solver sources directly (no GL context), so solver tests run
  headless. Add new solver sources to its `premake5.lua` `files` list. `src/eigen_test_config.h` is force-included in
  every test source: it defines `EIGEN_RUNTIME_NO_MALLOC` and an `eigen_assert` that only fires while
  `pd_test::g_eigen_asserts` is set (used by `StepDoesNotAllocate`). Benchmarks are `DISABLED_` tests in
  `benchmark_test.cpp`; run them with `--gtest_also_run_disabled_tests --gtest_filter=*Benchmark*`.
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

### Phase 5: Performance (done, except the optional GPU backend)
- [x] Flatten constraints into per-type arrays; fold constant RHS terms into `b0`; no per-step allocations
      (test: `StepDoesNotAllocate`). Strain + volume pairs on one tet share their SVD.
- [x] Chebyshev acceleration (Wang 2015) and Anderson acceleration (Peng et al. 2018), both safeguarded so the
      accepted energy never rises; the last iteration of a step is always plain.
- [x] Evaluate the OpenMP toggle at scale and error versus wall-clock time (`ConvergenceVsTime` benchmark).
- [ ] Optional GPU backend behind the same toggle: local step and a Jacobi/Chebyshev global solve both on the GPU
      (moving only the local step would be dominated by per-iteration transfers). Needs the flattened arrays above.
- Notes (Release, 20 threads; this machine's timings vary by about 10%):
  - ms/step serial / OpenMP, Phase 4 -> Phase 5: cantilever 1920 tets 38.5 / 11.4 -> 20.0 / 7.3; 15360 tets
    373 / 110 -> 172 / 75; muscle 1152 tets 25.5 / 6.4 -> 12.9 / 4.7. At 15360 tets with OpenMP the serial global
    solve (~35 ms) is half the step; it also gets ~20% slower while OpenMP threads spin after the local step.
  - 3x3 SVD alternatives measured on simulation data (~600 ns each): Eigen JacobiSVD, a one-sided cyclic Jacobi on
    F^T F (warm-started or not: near rest V is ill-conditioned, so warm starts barely help) and Eigen's closed-form
    `computeDirect` (fast but only ~sqrt(eps) accurate near repeated singular values). Kept JacobiSVD; the win
    came from flattening and sharing the SVD between paired constraints.
  - Accelerations cost no measurable time per iteration. Anderson (m = 5) reaches with 10 iterations the error
    plain PD has after 40 (muscle) or well beyond 40 (cantilever), so the app defaults to it; the solver default
    stays plain. Chebyshev with rho = 0.9 never triggers the safeguard and gives 2-3x lower error; larger rho
    overshoots and gets rejected.

### Phase 6: Contact (ground and bones done; self contact open)
- [x] Ground plane and bone collisions as fixed-weight slots (no refactor). Once a collider exists, every contact
      vertex owns one row `G = e_i` with a constant weight; its projection is the closest point outside all colliders
      (itself when outside), so posing colliders is RHS only and the energy stays monotone. After the step, contact
      vertices inside a collider lose their approaching normal velocity and get Coulomb friction with normal force
      `w * depth`, relative to the collider's rigid velocity. Demos: hair over a sphere onto the ground, a soft beam
      draped over a sphere (lower its strain stiffness), the biceps on capsule bones. Tests: `contact_test.cpp`.
- [ ] Muscle-muscle / self contact (later; see Ly et al. 2020 for contact and friction in PD).
- Notes:
  - Slots are vertex-only: edges and faces can still cut through a collider between vertices (thin colliders,
    coarse meshes). Overlapping colliders (e.g. the elbow capsules beyond 120 deg) are resolved approximately.
  - A slot whose vertex is outside still adds `w ||x - x_k||^2` to every global step, a proximal term that slows
    convergence of the slotted vertices. The weight is per scene (`Contact stiffness`, default 1e5): stiffer gives
    less penetration but slower plain PD. Anderson absorbs most of the lag.
  - Muscle: the passive belly is soft, so under gravity it slumps onto the upper arm and drapes around it; with
    slots on interior vertices too (default) and Anderson, volume stays within ~2% and boundary penetration ~6 mm at
    full activation and 120 deg (without contact the belly sinks ~4 cm into the bone). Plain PD at 10 iterations
    with contact loses up to ~15% volume and inverts some crushed tets; the Phase 4 muscle tests therefore run
    with bone contact off. Box bones folded the tets over their edges; capsules drape better.

### Phase 7: Asset pipeline (Blender -> baked asset -> running sim)
Goal: an arm modelled in Blender, exported as `.glb`, baked by `pdbake` and animated in the sandbox with simulated
muscles and a skin corrective, with no manual steps beyond tagging objects. Work it in chunks, one session each,
committing and updating this file at the end of each: **7a** (7.1, 7.10 format) first; then **7b** (7.2, 7.3, 7.13),
**7c** (7.4) and **7d** (7.5-7.9), which are independent; then **7e** (7.10 tool, 7.11, 7.12).

Authoring conventions (Blender side): objects carry a custom property `pd_role = "muscle" | "bone" | "skin"`
(exported as glTF `extras`); muscles may override `pd_stiffness`, `pd_max_contraction`, `pd_resolution`. Optional
vertex groups on muscles: `pd_origin` / `pd_insertion` (attachment regions), `pd_tendon` (0-1 tendon weight),
`pd_activation_group`. Bone meshes are rigidly parented to armature bones. Apply scale and modifiers before export.
Everything is converted to one internal convention at import: meters, Y-up, right-handed (glTF's own).

- [ ] 7.1 Core library split: solver, mesh, constraint and collider code into a GL-free `PD-Core` static library
      linked by the sandbox, the tests and `pdbake`. Keep its public headers small (future SDK C API base).
- [ ] 7.2 glTF import: meshes, armature, skin weights, animations, `extras`. Prefer cgltf (MIT, single header) over
      the vendored assimp unless assimp exposes `extras` cleanly.
- [ ] 7.3 Validation: watertight, no self-intersection, no muscle-bone intersection at rest, scale sanity, degenerate
      triangles. Errors name the Blender object and suggest a fix ("Apply Scale", "mesh has 3 holes").
- [ ] 7.4 Tetrahedralizer: surface -> tets with a resolution parameter, by isosurface stuffing on a BCC lattice
      (Labelle and Shewchuk 2007): bounded dihedral angles, no dependency, no licence issue (TetGen is AGPL; keep the
      `.node/.ele` loader but never ship TetGen). fTetWild (MPL-2.0) only as an optional fallback for bad meshes.
      Tests: tet volume ~ surface volume, all tets positive, minimum dihedral angle bound.
- [ ] 7.5 Render-mesh embedding: each render vertex -> (tet, 4 barycentric weights); vertices outside use the nearest
      tet with clamped weights. Runtime deform + normals. Tests: zero displacement at rest, rigid tets -> rigid mesh.
- [ ] 7.6 Attachments: auto pins by distance to bone meshes, vertex groups override. Stored as bone name + bone-local
      position; runtime `SetFixedTarget(bone_world * local)` (RHS only).
- [ ] 7.7 Fibre field: Laplace solve on the tets (0 at origin, 1 at insertion), fibre = normalised gradient per tet
      (Choi and Blemker 2013; Saito et al. 2015). Test: matches the hand-built spindle fibres.
- [ ] 7.8 Materials: tendon vs belly stiffness per tet from `pd_tendon`, else from the fibre-field coordinate.
- [ ] 7.9 Colliders: fit a capsule per bone mesh (principal axis + radius), attached to its armature bone.
- [ ] 7.10 `.pdasset`: versioned binary (header with version, units, axes; tets, rest positions, attachments, fibres,
      materials, embedding, colliders, skin bindings), references glTF bones by name; round-trip test.
      `pdbake input.glb -o out.pdasset [--resolution]` CLI with a report (tet count, pins per bone, warnings).
- [ ] 7.11 Cheap skin corrective: bake each skin vertex to its nearest muscle-surface point plus a falloff; runtime
      skin = LBS position + falloff * (simulated muscle point - skinned rest muscle point). Test: a muscle at rest
      relative to its bones leaves the skin equal to LBS. (A simulated fat layer is out of scope; see stretch goals.)
- [ ] 7.12 Sandbox `AssetScene`: loads `.glb` + `.pdasset`, plays the glTF animation, poses bones/colliders, runs
      the sim, draws bones, muscles (visible-muscle mode) and skin (cheap-skin mode) with toggles. Activation per
      muscle from a slider, an animated custom property or a joint-angle curve. Several muscles in one solver sharing
      bones (biceps + triceps); without muscle-muscle contact they may overlap.
- [ ] 7.13 Test asset: a Blender Python script that builds the arm (bones, biceps, triceps, skin, curl animation),
      plus the exported `.glb` checked in so CI never needs Blender. If Blender is unavailable, generate the `.glb`
      from the procedural muscle code instead.

### Phase 8: Game budget, robustness, comparison
Goal: <= 1 ms per hero character per frame on 4 threads at 60 fps (visible muscle + cheap skin), no blow-ups on
gameplay animation, and side-by-side footage against existing techniques. Starting point (Phase 5): ~11 us per tet
per frame serial at h = 1/60, 10 iterations, double precision; the per-tet SVD dominates.

- [ ] 8.1 `GameBudget` benchmark on the Phase 7 arm asset (ms/frame vs tets, iterations, threads). Baseline first;
      every later item is measured against it.
- [ ] 8.2 `ScalarType = float` builds and passes the tests. Watch the LLT with stiff contact weights (1e5); fallback:
      factorize in double, everything else in float.
- [ ] 8.3 SIMD 3x3 SVD (McAdams et al. 2011), accuracy vs JacobiSVD on simulation data, including near-repeated
      singular values (where `computeDirect` failed).
- [ ] 8.4 Replace OpenMP with a host-supplied `parallel_for` callback (serial / built-in pool fallback) so Unreal's
      task graph, Unity jobs and Godot's WorkerThreadPool can drive the solver.
- [ ] 8.5 Quality presets: tets per muscle (200-400), iterations (5-10 with Anderson), contact settings; bake-time
      warning when an asset exceeds the budget.
- [ ] 8.6 LOD: full sim / reduced rate (30 Hz, interpolated) / off (plain LBS), with blending and a clean wake-up
      (reset to the skinned pose, zero relative velocity).
- [ ] 8.7 Robustness: teleport / large-jump reset, NaN and extreme-strain guards with per-frame fallback to LBS.
      Fix known issue 1 (damping must be relative to bone motion or a running character's muscles lag).
- [ ] 8.8 Art controls per muscle: corrective strength, jiggle amount, activation scale, max bulge clamp (RHS or
      post-process only, never the matrix).
- [ ] 8.9 Comparison scene: same arm as LBS, LBS + jiggle bone, LBS + PD muscles + cheap skin; clips: slow curl,
      fast punch, holding a weight (same pose, different activation), landing. Capture video and timings.
- Positioning: complement existing skinning/correctives (keep them for joint folds and art-directed poses); PD adds
  activation-driven bulge, secondary motion, volume preservation and contact, which pose-driven methods cannot.

### Stretch goals
- Quasi-Newton hyperelastic materials (Liu, Bouaziz, Kavan 2017). GPU local step / Jacobi global solve.
  Real anatomical meshes.
- Engine SDK (Godot, Unity, Unreal) over a plain C API on `PD-Core`, loading `.pdasset`.
- Simulated fat layer (tets between skin, muscles and bones): ~1-3 ms GPU per hero character, needs self contact.
- ML-trained deformer: simulate offline with this solver, train pose (+ activation, velocity) -> PCA skin
  correction, evaluate at ~0.1-0.5 ms at runtime (as Ziva Real-Time and Unreal's ML Deformer do).

## References

- Liu et al. 2013, *Fast Simulation of Mass-Spring Systems*. Bouaziz et al. 2014, *Projective Dynamics*.
- Wang 2015, *A Chebyshev Semi-Iterative Approach for Accelerating Projective and Position-based Dynamics*.
- Liu, Bouaziz, Kavan 2017, *Quasi-Newton Methods for Real-time Simulation of Hyperelastic Materials*.
- Peng et al. 2018, *Anderson Acceleration for Geometry Optimization and Physics Simulation*.
- Saito, Zhou, Kavan 2015, *Computational Bodybuilding*. Ly et al. 2020, *Projective Dynamics with Dry Frictional Contact*.
