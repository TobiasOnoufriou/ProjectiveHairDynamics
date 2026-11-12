#pragma once
#include <memory>
#include <vector>

#include "Constraint.h"
#include "math_headers.h"

// Diagnostics of the most recent Step().
struct PDStepStats {
  double step_ms = 0.0;             // Wall-clock time of the step.
  double local_ms = 0.0;            // ... of which local steps (J x and projections).
  double global_ms = 0.0;           // ... of which gathers and global solves.
  // Total PD energy E(x_k) of the accepted iterate before iteration k and after
  // the last one (iterations + 1 values). Only filled when m_track_energy is set.
  std::vector<ScalarType> energy;
  ScalarType max_strain = 0.0;        // Max Constraint::Strain after the step.
  ScalarType max_volume_error = 0.0;  // Max Constraint::VolumeError after the step.
  int rejected_iterations = 0;        // Accelerated iterates that raised E (see PDAcceleration).
};

// How the local/global iterates are extrapolated. Every scheme is safeguarded:
// an extrapolated iterate whose energy is higher than the last accepted one is
// replaced by the plain PD iterate from that accepted one (which never raises
// E) and the acceleration restarts. So E(x_k) is non-increasing in every mode,
// and the last iteration of a step is always a plain one.
enum class PDAcceleration {
  kNone,
  // Chebyshev semi-iterative method (Wang 2015):
  //   x_{k+1} = w_{k+1} (gamma (G(x_k) - x_k) + x_k - x_{k-1}) + x_{k-1},
  //   w_1 = 1, w_2 = 2 / (2 - rho^2), w_{k+1} = 4 / (4 - rho^2 w_k),
  // with rho an estimate of the spectral radius of the plain iteration.
  kChebyshev,
  // Anderson acceleration (Peng et al. 2018) of the fixed-point map G (one
  // local + global step): with residuals f_j = G(x_j) - x_j over the last m
  // iterates, x_{k+1} = G(x_k) - dG theta, theta = argmin ||f_k - dF theta||.
  kAnderson,
};

// Generic Projective Dynamics solver (Bouaziz et al. 2014) on a set of
// vertices with per-vertex (lumped) masses and an arbitrary constraint list.
// Scene builders (e.g. HairScene) fill it via SetVertices / SetFixed /
// AddConstraint and then call Finalize().
//
// One step (h = fixed time step, M = diag(m_i)):
//   y  = x + h v + h^2 g                     (inertia + gravity, free vertices)
//   b0 = M/h^2 y - A_fc x_c                  (constant over the iterations)
//   repeat iterations:
//     local:  P = project(J x)               (each constraint its own rows of P)
//     gather: b = b0 + J^T W P               (one sparse product)
//     global: A_ff x_f = b                   (prefactored)
//     optionally extrapolate x_f             (m_acceleration)
//   v = (x_new - x_old) / h, then damped.
// J stacks every constraint's G_i, W holds their weights w_i, and
// A = M/h^2 + J^T W J. Fixed (Dirichlet) vertices x_c are eliminated from the
// system: only the free block A_ff is factored, and the coupling A_fc moves
// their known positions to the RHS. Moving a fixed vertex (e.g. a bone
// attachment) therefore only changes the RHS.
//
// Local step layout: RebuildSystemMatrix() sorts the constraints by
// ConstraintKind and copies the data of the built-in kinds into flat per-type
// arrays (ConstraintBatch), so the hot loop calls inlined kernels on contiguous
// rows instead of one virtual Project() per constraint. Each constraint
// projects into its own rows of the buffer P, so the loop is race-free and
// optionally runs in parallel (m_parallel_local_step, OpenMP); serial and
// parallel give identical positions. Step() makes no heap allocations once the
// buffers are sized.
//
// Each iteration minimises g(x, p) = 1/(2h^2) ||x - y||_M^2 + sum_i w_i/2
// ||G_i x - p_i||^2 alternately over p and x, so the energy
// E(x) = min_p g(x, p) never increases across plain iterations.
using RowSparseMatrix = Eigen::SparseMatrix<ScalarType, Eigen::RowMajor>;

class PDSolver {
 public:
  // Largest Anderson window (history kept in fixed-size small matrices).
  static constexpr int kMaxAndersonWindow = 10;

  PDSolver();
  ~PDSolver();

  // --- Scene construction ----------------------------------------------------
  // Removes all vertices and constraints.
  void Clear();
  // positions: N x 3, masses: N lumped vertex masses. Marks every vertex free.
  void SetVertices(const EigenMatrix3x& positions, const VectorX& masses);
  // Pins a vertex exactly at its current position (Dirichlet). Changes the
  // system's structure, so RebuildSystemMatrix() (or Finalize) must follow.
  void SetFixed(unsigned int vertex, bool fixed = true);
  // Built-in constraint kinds are copied into flat arrays on the next rebuild;
  // see Constraint::Kind().
  void AddConstraint(std::unique_ptr<Constraint> constraint);
  // Stores the current state as the reset pose and factors the system.
  void Finalize();

  // Moves a fixed vertex (kinematic target). RHS only, no refactorization.
  void SetFixedTarget(unsigned int vertex, const EigenVector3& target);

  // --- Time stepping ---------------------------------------------------------
  // Advances by frame time dt (seconds) using fixed sub-steps of m_h, so the
  // simulation speed does not depend on the frame rate. If steps run slower
  // than real time (and m_drop_backlog_when_slow is set), the remaining
  // backlog is dropped after the first such step: the simulation then runs in
  // slow motion at one step per frame instead of spiralling into ever longer
  // frames of m_max_substeps_per_frame steps each.
  void Update(float dt);
  // A single fixed time step.
  void Step();
  // Puts every vertex back in the pose stored by Finalize().
  void Reset();
  // Overwrites positions and velocities (N x 3 each), e.g. to replay a saved
  // state. Fixed vertices take the given positions as targets; RHS only.
  void SetState(const EigenMatrix3x& positions, const EigenMatrix3x& velocities);
  // Must be called after changing any stiffness, mass, pinning or m_h.
  void RebuildSystemMatrix();

  // Mouse dragging: a temporary soft attachment pulling one free vertex towards
  // a target point. Begin/End refactor the system matrix; moving the target only
  // changes the right-hand side, so UpdateDrag is cheap.
  void BeginDrag(unsigned int vertex, const EigenVector3& target);
  void UpdateDrag(const EigenVector3& target);
  void EndDrag();
  bool IsDragging() const { return m_drag_constraint != nullptr; }
  int DraggedVertex() const;

  // --- Queries -----------------------------------------------------------------
  unsigned int NumVertices() const { return static_cast<unsigned int>(m_positions.rows()); }
  unsigned int NumFixed() const { return static_cast<unsigned int>(m_fixed_vertices.size()); }
  size_t NumConstraints() const { return m_constraints.size(); }
  bool IsFixed(unsigned int vertex) const { return m_free_index[vertex] < 0; }
  bool IsReady() const { return m_system_ready; }
  const EigenMatrix3x& GetPositions() const { return m_positions; }
  const EigenMatrix3x& GetVelocities() const { return m_velocities; }
  const VectorX& GetMasses() const { return m_masses; }
  const PDStepStats& GetStats() const { return m_stats; }

  // The following need an up-to-date system (Finalize / RebuildSystemMatrix).
  // Sum of the constraint energies w_i/2 ||G_i x - p_i(x)||^2 (drag included).
  ScalarType ConstraintEnergy(const EigenMatrix3x& x) const;
  // Largest Constraint::Strain over all constraints at x.
  ScalarType MaxStrain(const EigenMatrix3x& x) const;
  // Largest Constraint::VolumeError over all constraints at x.
  ScalarType MaxVolumeError(const EigenMatrix3x& x) const;

  // Threads available to the parallel local step (1 without OpenMP).
  static int MaxThreads();

  // --- Tunable parameters ------------------------------------------------------
  // Changing m_h or m_stiffness_drag (while dragging) requires RebuildSystemMatrix().
  ScalarType m_h;                    // Fixed time step (s).
  ScalarType m_gravity_constant;     // Acceleration along -y (m/s^2).
  ScalarType m_stiffness_drag;
  ScalarType m_damping_coefficient;  // Fraction of velocity removed per step.
  int m_iterations_per_step;
  int m_max_substeps_per_frame;
  bool m_track_energy;               // Record PDStepStats::energy (costs one extra local pass).
  bool m_parallel_local_step;        // Project constraints (and metrics) with OpenMP threads.
  bool m_drop_backlog_when_slow;     // See Update(). Off for deterministic sub-stepping.

  PDAcceleration m_acceleration;
  ScalarType m_chebyshev_rho;        // Spectral radius estimate, in (0, 1).
  ScalarType m_chebyshev_gamma;      // Under-relaxation of the plain update, in (0, 1].
  int m_anderson_window;             // History length m, in [1, kMaxAndersonWindow].

 private:
  // A run of constraints of one kind, stored as flat arrays. Element e owns
  // rows [first_row + e * rows_per_element, ... + rows_per_element) of J, gx and
  // P; generic constraints use generic_rows instead.
  // Batch types: one per built-in ConstraintKind, plus strain + volume pairs on
  // the same tet. Both have G = Dm^-T E, hence the same F, so a pair shares one
  // SVD (the most expensive part of the local step).
  enum class BatchKind { kSpring, kTetStrainVolume, kTetStrain, kTetVolume, kTetFibre, kGeneric };
  struct ConstraintBatch {
    BatchKind kind = BatchKind::kGeneric;
    int first_row = 0;
    int rows_per_element = 0;                 // 0 for kGeneric; 6 for a pair (strain rows, then volume rows).
    std::vector<ScalarType> weight;           // w_e (the strain weight of a pair).
    std::vector<ScalarType> weight2;          // kTetStrainVolume: volume weight.
    std::vector<ScalarType> rest_length;      // kSpring.
    // Strain/volume: admissible range; kTetFibre: activation and max contraction;
    // kTetStrainVolume: strain range in param0/1, volume range in param2/3.
    std::vector<const ScalarType*> param0, param1, param2, param3;
    std::vector<EigenVector3> fibre;          // kTetFibre: unit rest fibre direction.
    std::vector<const Constraint*> generic;   // kGeneric: projected through the virtual interface.
    std::vector<int> generic_rows;            // kGeneric: first row of each element, plus the end.
    int Size() const { return static_cast<int>(weight.size()); }
  };

  // Inertial part 1/(2h^2) sum_free m_i ||x_i - y_i||^2.
  ScalarType inertiaEnergy(const EigenMatrix3x& x, const EigenMatrix3x& y) const;
  // Local step: gx = J x, then every constraint projects its rows of gx into
  // the same rows of p. Returns the summed constraint energy.
  ScalarType localStep(const EigenMatrix3x& x, RowMatrix3x* gx, RowMatrix3x* p) const;
  // Projects one batch (inside an OpenMP parallel region if one is active) and
  // returns this thread's share of its energy.
  static ScalarType projectBatch(const ConstraintBatch& batch, const RowMatrix3x& gx, RowMatrix3x* p);
  // Raises *strain / *volume_error to this thread's share of the batch's maxima.
  static void metricsBatch(const ConstraintBatch& batch, const RowMatrix3x& gx, ScalarType* strain,
                           ScalarType* volume_error);
  // Max Constraint::Strain and VolumeError at x in one pass (parallel with
  // m_parallel_local_step), using gx as scratch. Either output may be null.
  void maxMetrics(const EigenMatrix3x& x, RowMatrix3x* gx, ScalarType* strain, ScalarType* volume_error) const;
  bool laidOutFor(const EigenMatrix3x& x) const { return m_j.cols() == x.rows() && m_j.rows() > 0; }
  // x = A_ff^-1 b through the stored factor, without temporaries.
  void globalSolve(const EigenMatrix3x& b, EigenMatrix3x* x);
  // Writes the free rows u into the full position matrix x.
  void scatterFree(const EigenMatrix3x& u, EigenMatrix3x* x) const;
  // Sizes every per-step buffer for the current system.
  void allocateBuffers();
  // Anderson: next iterate from the plain iterate g and residual f = g - u.
  void andersonUpdate(const EigenMatrix3x& g, const EigenMatrix3x& f, EigenMatrix3x* u_next);

  std::vector<std::unique_ptr<Constraint>> m_constraints;
  std::unique_ptr<AttachmentConstraint> m_drag_constraint;  // Null when nothing is dragged.

  EigenMatrix3x m_positions;          // N x 3, row i = vertex i.
  EigenMatrix3x m_velocities;         // N x 3.
  EigenMatrix3x m_initial_positions;
  VectorX m_masses;                   // N lumped masses.

  // Dirichlet bookkeeping: m_free_index[i] is i's row in the reduced system,
  // or -1 if vertex i is fixed.
  std::vector<int> m_free_index;
  std::vector<unsigned int> m_free_vertices;
  std::vector<unsigned int> m_fixed_vertices;

  std::vector<ConstraintBatch> m_batches;     // Rows in batch order.
  RowSparseMatrix m_j;                        // J (rows x N).
  RowSparseMatrix m_gather;                   // Free rows of J^T W (free x rows).

  Eigen::SimplicialLLT<SparseMatrix> m_llt;  // Factor of A_ff.
  SparseMatrix m_coupling;                   // A_fc (free x fixed).

  // Per-step buffers, sized by allocateBuffers() so Step() never allocates.
  RowMatrix3x m_gx, m_projections;           // rows x 3.
  EigenMatrix3x m_x, m_y;                    // N x 3: current iterate, inertial target.
  EigenMatrix3x m_x_fixed;                   // Fixed positions.
  EigenMatrix3x m_b0, m_b, m_solve_tmp;      // Free x 3.
  EigenMatrix3x m_u, m_u_prev;               // Free iterate x_k and x_{k-1}.
  EigenMatrix3x m_g, m_g_accepted;           // Plain iterate G(x_k); G of the last accepted x_k.
  EigenMatrix3x m_f, m_f_prev, m_g_prev;     // Anderson residuals and previous plain iterate.
  std::vector<EigenMatrix3x> m_aa_df, m_aa_dg;  // Anderson ring buffers of dF and dG columns.
  int m_aa_count;                               // dF/dG pairs pushed since the last restart.
  bool m_aa_have_prev;                          // m_f_prev / m_g_prev are valid.

  bool m_system_ready;
  double m_time_accumulator;
  PDStepStats m_stats;
};
