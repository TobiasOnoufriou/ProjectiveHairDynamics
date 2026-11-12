#pragma once
#include <cmath>
#include <vector>
#include "math_headers.h"

// Projective Dynamics constraints, following Liu et al. 2013
// ("Fast Simulation of Mass-Spring Systems") / Bouaziz et al. 2014.
//
// Every constraint here is separable per coordinate (its matrices are of the
// form G (x) I3), so the global system is assembled as an N x N matrix and solved
// once per coordinate (x, y, z), with positions stored as an N x 3 matrix where
// row i is particle i.
//
// Constraint i owns NumRows() rows of the stacked selection matrix J (its G_i)
// and the same rows of the projection buffer P. Its energy is
//   E_i(x) = w_i/2 * ||G_i x - p_i||^2,
// with p_i = Project(G_i x) the closest point of the constraint set. The solver
// assembles from these:
//   LHS (global, constant): sum_i w_i G_i^T G_i = J^T W J
//   RHS (local, per step):  sum_i w_i G_i^T p_i = J^T W P
// Project() only reads G_i x and writes its own rows of P, so the local step is
// race-free and can run in parallel; the gather J^T W P is one sparse product.

// Row-major N x 3 buffer: the rows of one constraint are contiguous.
using RowMatrix3x = Eigen::Matrix<ScalarType, Eigen::Dynamic, 3, Eigen::RowMajor>;
// A constraint's slice (NumRows() x 3) of the G x and P buffers.
using ConstRows = Eigen::Ref<const RowMatrix3x>;
using Rows = Eigen::Ref<RowMatrix3x>;

// Constraint types the solver flattens into per-type arrays with an inlined
// (non-virtual) projection kernel. Anything else, e.g. the drag attachment or
// a user-defined constraint, is projected through the virtual Project().
enum class ConstraintKind { kGeneric, kSpring, kTetStrain, kTetVolume, kTetFibre };

class Constraint {
 public:
  explicit Constraint(const ScalarType* stiffness) : m_stiffness(stiffness) {}
  virtual ~Constraint() {}

  ScalarType Stiffness() const { return *m_stiffness; }
  // Energy weight w. Defaults to the stiffness; volumetric constraints scale it
  // by the rest volume so the material does not depend on the mesh resolution.
  virtual ScalarType Weight() const { return Stiffness(); }

  // Number of rows of G (each row is one 3-vector of G x).
  virtual int NumRows() const = 0;
  // Appends the (unweighted) entries of G, shifted down to start at first_row.
  virtual void AddG(int first_row, std::vector<SparseMatrixTriplet>* triplets) const = 0;
  // Local step: given gx = G x, writes the projection p onto the constraint set.
  // Must only depend on gx (and the constraint's own data).
  virtual void Project(const ConstRows& gx, Rows p) const = 0;

  // Dimensionless violation used for metrics (e.g. relative spring strain).
  virtual ScalarType Strain(const ConstRows& /*gx*/) const { return 0.0; }
  // Relative volume change |det F - 1| (volumetric constraints only).
  virtual ScalarType VolumeError(const ConstRows& /*gx*/) const { return 0.0; }

  // Non-generic kinds are final classes whose data the solver copies into flat
  // arrays in RebuildSystemMatrix(); later changes to such an object are not
  // seen until the next rebuild (values behind pointers, e.g. activation, are).
  virtual ConstraintKind Kind() const { return ConstraintKind::kGeneric; }

 protected:
  // Points at a parameter owned by the scene or solver so stiffness can be
  // changed at runtime (the system matrix must then be rebuilt).
  const ScalarType* m_stiffness;
};

// Soft pin of a particle to a target point (used for mouse dragging). Hard pins
// such as hair roots use PDSolver::SetFixed instead (exact Dirichlet condition).
// G = e_p0^T, p = target.
class AttachmentConstraint : public Constraint {
 public:
  AttachmentConstraint(const ScalarType* stiffness, unsigned int p0,
                       const EigenVector3& fixed_point)
      : Constraint(stiffness), m_p0(p0), m_fixed_point(fixed_point) {}

  unsigned int GetConstrainedVertexIndex() const { return m_p0; }
  EigenVector3 GetFixedPoint() const { return m_fixed_point; }
  // Only moves p: right-hand side only.
  void SetFixedPoint(const EigenVector3& target) { m_fixed_point = target; }

  int NumRows() const override { return 1; }
  void AddG(int first_row, std::vector<SparseMatrixTriplet>* triplets) const override {
    triplets->emplace_back(first_row, m_p0, 1.0);
  }
  void Project(const ConstRows& /*gx*/, Rows p) const override {
    p.row(0) = m_fixed_point.transpose();
  }

 private:
  unsigned int m_p0;
  EigenVector3 m_fixed_point;
};

// Keeps two particles at a rest distance. Used for both stretch (i, i+1) and
// bending (i, i+2) springs along a strand. G = e_p0^T - e_p1^T.
class SpringConstraint final : public Constraint {
 public:
  SpringConstraint(const ScalarType* stiffness, unsigned int p0, unsigned int p1,
                   ScalarType rest_length)
      : Constraint(stiffness), m_p0(p0), m_p1(p1), m_rest_length(rest_length) {}

  ScalarType GetRestLength() const { return m_rest_length; }
  unsigned int GetConstrainedVertexIndex1() const { return m_p0; }
  unsigned int GetConstrainedVertexIndex2() const { return m_p1; }

  int NumRows() const override { return 1; }
  ConstraintKind Kind() const override { return ConstraintKind::kSpring; }
  void AddG(int first_row, std::vector<SparseMatrixTriplet>* triplets) const override {
    triplets->emplace_back(first_row, m_p0, 1.0);
    triplets->emplace_back(first_row, m_p1, -1.0);
  }
  void Project(const ConstRows& gx, Rows p) const override {
    p.row(0) = ProjectSpring(gx.row(0).transpose(), m_rest_length).transpose();
  }
  ScalarType Strain(const ConstRows& gx) const override {
    return SpringStrain(gx.row(0).transpose(), m_rest_length);
  }

  // Kernels shared with the solver's flattened loop. d = x0 - x1.
  // Projection: the spring vector with its length reset to the rest length,
  // so the energy is w/2 (|x0 - x1| - L)^2. A degenerate (zero-length) spring
  // has no direction; any vector of length L is then a closest point, so pick x.
  static EigenVector3 ProjectSpring(const EigenVector3& d, ScalarType rest_length) {
    const ScalarType len = d.norm();
    if (len < 1e-12) return EigenVector3(rest_length, 0.0, 0.0);
    return (rest_length / len) * d;
  }
  static ScalarType SpringStrain(const EigenVector3& d, ScalarType rest_length) {
    return std::abs(d.norm() - rest_length) / rest_length;
  }

 private:
  unsigned int m_p0, m_p1;
  ScalarType m_rest_length;
};
