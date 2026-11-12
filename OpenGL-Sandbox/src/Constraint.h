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
// Each constraint i has the PD energy  E_i(x) = w_i/2 * ||G_i x - p_i(x)||^2,
// where p_i(x) is the projection of G_i x onto the constraint set, and contributes:
//   LHS (global, constant): w_i * G_i^T G_i
//   RHS (local, per step):  w_i * G_i^T p_i

class Constraint {
 public:
  explicit Constraint(const ScalarType* stiffness) : m_stiffness(stiffness) {}
  virtual ~Constraint() {}

  ScalarType Stiffness() const { return *m_stiffness; }

  // Adds w * G^T G to the N x N system matrix.
  virtual void AddToLHS(std::vector<SparseMatrixTriplet>* triplets) const = 0;
  // Local step: projects the current positions x onto the constraint, adds
  // w * G^T p to the right-hand side b (both N x 3) and returns the constraint
  // energy w/2 ||G x - p||^2 at x (free, since the projection is computed anyway).
  virtual ScalarType AddToRHS(const EigenMatrix3x& x, EigenMatrix3x* b) const = 0;
  // Dimensionless violation used for metrics (e.g. relative spring strain).
  virtual ScalarType Strain(const EigenMatrix3x& /*x*/) const { return 0.0; }

 protected:
  // Points at a parameter owned by the scene or solver so stiffness can be
  // changed at runtime (the system matrix must then be rebuilt).
  const ScalarType* m_stiffness;
};

// Soft pin of a particle to a target point (used for mouse dragging). Hard pins
// such as hair roots use PDSolver::SetFixed instead (exact Dirichlet condition).
class AttachmentConstraint : public Constraint {
 public:
  AttachmentConstraint(const ScalarType* stiffness, unsigned int p0,
                       const EigenVector3& fixed_point)
      : Constraint(stiffness), m_p0(p0), m_fixed_point(fixed_point) {}

  unsigned int GetConstrainedVertexIndex() const { return m_p0; }
  EigenVector3 GetFixedPoint() const { return m_fixed_point; }
  void SetFixedPoint(const EigenVector3& target) { m_fixed_point = target; }

  void AddToLHS(std::vector<SparseMatrixTriplet>* triplets) const override {
    triplets->emplace_back(m_p0, m_p0, Stiffness());
  }

  ScalarType AddToRHS(const EigenMatrix3x& x, EigenMatrix3x* b) const override {
    const ScalarType w = Stiffness();
    b->row(m_p0) += w * m_fixed_point.transpose();
    return 0.5 * w * (x.row(m_p0).transpose() - m_fixed_point).squaredNorm();
  }

 private:
  unsigned int m_p0;
  EigenVector3 m_fixed_point;
};

// Keeps two particles at a rest distance. Used for both stretch (i, i+1) and
// bending (i, i+2) springs along a strand.
class SpringConstraint : public Constraint {
 public:
  SpringConstraint(const ScalarType* stiffness, unsigned int p0, unsigned int p1,
                   ScalarType rest_length)
      : Constraint(stiffness), m_p0(p0), m_p1(p1), m_rest_length(rest_length) {}

  ScalarType GetRestLength() const { return m_rest_length; }
  unsigned int GetConstrainedVertexIndex1() const { return m_p0; }
  unsigned int GetConstrainedVertexIndex2() const { return m_p1; }

  // G = e_p0 - e_p1, so G^T G has +w on both diagonals and -w off-diagonal.
  void AddToLHS(std::vector<SparseMatrixTriplet>* triplets) const override {
    const ScalarType w = Stiffness();
    triplets->emplace_back(m_p0, m_p0, w);
    triplets->emplace_back(m_p1, m_p1, w);
    triplets->emplace_back(m_p0, m_p1, -w);
    triplets->emplace_back(m_p1, m_p0, -w);
  }

  // Projection: the spring vector with its length reset to the rest length.
  // Energy: w/2 (|x0 - x1| - L)^2.
  ScalarType AddToRHS(const EigenMatrix3x& x, EigenMatrix3x* b) const override {
    const EigenVector3 diff = (x.row(m_p0) - x.row(m_p1)).transpose();
    const ScalarType len = diff.norm();
    const ScalarType w = Stiffness();
    // Degenerate: no defined direction, so use p = 0 (adds nothing to b). The
    // returned energy is still the true distance to the constraint set, which
    // is >= the energy with p = 0, so iterations still never increase energy.
    if (len < 1e-12) return 0.5 * w * m_rest_length * m_rest_length;
    const EigenVector3 d = (m_rest_length / len) * diff;
    b->row(m_p0) += w * d.transpose();
    b->row(m_p1) -= w * d.transpose();
    return 0.5 * w * (len - m_rest_length) * (len - m_rest_length);
  }

  ScalarType Strain(const EigenMatrix3x& x) const override {
    const ScalarType len = (x.row(m_p0) - x.row(m_p1)).norm();
    return std::abs(len - m_rest_length) / m_rest_length;
  }

 private:
  unsigned int m_p0, m_p1;
  ScalarType m_rest_length;
};
