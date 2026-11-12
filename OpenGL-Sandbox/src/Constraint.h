#pragma once
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
// Each constraint contributes:
//   LHS (global, constant): w_i * G_i^T G_i
//   RHS (local, per step):  w_i * G_i^T p_i, with p_i the projection of the
//                           current positions onto the constraint set.

enum ConstraintType {
  AttConstraint,
  SpriConstraint,
};

class Constraint {
 public:
  Constraint(const ScalarType* stiffness, ConstraintType type)
      : m_stiffness(stiffness), m_ctype(type) {}
  virtual ~Constraint() {}

  ConstraintType GetConstraintType() const { return m_ctype; }
  ScalarType Stiffness() const { return *m_stiffness; }

  // Adds w * G^T G to the N x N system matrix.
  virtual void AddToLHS(std::vector<SparseMatrixTriplet>* triplets) const = 0;
  // Local step: projects the current positions x onto the constraint and adds
  // w * G^T p to the right-hand side b (both N x 3).
  virtual void AddToRHS(const EigenMatrix3x& x, EigenMatrix3x* b) const = 0;

 protected:
  // Points at the owning simulation's parameter so stiffness can be changed
  // at runtime (the system matrix must then be rebuilt).
  const ScalarType* m_stiffness;
  ConstraintType m_ctype;
};

// Pins a particle to a fixed point in space (e.g. a hair root).
class AttachmentConstraint : public Constraint {
 public:
  AttachmentConstraint(const ScalarType* stiffness, unsigned int p0,
                       const EigenVector3& fixed_point)
      : Constraint(stiffness, AttConstraint), m_p0(p0), m_fixed_point(fixed_point) {}

  unsigned int GetConstrainedVertexIndex() const { return m_p0; }
  EigenVector3 GetFixedPoint() const { return m_fixed_point; }
  void SetFixedPoint(const EigenVector3& target) { m_fixed_point = target; }

  void AddToLHS(std::vector<SparseMatrixTriplet>* triplets) const override {
    triplets->emplace_back(m_p0, m_p0, Stiffness());
  }

  void AddToRHS(const EigenMatrix3x& /*x*/, EigenMatrix3x* b) const override {
    b->row(m_p0) += Stiffness() * m_fixed_point.transpose();
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
      : Constraint(stiffness, SpriConstraint), m_p0(p0), m_p1(p1),
        m_rest_length(rest_length) {}

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
  void AddToRHS(const EigenMatrix3x& x, EigenMatrix3x* b) const override {
    EigenVector3 diff = (x.row(m_p0) - x.row(m_p1)).transpose();
    const ScalarType len = diff.norm();
    if (len < 1e-12) return;  // Degenerate: no defined direction, skip.
    const EigenVector3 d = (m_rest_length / len) * diff;
    const ScalarType w = Stiffness();
    b->row(m_p0) += w * d.transpose();
    b->row(m_p1) -= w * d.transpose();
  }

 private:
  unsigned int m_p0, m_p1;
  ScalarType m_rest_length;
};
