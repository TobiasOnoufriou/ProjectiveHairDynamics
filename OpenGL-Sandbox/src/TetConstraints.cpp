#include "TetConstraints.h"

#include <algorithm>
#include <cmath>
#include <vector>

TetConstraint::TetConstraint(const ScalarType* stiffness, const std::array<unsigned int, 4>& tet,
                             const EigenMatrix3x& rest)
    : Constraint(stiffness), m_tet(tet) {
  const EigenVector3 x0 = rest.row(tet[0]).transpose();
  EigenMatrix3 dm;
  for (int k = 0; k < 3; ++k) dm.col(k) = rest.row(tet[k + 1]).transpose() - x0;
  m_rest_volume = std::abs(dm.determinant()) / 6.0;

  Eigen::Matrix<ScalarType, 3, 4> e;
  e << -1, 1, 0, 0,
       -1, 0, 1, 0,
       -1, 0, 0, 1;
  m_g = dm.inverse().transpose() * e;
}

void TetConstraint::AddG(int first_row, std::vector<SparseMatrixTriplet>* triplets) const {
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 4; ++c) triplets->emplace_back(first_row + r, m_tet[c], m_g(r, c));
  }
}

void TetConstraint::SignedSVD(const EigenMatrix3& f, EigenMatrix3* u, EigenVector3* s,
                              EigenMatrix3* v) {
  Eigen::JacobiSVD<EigenMatrix3> svd(f, Eigen::ComputeFullU | Eigen::ComputeFullV);
  *u = svd.matrixU();
  *v = svd.matrixV();
  *s = svd.singularValues();  // Descending, so s(2) is the smallest.
  // Turn reflections into rotations by negating the smallest singular value.
  if (u->determinant() < 0.0) {
    u->col(2) *= -1.0;
    (*s)(2) *= -1.0;
  }
  if (v->determinant() < 0.0) {
    v->col(2) *= -1.0;
    (*s)(2) *= -1.0;
  }
}

ScalarType TetConstraint::Strain(const ConstRows& gx) const {
  Eigen::JacobiSVD<EigenMatrix3> svd(DeformationGradient(gx));
  EigenVector3 s = svd.singularValues();
  if (DeformationGradient(gx).determinant() < 0.0) s(2) = -s(2);
  return (s.array() - 1.0).abs().maxCoeff();
}

ScalarType TetConstraint::VolumeError(const ConstRows& gx) const {
  return std::abs(DeformationGradient(gx).determinant() - 1.0);
}

void TetStrainConstraint::Project(const ConstRows& gx, Rows p) const {
  EigenMatrix3 u, v;
  EigenVector3 s;
  SignedSVD(DeformationGradient(gx), &u, &s, &v);
  // Clamping each singular value independently is the closest point in
  // Frobenius norm, since the admissible set is a box in singular-value space.
  for (int k = 0; k < 3; ++k) s(k) = std::min(std::max(s(k), *m_min), *m_max);
  p = (u * s.asDiagonal() * v.transpose()).transpose();
}

void TetVolumeConstraint::Project(const ConstRows& gx, Rows p) const {
  EigenMatrix3 u, v;
  EigenVector3 s0;
  SignedSVD(DeformationGradient(gx), &u, &s0, &v);

  // Gauss-Newton on min ||d||^2 s.t. c(s0 + d) = det - v_target = 0: linearise
  // c at the current s = s0 + d_old (gradient f = (s1 s2, s0 s2, s0 s1)) and take
  // the minimum-norm d solving c(s) + f . (d - d_old) = 0.
  EigenVector3 s = s0;
  EigenVector3 d = EigenVector3::Zero();
  // Converges quadratically near the solution; far-off starts (e.g. inverted
  // tets) need more iterations, but the loop usually exits after 2-4.
  for (int it = 0; it < 30; ++it) {
    const ScalarType det = s(0) * s(1) * s(2);
    const ScalarType target = std::min(std::max(det, *m_min), *m_max);
    if (std::abs(det - target) < 1e-12) break;
    const EigenVector3 f(s(1) * s(2), s(0) * s(2), s(0) * s(1));
    const ScalarType f2 = f.squaredNorm();
    if (f2 < 1e-24) break;  // Collapsed to a line or point: no usable gradient.
    d = -((det - target) - f.dot(d)) / f2 * f;
    s = s0 + d;
  }
  p = (u * s.asDiagonal() * v.transpose()).transpose();
}
