#include "TetConstraints.h"

#include <algorithm>
#include <cmath>
#include <vector>

TetConstraint::TetConstraint(const ScalarType* stiffness, const std::array<unsigned int, 4>& tet,
                             const EigenMatrix3x& rest)
    : Constraint(stiffness), m_tet(tet) {
  m_g = RestG(tet, rest, &m_rest_volume);
}

Eigen::Matrix<ScalarType, 3, 4> TetConstraint::RestG(const std::array<unsigned int, 4>& tet,
                                                     const EigenMatrix3x& rest, ScalarType* rest_volume) {
  const EigenVector3 x0 = rest.row(tet[0]).transpose();
  EigenMatrix3 dm;
  for (int k = 0; k < 3; ++k) dm.col(k) = rest.row(tet[k + 1]).transpose() - x0;
  *rest_volume = std::abs(dm.determinant()) / 6.0;

  Eigen::Matrix<ScalarType, 3, 4> e;
  e << -1, 1, 0, 0,
       -1, 0, 1, 0,
       -1, 0, 0, 1;
  return dm.inverse().transpose() * e;
}

void TetConstraint::AddG(int first_row, std::vector<SparseMatrixTriplet>* triplets) const {
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 4; ++c) triplets->emplace_back(first_row + r, m_tet[c], m_g(r, c));
  }
}

void TetConstraint::SignedSVD(const EigenMatrix3& f, EigenMatrix3* u, EigenVector3* s, EigenMatrix3* v) {
  Eigen::JacobiSVD<EigenMatrix3, Eigen::NoQRPreconditioner> svd(f, Eigen::ComputeFullU | Eigen::ComputeFullV);
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

EigenVector3 TetConstraint::SignedSingularValues(const EigenMatrix3& f) {
  Eigen::JacobiSVD<EigenMatrix3, Eigen::NoQRPreconditioner> svd(f);  // Values only.
  EigenVector3 s = svd.singularValues();
  if (f.determinant() < 0.0) s(2) = -s(2);
  return s;
}

ScalarType TetConstraint::Strain(const ConstRows& gx) const { return StrainOfF(DeformationGradient(gx)); }

ScalarType TetConstraint::VolumeError(const ConstRows& gx) const {
  return VolumeErrorOfF(DeformationGradient(gx));
}

EigenVector3 TetConstraint::ClampStretch(const EigenVector3& s, ScalarType min_stretch, ScalarType max_stretch) {
  // Clamping each singular value independently is the closest point in
  // Frobenius norm, since the admissible set is a box in singular-value space.
  return s.cwiseMax(min_stretch).cwiseMin(max_stretch);
}

EigenVector3 TetConstraint::ProjectVolumeSingularValues(const EigenVector3& s0, ScalarType min_volume,
                                                        ScalarType max_volume) {
  // Gauss-Newton on min ||d||^2 s.t. c(s0 + d) = det - v_target = 0: linearise
  // c at the current s = s0 + d_old (gradient g = (s1 s2, s0 s2, s0 s1)) and take
  // the minimum-norm d solving c(s) + g . (d - d_old) = 0.
  EigenVector3 s = s0;
  EigenVector3 d = EigenVector3::Zero();
  // Converges quadratically near the solution; far-off starts (e.g. inverted
  // tets) need more iterations, but the loop usually exits after 2-4.
  for (int it = 0; it < 30; ++it) {
    const ScalarType det = s(0) * s(1) * s(2);
    const ScalarType target = std::min(std::max(det, min_volume), max_volume);
    if (std::abs(det - target) < 1e-12) break;
    const EigenVector3 grad(s(1) * s(2), s(0) * s(2), s(0) * s(1));
    const ScalarType g2 = grad.squaredNorm();
    if (g2 < 1e-24) break;  // Collapsed to a line or point: no usable gradient.
    d = -((det - target) - grad.dot(d)) / g2 * grad;
    s = s0 + d;
  }
  return s;
}

EigenMatrix3 TetConstraint::ProjectStrainF(const EigenMatrix3& f, ScalarType min_stretch, ScalarType max_stretch) {
  EigenMatrix3 u, v;
  EigenVector3 s;
  SignedSVD(f, &u, &s, &v);
  return u * ClampStretch(s, min_stretch, max_stretch).asDiagonal() * v.transpose();
}

EigenMatrix3 TetConstraint::ProjectVolumeF(const EigenMatrix3& f, ScalarType min_volume, ScalarType max_volume) {
  EigenMatrix3 u, v;
  EigenVector3 s;
  SignedSVD(f, &u, &s, &v);
  return u * ProjectVolumeSingularValues(s, min_volume, max_volume).asDiagonal() * v.transpose();
}

void TetStrainConstraint::Project(const ConstRows& gx, Rows p) const {
  p = ProjectStrainF(DeformationGradient(gx), *m_min, *m_max).transpose();
}

void TetVolumeConstraint::Project(const ConstRows& gx, Rows p) const {
  p = ProjectVolumeF(DeformationGradient(gx), *m_min, *m_max).transpose();
}

TetFibreConstraint::TetFibreConstraint(const ScalarType* stiffness, const std::array<unsigned int, 4>& tet,
                                       const EigenMatrix3x& rest, const EigenVector3& fibre,
                                       const ScalarType* activation, const ScalarType* max_contraction)
    : Constraint(stiffness),
      m_tet(tet),
      m_fibre(fibre.normalized()),
      m_activation(activation),
      m_max_contraction(max_contraction) {
  // (F a)^T = a^T F^T = a^T (Dm^-T E X_t).
  m_g = m_fibre.transpose() * TetConstraint::RestG(tet, rest, &m_rest_volume);
}

void TetFibreConstraint::AddG(int first_row, std::vector<SparseMatrixTriplet>* triplets) const {
  for (int c = 0; c < 4; ++c) triplets->emplace_back(first_row, m_tet[c], m_g(0, c));
}

EigenVector3 TetFibreConstraint::ProjectFibre(const EigenVector3& fa, const EigenVector3& fibre,
                                              ScalarType target) {
  const ScalarType len = fa.norm();
  // Fibre collapsed to a point: every direction is equally close.
  if (len < 1e-12) return target * fibre;
  return (target / len) * fa;
}

void TetFibreConstraint::Project(const ConstRows& gx, Rows p) const {
  p.row(0) = ProjectFibre(gx.row(0).transpose(), m_fibre, TargetStretch()).transpose();
}

ScalarType TetFibreConstraint::Strain(const ConstRows& gx) const {
  return FibreStrain(gx.row(0).transpose(), TargetStretch());
}
