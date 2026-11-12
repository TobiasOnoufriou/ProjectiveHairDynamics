#pragma once
#include <array>
#include <vector>

#include "Constraint.h"
#include "math_headers.h"

// Volumetric PD constraints on one tetrahedron (Bouaziz et al. 2014, Sec. 5).
//
// With rest edge matrix Dm = [X1 - X0, X2 - X0, X3 - X0] and current edge
// matrix Ds, the deformation gradient is F = Ds Dm^-1. In the N x 3 row layout
// (row i = vertex i) the 4 vertex rows X_t of the tet give
//   Ds^T = E X_t,  E = [-1 1 0 0; -1 0 1 0; -1 0 0 1],
// so F^T = Dm^-T E X_t: G = Dm^-T E is a 3 x 4 matrix acting on the vertex
// rows, i.e. separable per coordinate, and G x = F^T. The projection p is the
// transpose of the closest admissible matrix to F, and the energy
//   w/2 ||F - P||_F^2,  w = stiffness * rest volume,
// is the discretised integral of the energy density over the tet.
class TetConstraint : public Constraint {
 public:
  // rest: N x 3 rest positions; tet must have positive rest volume.
  TetConstraint(const ScalarType* stiffness, const std::array<unsigned int, 4>& tet,
                const EigenMatrix3x& rest);

  ScalarType Weight() const override { return Stiffness() * m_rest_volume; }
  int NumRows() const override { return 3; }
  void AddG(int first_row, std::vector<SparseMatrixTriplet>* triplets) const override;

  // max_i |sigma_i - 1| over the signed singular values of F.
  ScalarType Strain(const ConstRows& gx) const override;
  // |det F - 1|: relative volume change of the tet.
  ScalarType VolumeError(const ConstRows& gx) const override;

  ScalarType RestVolume() const { return m_rest_volume; }

  // F = U diag(s) V^T with U, V proper rotations. If det F < 0 the smallest
  // singular value s(2) is negative (Irving et al. 2004), so inverted tets are
  // pushed back out instead of being projected to a reflection.
  static void SignedSVD(const EigenMatrix3& f, EigenMatrix3* u, EigenVector3* s, EigenMatrix3* v);
  // G = Dm^-T E of a tet (so G x = F^T) and its rest volume.
  static Eigen::Matrix<ScalarType, 3, 4> RestG(const std::array<unsigned int, 4>& tet, const EigenMatrix3x& rest,
                                               ScalarType* rest_volume);

 protected:
  static EigenMatrix3 DeformationGradient(const ConstRows& gx) { return gx.transpose(); }

  std::array<unsigned int, 4> m_tet;
  Eigen::Matrix<ScalarType, 3, 4> m_g;  // Dm^-T E.
  ScalarType m_rest_volume;
};

// Corotational strain constraint: F is projected onto matrices whose singular
// values lie in [*min_stretch, *max_stretch]. With the range [1, 1] this is the
// closest rotation (as-rigid-as-possible); a wider range lets the material
// stretch or compress freely by that much (Bouaziz et al. 2014, Sec. 5.1).
class TetStrainConstraint : public TetConstraint {
 public:
  TetStrainConstraint(const ScalarType* stiffness, const std::array<unsigned int, 4>& tet,
                      const EigenMatrix3x& rest, const ScalarType* min_stretch,
                      const ScalarType* max_stretch)
      : TetConstraint(stiffness, tet, rest), m_min(min_stretch), m_max(max_stretch) {}

  void Project(const ConstRows& gx, Rows p) const override;

 private:
  const ScalarType* m_min;
  const ScalarType* m_max;
};

// Volume preservation: F is projected onto matrices with det F in
// [*min_volume, *max_volume] (relative to the rest volume). The closest such
// matrix shares F's singular vectors; its singular values solve
//   min ||d - s||^2  s.t.  d0 d1 d2 = v,
// which is solved by a few Gauss-Newton (linearised constraint) iterations
// (Bouaziz et al. 2014, Sec. 5.2; as in their ShapeOp library).
class TetVolumeConstraint : public TetConstraint {
 public:
  TetVolumeConstraint(const ScalarType* stiffness, const std::array<unsigned int, 4>& tet,
                      const EigenMatrix3x& rest, const ScalarType* min_volume,
                      const ScalarType* max_volume)
      : TetConstraint(stiffness, tet, rest), m_min(min_volume), m_max(max_volume) {}

  void Project(const ConstRows& gx, Rows p) const override;

 private:
  const ScalarType* m_min;
  const ScalarType* m_max;
};

// Active muscle fibre (Saito, Zhou, Kavan 2015, "Computational Bodybuilding";
// Bouaziz et al. 2014). With the unit rest fibre direction a, F a is the
// deformed fibre and |F a| its stretch. The energy
//   w/2 ||F a - l(alpha) F a / |F a| ||^2,  w = stiffness * rest volume,
// pulls the fibre stretch towards the target l(alpha) = 1 - alpha * c_max, where
// alpha in [0, 1] is the activation and c_max the maximum contraction. This is a
// spring in material space: G = a^T Dm^-T E is one row (G x = (F a)^T), and the
// projection rescales F a to length l. Activation only moves p, so changing it
// never touches the system matrix.
class TetFibreConstraint : public Constraint {
 public:
  // fibre: rest-space direction (normalised here, must be non-zero).
  TetFibreConstraint(const ScalarType* stiffness, const std::array<unsigned int, 4>& tet,
                     const EigenMatrix3x& rest, const EigenVector3& fibre, const ScalarType* activation,
                     const ScalarType* max_contraction);

  ScalarType Weight() const override { return Stiffness() * m_rest_volume; }
  int NumRows() const override { return 1; }
  void AddG(int first_row, std::vector<SparseMatrixTriplet>* triplets) const override;
  void Project(const ConstRows& gx, Rows p) const override;
  // Relative deviation of the fibre stretch from its activated target.
  ScalarType Strain(const ConstRows& gx) const override;

  // Target fibre stretch l(alpha), activation clamped to [0, 1].
  ScalarType TargetStretch() const;
  const EigenVector3& Fibre() const { return m_fibre; }

 private:
  std::array<unsigned int, 4> m_tet;
  Eigen::Matrix<ScalarType, 1, 4> m_g;  // a^T Dm^-T E.
  EigenVector3 m_fibre;
  ScalarType m_rest_volume;
  const ScalarType* m_activation;
  const ScalarType* m_max_contraction;
};
