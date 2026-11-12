#include "Collider.h"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace {

constexpr double kPi = 3.14159265358979323846;

// Unit vector along q, or `fallback` when q is (numerically) zero.
EigenVector3 Direction(const EigenVector3& q, const EigenVector3& fallback) {
  const ScalarType len = q.norm();
  return len > 1e-12 ? EigenVector3(q / len) : fallback;
}

}  // namespace

Collider Collider::Plane(const EigenVector3& point, const EigenVector3& normal) {
  Collider c;
  c.shape = Shape::kPlane;
  c.translation = point;
  // Any rotation taking local +y to the normal.
  c.rotation = Eigen::Quaternion<ScalarType>::FromTwoVectors(EigenVector3::UnitY(), normal.normalized())
                   .toRotationMatrix();
  return c;
}

Collider Collider::Sphere(const EigenVector3& centre, ScalarType radius) {
  Collider c;
  c.shape = Shape::kSphere;
  c.translation = centre;
  c.radius = radius;
  return c;
}

Collider Collider::Capsule(const EigenVector3& a, const EigenVector3& b, ScalarType radius) {
  Collider c;
  c.shape = Shape::kCapsule;
  c.translation = 0.5 * (a + b);
  c.half_length = 0.5 * (b - a).norm();
  c.radius = radius;
  c.rotation = Eigen::Quaternion<ScalarType>::FromTwoVectors(EigenVector3::UnitX(),
                                                             Direction(b - a, EigenVector3::UnitX()))
                   .toRotationMatrix();
  return c;
}

Collider Collider::Box(const EigenVector3& centre, const EigenVector3& half_extents, const EigenMatrix3& rotation) {
  Collider c;
  c.shape = Shape::kBox;
  c.translation = centre;
  c.half_extents = half_extents;
  c.rotation = rotation;
  return c;
}

ScalarType Collider::SignedDistance(const EigenVector3& p, EigenVector3* normal) const {
  const EigenVector3 q = rotation.transpose() * (p - translation);
  EigenVector3 n_local;
  ScalarType d = 0.0;
  switch (shape) {
    case Shape::kPlane:
      n_local = EigenVector3::UnitY();
      d = q.y();
      break;
    case Shape::kSphere:
    case Shape::kCapsule: {
      // Distance to the core point (the centre, or the closest point of the
      // capsule's segment) minus the radius.
      const ScalarType hl = shape == Shape::kCapsule ? half_length : 0.0;
      const EigenVector3 r = q - EigenVector3(std::min(std::max(q.x(), -hl), hl), 0.0, 0.0);
      n_local = Direction(r, EigenVector3::UnitY());
      d = r.norm() - radius;
      break;
    }
    case Shape::kBox: {
      const EigenVector3 excess = q.cwiseAbs() - half_extents;
      // Outside: distance to the clamped point. A point on the surface up to
      // round-off takes the inside branch, which has a well-defined normal.
      const EigenVector3 r = q - q.cwiseMax(-half_extents).cwiseMin(half_extents);
      if (r.norm() > 1e-12) {
        n_local = r.normalized();
        d = r.norm();
      } else {
        // Inside: the nearest face is the one with the least negative excess.
        int k = 0;
        d = excess.maxCoeff(&k);
        n_local = EigenVector3::Zero();
        n_local[k] = q[k] < 0.0 ? -1.0 : 1.0;
      }
      break;
    }
  }
  if (normal) *normal = rotation * n_local;
  return d;
}

int DeepestPenetration(const std::vector<Collider>& colliders, const EigenVector3& p, ScalarType* depth,
                       EigenVector3* normal) {
  int best = -1;
  ScalarType best_depth = 0.0;
  for (size_t k = 0; k < colliders.size(); ++k) {
    if (!colliders[k].enabled) continue;
    EigenVector3 n;
    const ScalarType d = colliders[k].SignedDistance(p, &n);
    if (-d > best_depth) {
      best = static_cast<int>(k);
      best_depth = -d;
      if (normal) *normal = n;
    }
  }
  if (depth) *depth = best_depth;
  return best;
}

ScalarType ProjectOutOfColliders(const std::vector<Collider>& colliders, EigenVector3* p, int* deepest) {
  ScalarType first_depth = 0.0;
  if (deepest) *deepest = -1;
  // One pass resolves a single (or the deepest) penetration exactly; the extra
  // passes only matter where colliders overlap.
  constexpr int kPasses = 4;
  for (int pass = 0; pass < kPasses; ++pass) {
    ScalarType depth = 0.0;
    EigenVector3 n;
    const int k = DeepestPenetration(colliders, *p, &depth, &n);
    if (k < 0) break;
    if (pass == 0) {
      first_depth = depth;
      if (deepest) *deepest = k;
    }
    *p += depth * n;  // Onto the surface: p - d n with d = -depth.
  }
  return first_depth;
}

void TriangulateCollider(const Collider& c, ScalarType plane_size, std::vector<EigenVector3>* positions,
                         std::vector<unsigned int>* indices) {
  const unsigned int base = static_cast<unsigned int>(positions->size());
  std::vector<EigenVector3> local;
  std::vector<unsigned int> tris;
  auto quad = [&](unsigned int a, unsigned int b, unsigned int cc, unsigned int d) {
    tris.insert(tris.end(), {a, b, cc, a, cc, d});
  };
  switch (c.shape) {
    case Collider::Shape::kPlane: {
      const ScalarType h = 0.5 * plane_size;
      local = {EigenVector3(-h, 0, -h), EigenVector3(-h, 0, h), EigenVector3(h, 0, h), EigenVector3(h, 0, -h)};
      quad(0, 1, 2, 3);
      break;
    }
    case Collider::Shape::kBox:
      // Four vertices per face so the faces shade flat.
      for (int axis = 0; axis < 3; ++axis) {
        const int b = (axis + 1) % 3, cc = (axis + 2) % 3;
        for (int side = 0; side < 2; ++side) {
          const unsigned int first = static_cast<unsigned int>(local.size());
          const int corners[4][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
          for (const auto& corner : corners) {
            EigenVector3 q;
            q[axis] = (side ? 1.0 : -1.0) * c.half_extents[axis];
            q[b] = corner[0] * c.half_extents[b];
            q[cc] = corner[1] * c.half_extents[cc];
            local.push_back(q);
          }
          quad(first, first + 1, first + 2, first + 3);
        }
      }
      break;
    case Collider::Shape::kSphere:
    case Collider::Shape::kCapsule: {
      // Latitude rings around local x: the +x hemisphere shifted by
      // +half_length, the -x one by -half_length; the strip between the two
      // equator rings is the capsule's cylinder (empty for a sphere).
      constexpr int kHalfRings = 8, kSegments = 24;
      const ScalarType hl = c.shape == Collider::Shape::kCapsule ? c.half_length : 0.0;
      int rings = 0;
      for (int hemi = 0; hemi < 2; ++hemi) {
        for (int i = 0; i <= kHalfRings; ++i, ++rings) {
          const ScalarType theta = 0.5 * kPi * (hemi == 0 ? i : kHalfRings + i) / kHalfRings;
          const ScalarType offset = hemi == 0 ? hl : -hl;
          for (int j = 0; j < kSegments; ++j) {
            const ScalarType phi = 2.0 * kPi * j / kSegments;
            local.emplace_back(offset + c.radius * std::cos(theta), c.radius * std::sin(theta) * std::cos(phi),
                               c.radius * std::sin(theta) * std::sin(phi));
          }
        }
      }
      for (int r = 0; r + 1 < rings; ++r) {
        for (int j = 0; j < kSegments; ++j) {
          const int jn = (j + 1) % kSegments;
          quad(r * kSegments + j, (r + 1) * kSegments + j, (r + 1) * kSegments + jn, r * kSegments + jn);
        }
      }
      break;
    }
  }
  for (const EigenVector3& q : local) positions->push_back(c.rotation * q + c.translation);
  // Orient every triangle outwards (along the surface normal at its centroid)
  // instead of tracking windings per shape. Degenerate (pole) triangles have
  // no area and do not matter.
  for (size_t t = 0; t + 2 < tris.size(); t += 3) {
    const EigenVector3& a = (*positions)[base + tris[t]];
    const EigenVector3& b = (*positions)[base + tris[t + 1]];
    const EigenVector3& d = (*positions)[base + tris[t + 2]];
    EigenVector3 n;
    c.SignedDistance((a + b + d) / 3.0, &n);
    if ((b - a).cross(d - a).dot(n) < 0.0) std::swap(tris[t + 1], tris[t + 2]);
    for (int k = 0; k < 3; ++k) indices->push_back(base + tris[t + k]);
  }
}
