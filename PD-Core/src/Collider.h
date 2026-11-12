#pragma once
#include <vector>

#include "math_headers.h"

// A rigid, kinematic collision shape (ground plane, bone, obstacle). Colliders
// are owned by PDSolver; moving one only changes the projection targets of the
// contact slots (RHS), never the system matrix.
//
// The shape is defined in local coordinates and placed by a rigid pose
// world = rotation * local + translation:
//   kPlane:   the half-space below the local x-z plane (normal = local +y).
//   kSphere:  a ball of `radius` about the local origin.
//   kCapsule: the points within `radius` of the segment from -half_length to
//             +half_length along local x.
//   kBox:     the box |local_k| <= half_extents_k.
// No GL here, so contact is tested headless.
struct Collider {
  enum class Shape { kPlane, kSphere, kCapsule, kBox };

  Shape shape = Shape::kPlane;
  EigenMatrix3 rotation = EigenMatrix3::Identity();  // Local -> world.
  EigenVector3 translation = EigenVector3::Zero();
  ScalarType radius = 0.0;                            // kSphere, kCapsule.
  ScalarType half_length = 0.0;                       // kCapsule.
  EigenVector3 half_extents = EigenVector3::Zero();   // kBox.

  // Rigid velocity of the collider, used only by friction (the contact
  // projection only needs the pose). Velocity of a world point p:
  // linear_velocity + angular_velocity x (p - translation).
  EigenVector3 linear_velocity = EigenVector3::Zero();
  EigenVector3 angular_velocity = EigenVector3::Zero();
  ScalarType friction = 0.5;  // Coulomb coefficient mu.
  bool enabled = true;        // Disabled colliders are skipped (RHS only).
  bool visible = true;        // Drawn by the app (off when the scene draws it itself).

  static Collider Plane(const EigenVector3& point, const EigenVector3& normal);
  static Collider Sphere(const EigenVector3& centre, ScalarType radius);
  static Collider Capsule(const EigenVector3& a, const EigenVector3& b, ScalarType radius);
  static Collider Box(const EigenVector3& centre, const EigenVector3& half_extents,
                      const EigenMatrix3& rotation = EigenMatrix3::Identity());

  // Signed distance from world point p to the surface (negative inside) and
  // the outward unit normal at the closest surface point, so p - d n is that
  // point. Inside a box the closest face wins.
  ScalarType SignedDistance(const EigenVector3& p, EigenVector3* normal) const;
  // Velocity of the collider's material at world point p.
  EigenVector3 PointVelocity(const EigenVector3& p) const {
    return linear_velocity + angular_velocity.cross(p - translation);
  }
};

// Contact projection of one point: the closest point outside every enabled
// collider. For a single collider (or disjoint ones) this is exact: a point
// inside is moved to the closest point on that collider's surface, a point
// outside stays. Where colliders overlap, the deepest penetration is resolved
// first and the result re-checked a few times, which is a close approximation.
// Returns the deepest penetration depth found (0 if p was outside all of them)
// and, if one was found, writes that collider's index.
ScalarType ProjectOutOfColliders(const std::vector<Collider>& colliders, EigenVector3* p, int* deepest = nullptr);

// Deepest penetrating enabled collider at p: returns its index (or -1) and
// writes the (positive) depth and its outward normal.
int DeepestPenetration(const std::vector<Collider>& colliders, const EigenVector3& p, ScalarType* depth,
                       EigenVector3* normal);

// Triangulates a collider for drawing (counter-clockwise from outside). Planes
// become a square of side `plane_size` about their translation. Appends to
// *positions and *indices.
void TriangulateCollider(const Collider& collider, ScalarType plane_size, std::vector<EigenVector3>* positions,
                         std::vector<unsigned int>* indices);
