#include "MuscleScene.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include "TetConstraints.h"

namespace {

constexpr double kPi = 3.14159265358979323846;
// Strain and volume ranges [1, 1]: ARAP strain and exact volume preservation.
constexpr ScalarType kOne = 1.0;

// Geometry in metres. The upper arm lies along x with its axis at y = z = 0 and
// the elbow at x = 1; y is up, so flexing the forearm moves it towards +y.
constexpr ScalarType kElbowX = 1.0;
constexpr ScalarType kBoneHalfWidth = 0.04;
constexpr ScalarType kJointGap = 0.05;       // Half the gap between the bones at the elbow.
constexpr ScalarType kForearmRadius = 0.85 * kBoneHalfWidth;
constexpr ScalarType kUpperArmStartX = -0.15;
constexpr ScalarType kForearmEndX = 1.95;
// Muscle: origin on the upper arm, insertion on the forearm past the elbow.
constexpr ScalarType kMuscleStartX = 0.1;
constexpr ScalarType kMuscleEndX = 1.15;
constexpr ScalarType kTendonRadius = 0.02;
constexpr ScalarType kBellyRadius = 0.065;
constexpr ScalarType kAttachY = kBoneHalfWidth + kTendonRadius;  // Tendons rest on the bones.
constexpr ScalarType kBellyY = 0.12;
// Tets whose centroid lies within this fraction of the length of either end are tendon.
constexpr ScalarType kTendonFraction = 0.12;
// Cells across the muscle in each transverse direction.
constexpr unsigned int kCrossCells = 4;

// Maps beam parameters q = (s, u, v) in [0, 1]^3 (s along the muscle, u and v
// across it) to the fusiform rest shape: the radius grows from the tendon
// radius at both ends to the belly radius in the middle, and the centre line
// arches from the bones up over the upper arm.
EigenVector3 SpindleMap(const EigenVector3& q) {
  const ScalarType s = q.x();
  const ScalarType sn = std::sin(kPi * s);
  const ScalarType radius = kTendonRadius + (kBellyRadius - kTendonRadius) * sn * sn;
  const ScalarType centre_y = kAttachY + (kBellyY - kAttachY) * sn;
  // Rounds the square cross-section [-1, 1]^2 towards a disc. This is the
  // elliptical grid mapping scaled by k < 1: k = 1 maps onto the disc exactly
  // but is singular at the corners, which would flatten the corner tets.
  constexpr ScalarType k = 0.75;
  const ScalarType a = 2.0 * q.y() - 1.0, b = 2.0 * q.z() - 1.0;
  const ScalarType da = a * std::sqrt(1.0 - 0.5 * k * b * b);
  const ScalarType db = b * std::sqrt(1.0 - 0.5 * k * a * a);
  return EigenVector3(kMuscleStartX + (kMuscleEndX - kMuscleStartX) * s, centre_y + radius * da, radius * db);
}

}  // namespace

EigenVector3 KinematicBone::Transform(const EigenVector3& rest_point, ScalarType angle) const {
  return pivot + Eigen::AngleAxis<ScalarType>(angle, EigenVector3::UnitZ()) * (rest_point - pivot);
}

void KinematicBone::Apply(PDSolver* solver, ScalarType angle) const {
  for (size_t k = 0; k < vertices.size(); ++k) {
    solver->SetFixedTarget(vertices[k], Transform(rest.row(k).transpose(), angle));
  }
}

MuscleScene::MuscleScene()
    : m_stiffness_strain(2e4),
      m_stiffness_tendon(1e7),
      m_stiffness_volume(1e6),
      m_stiffness_fibre(2e5),
      m_activation(0.0),
      m_max_contraction(0.4),
      m_elbow_angle_deg(0.0),
      m_cycle_seconds(4.0),
      m_density(1060.0),
      m_resolution(12.0),
      m_friction(0.3),
      m_stiffness_contact(1e5),
      m_animate(true),
      m_bone_contact(true),
      m_interior_contact(true),
      m_show_wireframe(false),
      m_phase(0.0) {}

ScalarType MuscleScene::ElbowAngleRadians() const { return m_elbow_angle_deg * kPi / 180.0; }

void MuscleScene::MakeMuscle() {
  const unsigned int nx = static_cast<unsigned int>(std::max(4.0, std::round(m_resolution)));
  // Build in parameter space [0, 1]^3, where the fibres run along s and the
  // tendons are the ends in s, then map the vertices to the spindle shape.
  m_mesh = MakeTetBeam(nx, kCrossCells, kCrossCells, EigenVector3(1.0, 1.0, 1.0));
  const TetMesh param = m_mesh;

  // Fibre field: tangent of the mapped s-lines (central difference), so the
  // fibres follow the spindle and converge into the tendons.
  m_mesh.SetFibresFromField([](const EigenVector3& q) -> EigenVector3 {
    const ScalarType eps = 1e-4;
    return SpindleMap(q + eps * EigenVector3::UnitX()) - SpindleMap(q - eps * EigenVector3::UnitX());
  });
  m_is_tendon.resize(param.NumTets());
  for (size_t t = 0; t < param.NumTets(); ++t) {
    const ScalarType s = param.Centroid(t).x();
    m_is_tendon[t] = s < kTendonFraction || s > 1.0 - kTendonFraction;
  }

  m_origin_vertices.clear();
  m_insertion_vertices.clear();
  for (unsigned int i = 0; i < param.NumVertices(); ++i) {
    const EigenVector3 q = param.vertices.row(i).transpose();
    m_mesh.vertices.row(i) = SpindleMap(q).transpose();
    if (q.x() <= 0.0) m_origin_vertices.push_back(i);
    if (q.x() >= 1.0) m_insertion_vertices.push_back(i);
  }
  m_mesh.FixOrientation();  // The map preserves orientation; this is a safeguard.
}

void MuscleScene::Build(PDSolver* solver) {
  MakeMuscle();

  m_surface_indices.clear();
  for (const TetMesh::Triangle& tri : m_mesh.BoundaryTriangles()) {
    m_surface_indices.insert(m_surface_indices.end(), tri.begin(), tri.end());
  }
  m_edge_indices.clear();
  for (const auto& e : m_mesh.Edges()) {
    m_edge_indices.push_back(e[0]);
    m_edge_indices.push_back(e[1]);
  }

  solver->Clear();
  solver->SetVertices(m_mesh.vertices, m_mesh.LumpedMasses(m_density));

  // Bones are kinematic: their tendon vertices are fixed (Dirichlet) and
  // posed with the bone.
  auto make_bone = [&](const std::vector<unsigned int>& tendon) {
    KinematicBone bone;
    bone.pivot = EigenVector3(kElbowX, 0.0, 0.0);
    bone.vertices = tendon;
    bone.rest.resize(bone.vertices.size(), 3);
    for (size_t k = 0; k < bone.vertices.size(); ++k) {
      bone.rest.row(k) = m_mesh.vertices.row(bone.vertices[k]);
      solver->SetFixed(bone.vertices[k]);
    }
    return bone;
  };
  m_upper_arm = make_bone(m_origin_vertices);
  m_forearm = make_bone(m_insertion_vertices);

  for (size_t t = 0; t < m_mesh.NumTets(); ++t) {
    const TetMesh::Tet& tet = m_mesh.tets[t];
    if (m_is_tendon[t]) {
      solver->AddConstraint(std::make_unique<TetStrainConstraint>(&m_stiffness_tendon, tet, m_mesh.vertices,
                                                                  &kOne, &kOne));
    } else {
      solver->AddConstraint(std::make_unique<TetStrainConstraint>(&m_stiffness_strain, tet, m_mesh.vertices,
                                                                  &kOne, &kOne));
      solver->AddConstraint(std::make_unique<TetFibreConstraint>(&m_stiffness_fibre, tet, m_mesh.vertices,
                                                                 m_mesh.fibres[t], &m_activation,
                                                                 &m_max_contraction));
    }
    solver->AddConstraint(
        std::make_unique<TetVolumeConstraint>(&m_stiffness_volume, tet, m_mesh.vertices, &kOne, &kOne));
  }

  // Bone colliders, also drawn as the bones: capsules along the bone axes (a
  // round bone lets the soft belly drape over it, where a box's edges fold the
  // tets). The upper arm ends kJointGap before the elbow, the forearm starts
  // kJointGap after it; they stay apart up to 120 degrees of flexion.
  // With m_interior_contact every free vertex gets a contact slot, so interior
  // vertices of crushed surface tets cannot slip into a bone either (fewer
  // inverted tets); otherwise only the boundary vertices do, which slows plain
  // PD less (each slot's proximal term w ||x - x_k||^2 lags the iteration).
  if (!m_interior_contact) {
    std::vector<unsigned int> boundary;
    for (const TetMesh::Triangle& tri : m_mesh.BoundaryTriangles()) {
      boundary.insert(boundary.end(), tri.begin(), tri.end());
    }
    solver->SetContactVertices(boundary);
  }
  solver->SetContactStiffness(&m_stiffness_contact);
  const Collider upper = Collider::Capsule(EigenVector3(kUpperArmStartX + kBoneHalfWidth, 0, 0),
                                           EigenVector3(kElbowX - kJointGap - kBoneHalfWidth, 0, 0), kBoneHalfWidth);
  const Collider fore = Collider::Capsule(EigenVector3(kElbowX + kJointGap + kForearmRadius, 0, 0),
                                          EigenVector3(kForearmEndX - kForearmRadius, 0, 0), kForearmRadius);
  m_forearm_rest_centre = fore.translation;
  m_upper_arm_collider = solver->AddCollider(upper);
  m_forearm_collider = solver->AddCollider(fore);

  solver->Finalize();
  m_posed_angle = ElbowAngleRadians();
  PoseForearm(solver, 0.0);
}

void MuscleScene::PoseForearm(PDSolver* solver, double dt) {
  const ScalarType angle = ElbowAngleRadians();
  m_forearm.Apply(solver, angle);
  for (int k : {m_upper_arm_collider, m_forearm_collider}) {
    solver->GetCollider(k).enabled = m_bone_contact;
    solver->GetCollider(k).friction = m_friction;
  }
  Collider& c = solver->GetCollider(m_forearm_collider);
  c.rotation = Eigen::AngleAxis<ScalarType>(angle, EigenVector3::UnitZ()).toRotationMatrix();
  c.translation = m_forearm.Transform(m_forearm_rest_centre, angle);
  // Rigid rotation about the elbow: omega = d angle / dt about +z, so the
  // centre moves with omega x (centre - pivot).
  c.angular_velocity = EigenVector3(0.0, 0.0, dt > 0.0 ? (angle - m_posed_angle) / dt : 0.0);
  c.linear_velocity = c.angular_velocity.cross(c.translation - m_forearm.pivot);
  m_posed_angle = angle;
}

void MuscleScene::Update(PDSolver* solver, double dt) {
  if (m_animate) {
    // Flex and extend; the biceps is active in proportion to the flexion.
    m_phase = std::fmod(m_phase + dt / std::max(m_cycle_seconds, 1e-3), 1.0);
    const ScalarType f = 0.5 * (1.0 - std::cos(2.0 * kPi * m_phase));
    m_elbow_angle_deg = f * kMaxAnimatedAngleDeg;
    m_activation = f;
  }
  PoseForearm(solver, dt);
}

const std::vector<unsigned int>& MuscleScene::LineIndices() const {
  static const std::vector<unsigned int> kNone;
  return m_show_wireframe ? m_edge_indices : kNone;
}

std::vector<SceneParam> MuscleScene::Params() {
  return {
      {"Activation", &m_activation, 0.0f, 1.0f, false, SceneParam::kNone},
      {"Max contraction", &m_max_contraction, 0.0f, 0.7f, false, SceneParam::kNone},
      {"Elbow angle (deg)", &m_elbow_angle_deg, 0.0f, 140.0f, false, SceneParam::kNone},
      {"Cycle (s)", &m_cycle_seconds, 1.0f, 10.0f, false, SceneParam::kNone},
      {"Bone friction", &m_friction, 0.0f, 2.0f, false, SceneParam::kNone},
      {"Contact stiffness", &m_stiffness_contact, 1e3f, 1e7f, true, SceneParam::kSystemMatrix},
      {"Fibre stiffness", &m_stiffness_fibre, 1e3f, 1e8f, true, SceneParam::kSystemMatrix},
      {"Passive stiffness", &m_stiffness_strain, 1e2f, 1e8f, true, SceneParam::kSystemMatrix},
      {"Tendon stiffness", &m_stiffness_tendon, 1e4f, 1e9f, true, SceneParam::kSystemMatrix},
      {"Volume stiffness", &m_stiffness_volume, 0.0f, 1e9f, true, SceneParam::kSystemMatrix},
      {"Density (kg/m^3)", &m_density, 10.0f, 1e4f, true, SceneParam::kRebuildScene},
      {"Cells along muscle", &m_resolution, 4.0f, 48.0f, false, SceneParam::kRebuildScene},
  };
}
