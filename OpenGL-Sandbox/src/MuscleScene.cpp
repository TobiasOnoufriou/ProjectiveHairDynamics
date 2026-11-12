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
constexpr ScalarType kJointGap = 0.05;       // Half the gap between the bone boxes at the elbow.
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
      m_animate(true),
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

std::vector<unsigned int> MuscleScene::AddBox(const EigenVector3& lo, const EigenVector3& hi,
                                              std::vector<EigenVector3>* positions) {
  std::vector<unsigned int> ids;
  const unsigned int base = m_mesh.NumVertices();
  for (int axis = 0; axis < 3; ++axis) {
    const int b = (axis + 1) % 3, c = (axis + 2) % 3;
    for (int side = 0; side < 2; ++side) {
      // Quad (0,0) (1,0) (1,1) (0,1) in the (b, c) plane winds counter-clockwise
      // around e_b x e_c = e_axis, so it faces +axis; reverse it for -axis.
      const int corners[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
      unsigned int quad[4];
      for (int k = 0; k < 4; ++k) {
        const int (&corner)[2] = corners[side == 1 ? k : 3 - k];
        EigenVector3 p;
        p[axis] = side == 1 ? hi[axis] : lo[axis];
        p[b] = corner[0] ? hi[b] : lo[b];
        p[c] = corner[1] ? hi[c] : lo[c];
        quad[k] = base + static_cast<unsigned int>(positions->size());
        positions->push_back(p);
        ids.push_back(quad[k]);
      }
      m_surface_indices.insert(m_surface_indices.end(), {quad[0], quad[1], quad[2], quad[0], quad[2], quad[3]});
    }
  }
  return ids;
}

void MuscleScene::Build(PDSolver* solver) {
  MakeMuscle();
  const unsigned int n = m_mesh.NumVertices();

  m_surface_indices.clear();
  for (const TetMesh::Triangle& tri : m_mesh.BoundaryTriangles()) {
    m_surface_indices.insert(m_surface_indices.end(), tri.begin(), tri.end());
  }
  m_edge_indices.clear();
  for (const auto& e : m_mesh.Edges()) {
    m_edge_indices.push_back(e[0]);
    m_edge_indices.push_back(e[1]);
  }

  // Bone boxes become extra fixed vertices after the muscle's, so they are drawn
  // with the surface and posed like any other Dirichlet vertex.
  std::vector<EigenVector3> bone_points;
  const EigenVector3 half(0.0, kBoneHalfWidth, kBoneHalfWidth);
  const std::vector<unsigned int> upper_box =
      AddBox(EigenVector3(kUpperArmStartX, 0, 0) - half, EigenVector3(kElbowX - kJointGap, 0, 0) + half,
             &bone_points);
  const std::vector<unsigned int> forearm_box =
      AddBox(EigenVector3(kElbowX + kJointGap, 0, 0) - 0.85 * half, EigenVector3(kForearmEndX, 0, 0) + 0.85 * half,
             &bone_points);

  EigenMatrix3x positions(n + bone_points.size(), 3);
  positions.topRows(n) = m_mesh.vertices;
  for (size_t k = 0; k < bone_points.size(); ++k) positions.row(n + k) = bone_points[k].transpose();
  VectorX masses(positions.rows());
  masses.head(n) = m_mesh.LumpedMasses(m_density);
  masses.tail(bone_points.size()).setOnes();  // Unused: bone vertices are fixed.

  solver->Clear();
  solver->SetVertices(positions, masses);

  auto make_bone = [&](const std::vector<unsigned int>& box, const std::vector<unsigned int>& tendon) {
    KinematicBone bone;
    bone.pivot = EigenVector3(kElbowX, 0.0, 0.0);
    bone.vertices = box;
    bone.vertices.insert(bone.vertices.end(), tendon.begin(), tendon.end());
    bone.rest.resize(bone.vertices.size(), 3);
    for (size_t k = 0; k < bone.vertices.size(); ++k) {
      bone.rest.row(k) = positions.row(bone.vertices[k]);
      solver->SetFixed(bone.vertices[k]);
    }
    return bone;
  };
  m_upper_arm = make_bone(upper_box, m_origin_vertices);
  m_forearm = make_bone(forearm_box, m_insertion_vertices);

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
  solver->Finalize();
  m_forearm.Apply(solver, ElbowAngleRadians());
}

void MuscleScene::Update(PDSolver* solver, double dt) {
  if (m_animate) {
    // Flex and extend; the biceps is active in proportion to the flexion.
    m_phase = std::fmod(m_phase + dt / std::max(m_cycle_seconds, 1e-3), 1.0);
    const ScalarType f = 0.5 * (1.0 - std::cos(2.0 * kPi * m_phase));
    m_elbow_angle_deg = f * kMaxAnimatedAngleDeg;
    m_activation = f;
  }
  m_forearm.Apply(solver, ElbowAngleRadians());
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
      {"Fibre stiffness", &m_stiffness_fibre, 1e3f, 1e8f, true, SceneParam::kSystemMatrix},
      {"Passive stiffness", &m_stiffness_strain, 1e2f, 1e8f, true, SceneParam::kSystemMatrix},
      {"Tendon stiffness", &m_stiffness_tendon, 1e4f, 1e9f, true, SceneParam::kSystemMatrix},
      {"Volume stiffness", &m_stiffness_volume, 0.0f, 1e9f, true, SceneParam::kSystemMatrix},
      {"Density (kg/m^3)", &m_density, 10.0f, 1e4f, true, SceneParam::kRebuildScene},
      {"Cells along muscle", &m_resolution, 4.0f, 48.0f, false, SceneParam::kRebuildScene},
  };
}
