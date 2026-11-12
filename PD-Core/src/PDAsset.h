// Copyright 2022 Tobias Onoufriou
// The baked simulation asset (`.pdasset`): everything `pdbake` derives from a
// glTF character (tet meshes, bone attachments, fibres, materials, render-mesh
// embeddings, colliders, skin bindings), in a form the runtime can load without
// redoing any geometry processing. GL-free, no dependency beyond Eigen.
//
// Conventions (glTF's own; fixed for the whole file): meters, +Y up,
// right-handed. Everything is in the model space of the glTF scene at its bind
// pose, except where a field says bone-local. Bones are referenced by glTF node
// name, so the asset binds to any glTF/engine skeleton with the same names.
//
// File layout (all little-endian, scalars are IEEE float64 so a round trip is
// exact):
//   header  "PDASSET\0" | u32 version | u32 flags (0) | f64 meters_per_unit |
//           u8 up_axis (1 = +Y) | u8 handedness (0 = right) | u16 reserved |
//           u32 chunk_count
//   chunks  u32 tag (FourCC) | u64 payload size | payload
// Strings are u32 length + UTF-8 bytes, arrays are u32 count + items. A reader
// skips chunks with unknown tags, so later versions can add chunks without
// breaking older readers; changing an existing chunk bumps kAssetVersion and
// readers refuse files newer than they know.
#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "Collider.h"
#include "TetMesh.h"

namespace pd {

inline constexpr uint32_t kAssetVersion = 1;

// A skeleton bone the asset refers to: a glTF node name and its world (model)
// transform at the bind pose the asset was baked in. The runtime takes the
// animated transform from the glTF / engine by name; the rest transform maps
// world-space rest data (e.g. a skinned rest point) into bone space.
struct AssetBone {
  std::string name;
  EigenMatrix3 rest_rotation = EigenMatrix3::Identity();
  EigenVector3 rest_translation = EigenVector3::Zero();
};

// Per-tet material (PD constraint weights per unit rest volume and their
// projection ranges, see TetConstraints.h). Tendon tets typically have a high
// strain stiffness and no fibre (fibre_stiffness = 0).
struct AssetMaterial {
  ScalarType density = 1060.0;          // kg / m^3 (muscle tissue).
  ScalarType strain_stiffness = 1e4;    // Corotational strain.
  ScalarType min_stretch = 1.0;         // Singular-value clamp of the strain projection.
  ScalarType max_stretch = 1.0;
  ScalarType volume_stiffness = 1e5;    // Volume preservation.
  ScalarType min_volume = 1.0;          // det F range of the volume projection.
  ScalarType max_volume = 1.0;
  ScalarType fibre_stiffness = 0.0;     // Active fibre constraint; 0 = none.
  ScalarType max_contraction = 0.3;     // Fibre stretch at full activation is 1 - this.
};

// An exact Dirichlet pin: muscle vertex `vertex` follows bone `bone`, at
// `local_position` in that bone's frame (runtime: SetFixedTarget(bone_world *
// local_position), RHS only).
struct AssetAttachment {
  uint32_t vertex = 0;
  uint32_t bone = 0;  // Index into Asset::bones.
  EigenVector3 local_position = EigenVector3::Zero();
};

// One simulated muscle: its own tet mesh (rest positions, positively oriented
// tets, optional unit fibre direction per tet) plus materials and pins.
struct AssetMuscle {
  std::string name;              // Source object name (Blender / glTF node).
  std::string activation_group;  // Muscles in one group share an activation; "" = own.
  TetMesh mesh;
  std::vector<AssetMaterial> materials;
  // Material index per tet. Empty = every tet uses materials[0].
  std::vector<uint32_t> tet_material;
  // Origin-to-insertion coordinate per vertex in [0, 1] (the Laplace field
  // the fibres are the gradient of). Empty if not computed.
  std::vector<ScalarType> fibre_coordinate;
  std::vector<AssetAttachment> attachments;
};

// A point inside (or, clamped, near) a tet of a muscle: position =
// sum_k weights[k] * x[tets[tet][k]].
struct AssetTetPoint {
  uint32_t tet = 0;
  std::array<ScalarType, 4> weights = {1.0, 0.0, 0.0, 0.0};
};

// A render mesh driven by a simulated muscle: one tet point per render vertex,
// in the vertex order of the glTF mesh `mesh`.
struct AssetEmbedding {
  std::string mesh;   // glTF mesh / node name.
  uint32_t muscle = 0;
  std::vector<AssetTetPoint> vertices;
};

// A kinematic collider rigidly attached to a bone (or the world): `shape`'s
// pose (rotation / translation) is bone-local. Only the shape fields and the
// friction coefficient are stored, not the runtime state (velocities, flags).
struct AssetCollider {
  int32_t bone = -1;  // Index into Asset::bones; -1 = static in world space.
  Collider shape;
};

// Cheap skin corrective for one skin vertex (Phase 7.11): skin = LBS +
// falloff * (simulated muscle point - skinned rest muscle point).
struct AssetSkinVertex {
  int32_t muscle = -1;  // -1 = plain LBS (no muscle near this vertex).
  AssetTetPoint point;  // Nearest muscle-surface point (tet of that muscle).
  ScalarType falloff = 0.0;
};

struct AssetSkinBinding {
  std::string mesh;  // glTF skin mesh / node name.
  std::vector<AssetSkinVertex> vertices;
};

struct Asset {
  std::string generator;  // Tool and version that wrote the file.
  std::string source;     // Source file name (e.g. the .glb).
  std::vector<AssetBone> bones;
  std::vector<AssetMuscle> muscles;
  std::vector<AssetCollider> colliders;
  std::vector<AssetEmbedding> embeddings;
  std::vector<AssetSkinBinding> skin_bindings;
};

// Checks indices and sizes (tets, fibres, materials, pins, tet points, bones),
// finite values and unique bone names. Returns false and sets *error (naming
// the muscle / mesh) on the first problem. Save and Load both call it.
bool ValidateAsset(const Asset& asset, std::string* error);

// In-memory form, e.g. for engines that read from their own packages.
// SerializeAsset expects a valid asset (see ValidateAsset).
std::vector<uint8_t> SerializeAsset(const Asset& asset);
bool DeserializeAsset(const uint8_t* data, size_t size, Asset* asset, std::string* error);

bool SaveAsset(const std::string& path, const Asset& asset, std::string* error);
bool LoadAsset(const std::string& path, Asset* asset, std::string* error);

}  // namespace pd
