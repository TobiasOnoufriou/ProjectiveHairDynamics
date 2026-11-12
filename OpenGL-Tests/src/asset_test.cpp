// Copyright 2022 Tobias Onoufriou
// Tests for the baked asset format (.pdasset): exact round trip, forward
// compatibility (unknown chunks are skipped), and rejection of corrupt,
// truncated, newer or inconsistent files.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

#include "PDAsset.h"
#include "gtest/gtest.h"

namespace {

// Offset of the chunk count in the header (magic, version, flags,
// meters_per_unit, up axis, handedness, reserved).
constexpr size_t kChunkCountOffset = 8 + 4 + 4 + 8 + 1 + 1 + 2;

// An asset that fills every field with distinct, non-default values.
pd::Asset MakeTestAsset() {
  pd::Asset asset;
  asset.generator = "asset_test";
  asset.source = "arm.glb";

  pd::AssetBone upper;
  upper.name = "upper_arm";
  upper.rest_translation = EigenVector3(0.0, 1.4, 0.0);
  pd::AssetBone fore;
  fore.name = "forearm";
  fore.rest_rotation = Eigen::AngleAxis<ScalarType>(0.3, EigenVector3(0.0, 0.0, 1.0)).toRotationMatrix();
  fore.rest_translation = EigenVector3(0.31, 1.1, -0.02);
  asset.bones = {upper, fore};

  pd::AssetMuscle muscle;
  muscle.name = "biceps";
  muscle.activation_group = "flexors";
  muscle.mesh = MakeTetBeam(3, 1, 1, EigenVector3(0.3, 0.05, 0.04), EigenVector3(0.01, 1.2, -0.02));
  pd::AssetMaterial tendon;
  tendon.strain_stiffness = 1e6;
  tendon.fibre_stiffness = 0.0;
  pd::AssetMaterial belly;
  belly.density = 1070.5;
  belly.strain_stiffness = 2e3;
  belly.min_stretch = 0.9;
  belly.max_stretch = 1.1;
  belly.min_volume = 0.98;
  belly.max_volume = 1.02;
  belly.fibre_stiffness = 5e4;
  belly.max_contraction = 0.25;
  muscle.materials = {tendon, belly};
  for (size_t t = 0; t < muscle.mesh.NumTets(); ++t) muscle.tet_material.push_back(t % 3 == 0 ? 0 : 1);
  for (unsigned int v = 0; v < muscle.mesh.NumVertices(); ++v)
    muscle.fibre_coordinate.push_back(v / static_cast<ScalarType>(muscle.mesh.NumVertices() - 1));
  muscle.attachments.push_back({0, 0, EigenVector3(0.01, -0.2, -0.02)});
  muscle.attachments.push_back({muscle.mesh.NumVertices() - 1, 1, EigenVector3(0.1 / 3.0, 0.0, 0.5)});
  asset.muscles.push_back(muscle);

  pd::AssetMuscle tiny;  // A second muscle without optional data.
  tiny.name = "triceps";
  tiny.mesh = MakeTetBeam(1, 1, 1, EigenVector3(0.1, 0.1, 0.1));
  tiny.mesh.fibres.clear();
  tiny.materials = {pd::AssetMaterial()};
  asset.muscles.push_back(tiny);

  pd::AssetCollider c;
  c.bone = 0;
  c.shape = Collider::Capsule(EigenVector3(0.0, 0.0, 0.0), EigenVector3(0.0, -0.3, 0.0), 0.035);
  c.shape.friction = 0.7;
  asset.colliders.push_back(c);
  c.bone = 1;
  c.shape = Collider::Box(EigenVector3(0.1, 0.0, 0.0), EigenVector3(0.12, 0.02, 0.03),
                          fore.rest_rotation);
  asset.colliders.push_back(c);
  c.bone = -1;
  c.shape = Collider::Plane(EigenVector3::Zero(), EigenVector3(0.0, 1.0, 0.0));
  asset.colliders.push_back(c);
  c.shape = Collider::Sphere(EigenVector3(1.0, 2.0, 3.0), 0.5);
  asset.colliders.push_back(c);

  pd::AssetEmbedding embedding;
  embedding.mesh = "biceps_render";
  embedding.muscle = 0;
  embedding.vertices.push_back({3, {0.1, 0.2, 0.3, 0.4}});
  embedding.vertices.push_back({0, {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.0}});
  asset.embeddings.push_back(embedding);

  pd::AssetSkinBinding skin;
  skin.mesh = "arm_skin";
  skin.vertices.push_back({-1, {}, 0.0});
  skin.vertices.push_back({0, {5, {0.25, 0.25, 0.25, 0.25}}, 0.8});
  skin.vertices.push_back({1, {2, {0.0, 0.0, 1.0, 0.0}}, 0.125});
  asset.skin_bindings.push_back(skin);
  return asset;
}

void ExpectSameTetPoint(const pd::AssetTetPoint& a, const pd::AssetTetPoint& b) {
  EXPECT_EQ(a.tet, b.tet);
  EXPECT_EQ(a.weights, b.weights);
}

void ExpectSameCollider(const Collider& a, const Collider& b) {
  EXPECT_EQ(a.shape, b.shape);
  EXPECT_EQ(a.rotation, b.rotation);
  EXPECT_EQ(a.translation, b.translation);
  EXPECT_EQ(a.radius, b.radius);
  EXPECT_EQ(a.half_length, b.half_length);
  EXPECT_EQ(a.half_extents, b.half_extents);
  EXPECT_EQ(a.friction, b.friction);
}

// Field-by-field, exact (scalars are stored as float64).
void ExpectSameAsset(const pd::Asset& a, const pd::Asset& b) {
  EXPECT_EQ(a.generator, b.generator);
  EXPECT_EQ(a.source, b.source);
  ASSERT_EQ(a.bones.size(), b.bones.size());
  for (size_t i = 0; i < a.bones.size(); ++i) {
    EXPECT_EQ(a.bones[i].name, b.bones[i].name);
    EXPECT_EQ(a.bones[i].rest_rotation, b.bones[i].rest_rotation);
    EXPECT_EQ(a.bones[i].rest_translation, b.bones[i].rest_translation);
  }
  ASSERT_EQ(a.muscles.size(), b.muscles.size());
  for (size_t i = 0; i < a.muscles.size(); ++i) {
    const pd::AssetMuscle& ma = a.muscles[i];
    const pd::AssetMuscle& mb = b.muscles[i];
    EXPECT_EQ(ma.name, mb.name);
    EXPECT_EQ(ma.activation_group, mb.activation_group);
    EXPECT_EQ(ma.mesh.vertices, mb.mesh.vertices);
    EXPECT_EQ(ma.mesh.tets, mb.mesh.tets);
    EXPECT_EQ(ma.mesh.fibres, mb.mesh.fibres);
    ASSERT_EQ(ma.materials.size(), mb.materials.size());
    for (size_t k = 0; k < ma.materials.size(); ++k) {
      const pd::AssetMaterial& x = ma.materials[k];
      const pd::AssetMaterial& y = mb.materials[k];
      EXPECT_EQ(x.density, y.density);
      EXPECT_EQ(x.strain_stiffness, y.strain_stiffness);
      EXPECT_EQ(x.min_stretch, y.min_stretch);
      EXPECT_EQ(x.max_stretch, y.max_stretch);
      EXPECT_EQ(x.volume_stiffness, y.volume_stiffness);
      EXPECT_EQ(x.min_volume, y.min_volume);
      EXPECT_EQ(x.max_volume, y.max_volume);
      EXPECT_EQ(x.fibre_stiffness, y.fibre_stiffness);
      EXPECT_EQ(x.max_contraction, y.max_contraction);
    }
    EXPECT_EQ(ma.tet_material, mb.tet_material);
    EXPECT_EQ(ma.fibre_coordinate, mb.fibre_coordinate);
    ASSERT_EQ(ma.attachments.size(), mb.attachments.size());
    for (size_t k = 0; k < ma.attachments.size(); ++k) {
      EXPECT_EQ(ma.attachments[k].vertex, mb.attachments[k].vertex);
      EXPECT_EQ(ma.attachments[k].bone, mb.attachments[k].bone);
      EXPECT_EQ(ma.attachments[k].local_position, mb.attachments[k].local_position);
    }
  }
  ASSERT_EQ(a.colliders.size(), b.colliders.size());
  for (size_t i = 0; i < a.colliders.size(); ++i) {
    EXPECT_EQ(a.colliders[i].bone, b.colliders[i].bone);
    ExpectSameCollider(a.colliders[i].shape, b.colliders[i].shape);
  }
  ASSERT_EQ(a.embeddings.size(), b.embeddings.size());
  for (size_t i = 0; i < a.embeddings.size(); ++i) {
    EXPECT_EQ(a.embeddings[i].mesh, b.embeddings[i].mesh);
    EXPECT_EQ(a.embeddings[i].muscle, b.embeddings[i].muscle);
    ASSERT_EQ(a.embeddings[i].vertices.size(), b.embeddings[i].vertices.size());
    for (size_t k = 0; k < a.embeddings[i].vertices.size(); ++k)
      ExpectSameTetPoint(a.embeddings[i].vertices[k], b.embeddings[i].vertices[k]);
  }
  ASSERT_EQ(a.skin_bindings.size(), b.skin_bindings.size());
  for (size_t i = 0; i < a.skin_bindings.size(); ++i) {
    EXPECT_EQ(a.skin_bindings[i].mesh, b.skin_bindings[i].mesh);
    ASSERT_EQ(a.skin_bindings[i].vertices.size(), b.skin_bindings[i].vertices.size());
    for (size_t k = 0; k < a.skin_bindings[i].vertices.size(); ++k) {
      const pd::AssetSkinVertex& x = a.skin_bindings[i].vertices[k];
      const pd::AssetSkinVertex& y = b.skin_bindings[i].vertices[k];
      EXPECT_EQ(x.muscle, y.muscle);
      ExpectSameTetPoint(x.point, y.point);
      EXPECT_EQ(x.falloff, y.falloff);
    }
  }
}

bool Deserialize(const std::vector<uint8_t>& bytes, pd::Asset* asset, std::string* error) {
  return pd::DeserializeAsset(bytes.data(), bytes.size(), asset, error);
}

void PutU32(std::vector<uint8_t>* bytes, size_t at, uint32_t v) {
  for (int i = 0; i < 4; ++i) (*bytes)[at + i] = static_cast<uint8_t>(v >> (8 * i));
}

}  // namespace

TEST(AssetTest, TestAssetIsValid) {
  std::string error;
  EXPECT_TRUE(pd::ValidateAsset(MakeTestAsset(), &error)) << error;
}

TEST(AssetTest, FileRoundTripIsExact) {
  const pd::Asset asset = MakeTestAsset();
  const std::string path = ::testing::TempDir() + "pd_asset_roundtrip.pdasset";
  std::string error;
  ASSERT_TRUE(pd::SaveAsset(path, asset, &error)) << error;
  pd::Asset loaded;
  ASSERT_TRUE(pd::LoadAsset(path, &loaded, &error)) << error;
  ExpectSameAsset(asset, loaded);
  std::remove(path.c_str());
}

TEST(AssetTest, SerializationIsStable) {
  // Serialize -> deserialize -> serialize reproduces the same bytes.
  const std::vector<uint8_t> bytes = pd::SerializeAsset(MakeTestAsset());
  pd::Asset loaded;
  std::string error;
  ASSERT_TRUE(Deserialize(bytes, &loaded, &error)) << error;
  EXPECT_EQ(pd::SerializeAsset(loaded), bytes);
}

TEST(AssetTest, EmptyAssetRoundTrips) {
  pd::Asset loaded = MakeTestAsset();
  std::string error;
  ASSERT_TRUE(Deserialize(pd::SerializeAsset(pd::Asset()), &loaded, &error)) << error;
  ExpectSameAsset(pd::Asset(), loaded);
}

TEST(AssetTest, SkipsUnknownChunks) {
  const pd::Asset asset = MakeTestAsset();
  std::vector<uint8_t> bytes = pd::SerializeAsset(asset);
  // Append a chunk with an unknown tag ("ZZZZ") and a 5-byte payload.
  const uint8_t chunk[] = {'Z', 'Z', 'Z', 'Z', 5, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5};
  bytes.insert(bytes.end(), std::begin(chunk), std::end(chunk));
  uint32_t count = 0;
  for (int i = 0; i < 4; ++i) count |= static_cast<uint32_t>(bytes[kChunkCountOffset + i]) << (8 * i);
  PutU32(&bytes, kChunkCountOffset, count + 1);

  pd::Asset loaded;
  std::string error;
  ASSERT_TRUE(Deserialize(bytes, &loaded, &error)) << error;
  ExpectSameAsset(asset, loaded);
}

TEST(AssetTest, RejectsBadMagicAndNewerVersion) {
  const std::vector<uint8_t> good = pd::SerializeAsset(MakeTestAsset());
  pd::Asset loaded;
  std::string error;

  std::vector<uint8_t> bytes = good;
  bytes[0] = 'X';
  EXPECT_FALSE(Deserialize(bytes, &loaded, &error));
  EXPECT_NE(error.find("magic"), std::string::npos) << error;

  bytes = good;
  PutU32(&bytes, 8, pd::kAssetVersion + 1);
  EXPECT_FALSE(Deserialize(bytes, &loaded, &error));
  EXPECT_NE(error.find("version"), std::string::npos) << error;

  bytes = good;
  bytes[8 + 4 + 4 + 8] = 2;  // Up axis +Z.
  EXPECT_FALSE(Deserialize(bytes, &loaded, &error));
}

TEST(AssetTest, RejectsEveryTruncation) {
  // Every proper prefix of a valid file must be rejected, never crash or read
  // out of bounds, and leave the output untouched.
  const std::vector<uint8_t> bytes = pd::SerializeAsset(MakeTestAsset());
  for (size_t n = 0; n < bytes.size(); ++n) {
    pd::Asset loaded;
    loaded.generator = "untouched";
    std::string error;
    const std::vector<uint8_t> prefix(bytes.begin(), bytes.begin() + n);
    ASSERT_FALSE(Deserialize(prefix, &loaded, &error)) << "prefix length " << n;
    EXPECT_FALSE(error.empty());
    EXPECT_EQ(loaded.generator, "untouched");
  }
}

TEST(AssetTest, RejectsHugeCounts) {
  // A corrupt array count must fail cleanly instead of allocating.
  std::vector<uint8_t> bytes = pd::SerializeAsset(pd::Asset());
  // First chunk is INFO: tag (4), size (8), then the generator string length.
  PutU32(&bytes, kChunkCountOffset + 4 + 4 + 8, 0xFFFFFFFFu);
  pd::Asset loaded;
  std::string error;
  EXPECT_FALSE(Deserialize(bytes, &loaded, &error));
}

TEST(AssetTest, ValidationCatchesBadReferences) {
  std::string error;
  auto expect_invalid = [&error](const pd::Asset& asset, const char* what) {
    EXPECT_FALSE(pd::ValidateAsset(asset, &error)) << what;
    EXPECT_FALSE(pd::SaveAsset(::testing::TempDir() + "pd_asset_invalid.pdasset", asset, &error)) << what;
  };

  pd::Asset asset = MakeTestAsset();
  asset.muscles[0].mesh.tets[0][2] = asset.muscles[0].mesh.NumVertices();
  expect_invalid(asset, "tet vertex out of range");
  EXPECT_NE(error.find("biceps"), std::string::npos) << error;

  asset = MakeTestAsset();
  std::swap(asset.muscles[0].mesh.tets[1][0], asset.muscles[0].mesh.tets[1][1]);
  expect_invalid(asset, "inverted tet");

  asset = MakeTestAsset();
  asset.muscles[0].attachments[0].bone = 2;
  expect_invalid(asset, "attachment bone out of range");

  asset = MakeTestAsset();
  asset.muscles[0].tet_material[0] = 2;
  expect_invalid(asset, "material index out of range");

  asset = MakeTestAsset();
  asset.muscles[0].mesh.fibres.pop_back();
  expect_invalid(asset, "fibre count");

  asset = MakeTestAsset();
  asset.bones[1].name = asset.bones[0].name;
  expect_invalid(asset, "duplicate bone name");

  asset = MakeTestAsset();
  asset.colliders[0].bone = 5;
  expect_invalid(asset, "collider bone out of range");

  asset = MakeTestAsset();
  asset.embeddings[0].vertices[0].tet = static_cast<uint32_t>(asset.muscles[0].mesh.NumTets());
  expect_invalid(asset, "embedding tet out of range");

  asset = MakeTestAsset();
  asset.skin_bindings[0].vertices[2].point.tet = 6;  // The second muscle has 6 tets.
  expect_invalid(asset, "skin tet out of range");

  asset = MakeTestAsset();
  asset.muscles[0].mesh.vertices(0, 1) = std::nan("");
  expect_invalid(asset, "NaN position");
}
