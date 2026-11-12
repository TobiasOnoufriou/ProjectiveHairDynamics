// Copyright 2022 Tobias Onoufriou
#include "PDAsset.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <iterator>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace pd {
namespace {

constexpr char kMagic[8] = {'P', 'D', 'A', 'S', 'S', 'E', 'T', '\0'};
constexpr uint8_t kUpAxisY = 1;
constexpr uint8_t kRightHanded = 0;
constexpr size_t kHeaderSize = 8 + 4 + 4 + 8 + 1 + 1 + 2 + 4;
constexpr size_t kChunkHeaderSize = 4 + 8;

constexpr uint32_t FourCC(const char (&s)[5]) {
  return static_cast<uint32_t>(static_cast<uint8_t>(s[0])) |
         static_cast<uint32_t>(static_cast<uint8_t>(s[1])) << 8 |
         static_cast<uint32_t>(static_cast<uint8_t>(s[2])) << 16 |
         static_cast<uint32_t>(static_cast<uint8_t>(s[3])) << 24;
}
constexpr uint32_t kTagInfo = FourCC("INFO");
constexpr uint32_t kTagBones = FourCC("BONE");
constexpr uint32_t kTagMuscle = FourCC("MUSC");  // One chunk per muscle, in order.
constexpr uint32_t kTagColliders = FourCC("COLL");
constexpr uint32_t kTagEmbedding = FourCC("EMBD");  // One chunk per embedding.
constexpr uint32_t kTagSkin = FourCC("SKIN");       // One chunk per skin binding.

// Collider shapes have explicit codes so reordering the enum cannot change
// the file format.
uint8_t ShapeCode(Collider::Shape shape) {
  switch (shape) {
    case Collider::Shape::kPlane: return 0;
    case Collider::Shape::kSphere: return 1;
    case Collider::Shape::kCapsule: return 2;
    case Collider::Shape::kBox: return 3;
  }
  return 0;
}
bool ShapeFromCode(uint8_t code, Collider::Shape* shape) {
  switch (code) {
    case 0: *shape = Collider::Shape::kPlane; return true;
    case 1: *shape = Collider::Shape::kSphere; return true;
    case 2: *shape = Collider::Shape::kCapsule; return true;
    case 3: *shape = Collider::Shape::kBox; return true;
    default: return false;
  }
}

// ---------------------------------------------------------------- Writing.

class Writer {
 public:
  std::vector<uint8_t> bytes;

  void U8(uint8_t v) { bytes.push_back(v); }
  void U16(uint16_t v) { UInt(v, 2); }
  void U32(uint32_t v) { UInt(v, 4); }
  void U64(uint64_t v) { UInt(v, 8); }
  void I32(int32_t v) { U32(static_cast<uint32_t>(v)); }
  void F64(double v) {
    uint64_t bits;
    std::memcpy(&bits, &v, sizeof(bits));
    U64(bits);
  }
  void Count(size_t n) { U32(static_cast<uint32_t>(n)); }
  void Str(const std::string& s) {
    Count(s.size());
    bytes.insert(bytes.end(), s.begin(), s.end());
  }
  void Vec3(const EigenVector3& v) {
    for (int i = 0; i < 3; ++i) F64(v[i]);
  }
  void Mat3(const EigenMatrix3& m) {  // Row-major.
    for (int r = 0; r < 3; ++r)
      for (int c = 0; c < 3; ++c) F64(m(r, c));
  }
  void TetPoint(const AssetTetPoint& p) {
    U32(p.tet);
    for (ScalarType w : p.weights) F64(w);
  }

  // Starts a chunk; EndChunk patches its payload size.
  size_t BeginChunk(uint32_t tag) {
    U32(tag);
    const size_t size_at = bytes.size();
    U64(0);
    return size_at;
  }
  void EndChunk(size_t size_at) {
    const uint64_t size = bytes.size() - size_at - 8;
    for (int i = 0; i < 8; ++i) bytes[size_at + i] = static_cast<uint8_t>(size >> (8 * i));
  }

 private:
  void UInt(uint64_t v, int n) {
    for (int i = 0; i < n; ++i) bytes.push_back(static_cast<uint8_t>(v >> (8 * i)));
  }
};

void WriteMuscle(const AssetMuscle& m, Writer* w) {
  w->Str(m.name);
  w->Str(m.activation_group);
  w->Count(m.mesh.vertices.rows());
  for (Eigen::Index i = 0; i < m.mesh.vertices.rows(); ++i) w->Vec3(m.mesh.vertices.row(i).transpose());
  w->Count(m.mesh.tets.size());
  for (const TetMesh::Tet& t : m.mesh.tets)
    for (unsigned int v : t) w->U32(v);
  w->Count(m.mesh.fibres.size());
  for (const EigenVector3& f : m.mesh.fibres) w->Vec3(f);
  w->Count(m.materials.size());
  for (const AssetMaterial& mat : m.materials) {
    for (ScalarType v : {mat.density, mat.strain_stiffness, mat.min_stretch, mat.max_stretch,
                         mat.volume_stiffness, mat.min_volume, mat.max_volume, mat.fibre_stiffness,
                         mat.max_contraction})
      w->F64(v);
  }
  w->Count(m.tet_material.size());
  for (uint32_t i : m.tet_material) w->U32(i);
  w->Count(m.fibre_coordinate.size());
  for (ScalarType u : m.fibre_coordinate) w->F64(u);
  w->Count(m.attachments.size());
  for (const AssetAttachment& a : m.attachments) {
    w->U32(a.vertex);
    w->U32(a.bone);
    w->Vec3(a.local_position);
  }
}

// ---------------------------------------------------------------- Reading.

// Bounds-checked cursor over [data, data + size). After the first failure every
// read returns zeros and ok() stays false, so parsers can read a whole record
// and check once.
class Reader {
 public:
  Reader(const uint8_t* data, size_t size) : m_data(data), m_size(size) {}

  bool ok() const { return m_ok; }
  const std::string& error() const { return m_error; }
  size_t remaining() const { return m_size - m_pos; }
  void Fail(const std::string& error) {
    if (m_ok) m_error = error;
    m_ok = false;
  }

  uint8_t U8() { return static_cast<uint8_t>(UInt(1)); }
  uint16_t U16() { return static_cast<uint16_t>(UInt(2)); }
  uint32_t U32() { return static_cast<uint32_t>(UInt(4)); }
  uint64_t U64() { return UInt(8); }
  int32_t I32() { return static_cast<int32_t>(U32()); }
  double F64() {
    const uint64_t bits = U64();
    double v;
    std::memcpy(&v, &bits, sizeof(v));
    return v;
  }
  // An array count, rejected if `count * min_item_bytes` cannot fit in the
  // remaining bytes (so a corrupt count never triggers a huge allocation).
  size_t Count(size_t min_item_bytes) {
    const uint32_t n = U32();
    if (!m_ok) return 0;
    if (min_item_bytes > 0 && n > remaining() / min_item_bytes) {
      Fail("array count exceeds the data left (truncated or corrupt file)");
      return 0;
    }
    return n;
  }
  std::string Str() {
    const size_t n = Count(1);
    if (!m_ok) return std::string();
    std::string s(reinterpret_cast<const char*>(m_data + m_pos), n);
    m_pos += n;
    return s;
  }
  bool Bytes(void* out, size_t n) {
    if (!Need(n)) return false;
    std::memcpy(out, m_data + m_pos, n);
    m_pos += n;
    return true;
  }
  EigenVector3 Vec3() {
    EigenVector3 v;
    for (int i = 0; i < 3; ++i) v[i] = F64();
    return v;
  }
  EigenMatrix3 Mat3() {
    EigenMatrix3 m;
    for (int r = 0; r < 3; ++r)
      for (int c = 0; c < 3; ++c) m(r, c) = F64();
    return m;
  }
  AssetTetPoint TetPoint() {
    AssetTetPoint p;
    p.tet = U32();
    for (ScalarType& w : p.weights) w = F64();
    return p;
  }
  // A sub-reader over the next n bytes, which this reader skips.
  Reader Sub(size_t n) {
    if (!Need(n)) return Reader(nullptr, 0);
    Reader sub(m_data + m_pos, n);
    m_pos += n;
    return sub;
  }

 private:
  bool Need(size_t n) {
    if (m_ok && n > remaining()) Fail("unexpected end of data (truncated file)");
    return m_ok;
  }
  uint64_t UInt(int n) {
    if (!Need(static_cast<size_t>(n))) return 0;
    uint64_t v = 0;
    for (int i = 0; i < n; ++i) v |= static_cast<uint64_t>(m_data[m_pos + i]) << (8 * i);
    m_pos += n;
    return v;
  }

  const uint8_t* m_data;
  size_t m_size;
  size_t m_pos = 0;
  bool m_ok = true;
  std::string m_error;
};

void ReadMuscle(Reader* r, AssetMuscle* m) {
  m->name = r->Str();
  m->activation_group = r->Str();
  const size_t nv = r->Count(3 * 8);
  m->mesh.vertices.resize(static_cast<Eigen::Index>(nv), 3);
  for (size_t i = 0; i < nv; ++i) m->mesh.vertices.row(static_cast<Eigen::Index>(i)) = r->Vec3().transpose();
  m->mesh.tets.resize(r->Count(4 * 4));
  for (TetMesh::Tet& t : m->mesh.tets)
    for (unsigned int& v : t) v = r->U32();
  m->mesh.fibres.resize(r->Count(3 * 8));
  for (EigenVector3& f : m->mesh.fibres) f = r->Vec3();
  m->materials.resize(r->Count(9 * 8));
  for (AssetMaterial& mat : m->materials) {
    for (ScalarType* v : {&mat.density, &mat.strain_stiffness, &mat.min_stretch, &mat.max_stretch,
                          &mat.volume_stiffness, &mat.min_volume, &mat.max_volume, &mat.fibre_stiffness,
                          &mat.max_contraction})
      *v = r->F64();
  }
  m->tet_material.resize(r->Count(4));
  for (uint32_t& i : m->tet_material) i = r->U32();
  m->fibre_coordinate.resize(r->Count(8));
  for (ScalarType& u : m->fibre_coordinate) u = r->F64();
  m->attachments.resize(r->Count(4 + 4 + 3 * 8));
  for (AssetAttachment& a : m->attachments) {
    a.vertex = r->U32();
    a.bone = r->U32();
    a.local_position = r->Vec3();
  }
}

// Parses one chunk payload into *asset. Unknown tags are ignored.
void ReadChunk(uint32_t tag, Reader* r, Asset* asset, bool* seen_info, bool* seen_bones,
               bool* seen_colliders) {
  auto once = [r](bool* seen, const char* name) {
    if (*seen) r->Fail(std::string("duplicate ") + name + " chunk");
    *seen = true;
  };
  if (tag == kTagInfo) {
    once(seen_info, "INFO");
    asset->generator = r->Str();
    asset->source = r->Str();
  } else if (tag == kTagBones) {
    once(seen_bones, "BONE");
    asset->bones.resize(r->Count(4 + 12 * 8));
    for (AssetBone& b : asset->bones) {
      b.name = r->Str();
      b.rest_rotation = r->Mat3();
      b.rest_translation = r->Vec3();
    }
  } else if (tag == kTagMuscle) {
    asset->muscles.emplace_back();
    ReadMuscle(r, &asset->muscles.back());
  } else if (tag == kTagColliders) {
    once(seen_colliders, "COLL");
    asset->colliders.resize(r->Count(4 + 1 + 18 * 8));
    for (AssetCollider& c : asset->colliders) {
      c.bone = r->I32();
      if (!ShapeFromCode(r->U8(), &c.shape.shape)) r->Fail("unknown collider shape");
      c.shape.rotation = r->Mat3();
      c.shape.translation = r->Vec3();
      c.shape.radius = r->F64();
      c.shape.half_length = r->F64();
      c.shape.half_extents = r->Vec3();
      c.shape.friction = r->F64();
    }
  } else if (tag == kTagEmbedding) {
    asset->embeddings.emplace_back();
    AssetEmbedding& e = asset->embeddings.back();
    e.mesh = r->Str();
    e.muscle = r->U32();
    e.vertices.resize(r->Count(4 + 4 * 8));
    for (AssetTetPoint& p : e.vertices) p = r->TetPoint();
  } else if (tag == kTagSkin) {
    asset->skin_bindings.emplace_back();
    AssetSkinBinding& s = asset->skin_bindings.back();
    s.mesh = r->Str();
    s.vertices.resize(r->Count(4 + 4 + 4 * 8 + 8));
    for (AssetSkinVertex& v : s.vertices) {
      v.muscle = r->I32();
      v.point = r->TetPoint();
      v.falloff = r->F64();
    }
  } else {
    return;  // Unknown chunk from a newer writer: skipped.
  }
  if (r->ok() && r->remaining() != 0) r->Fail("chunk has trailing bytes");
}

// ---------------------------------------------------------------- Validation.

bool Finite(ScalarType v) { return std::isfinite(v); }
template <typename Derived>
bool Finite(const Eigen::MatrixBase<Derived>& m) {
  return m.allFinite();
}

bool Fail(std::string* error, const std::string& message) {
  if (error) *error = message;
  return false;
}

bool ValidTetPoint(const AssetTetPoint& p, const AssetMuscle& muscle) {
  if (p.tet >= muscle.mesh.tets.size()) return false;
  for (ScalarType w : p.weights)
    if (!Finite(w)) return false;
  return true;
}

bool ValidateMuscle(const AssetMuscle& m, size_t num_bones, std::string* error) {
  const std::string where = "muscle '" + m.name + "': ";
  const size_t nv = static_cast<size_t>(m.mesh.vertices.rows());
  if (!Finite(m.mesh.vertices)) return Fail(error, where + "non-finite vertex position");
  for (size_t t = 0; t < m.mesh.tets.size(); ++t) {
    const TetMesh::Tet& tet = m.mesh.tets[t];
    for (unsigned int v : tet)
      if (v >= nv) return Fail(error, where + "tet " + std::to_string(t) + " has an out-of-range vertex");
    if (!(m.mesh.SignedVolume(t) > 0.0))
      return Fail(error, where + "tet " + std::to_string(t) + " is inverted or degenerate");
  }
  if (!m.mesh.fibres.empty()) {
    if (m.mesh.fibres.size() != m.mesh.tets.size())
      return Fail(error, where + "fibre count does not match the tet count");
    for (const EigenVector3& f : m.mesh.fibres)
      if (!Finite(f)) return Fail(error, where + "non-finite fibre direction");
  }
  if (!m.mesh.tets.empty() && m.materials.empty()) return Fail(error, where + "has tets but no material");
  for (const AssetMaterial& mat : m.materials) {
    const ScalarType values[] = {mat.density, mat.strain_stiffness, mat.min_stretch, mat.max_stretch,
                                 mat.volume_stiffness, mat.min_volume, mat.max_volume, mat.fibre_stiffness,
                                 mat.max_contraction};
    for (ScalarType v : values)
      if (!Finite(v)) return Fail(error, where + "non-finite material value");
    if (!(mat.density > 0.0) || mat.strain_stiffness < 0.0 || mat.volume_stiffness < 0.0 ||
        mat.fibre_stiffness < 0.0)
      return Fail(error, where + "material needs a positive density and non-negative stiffnesses");
    if (mat.min_stretch > mat.max_stretch || mat.min_volume > mat.max_volume)
      return Fail(error, where + "material has an empty stretch or volume range");
  }
  if (!m.tet_material.empty()) {
    if (m.tet_material.size() != m.mesh.tets.size())
      return Fail(error, where + "per-tet material count does not match the tet count");
    for (uint32_t i : m.tet_material)
      if (i >= m.materials.size()) return Fail(error, where + "per-tet material index out of range");
  }
  if (!m.fibre_coordinate.empty()) {
    if (m.fibre_coordinate.size() != nv)
      return Fail(error, where + "fibre coordinate count does not match the vertex count");
    for (ScalarType u : m.fibre_coordinate)
      if (!Finite(u)) return Fail(error, where + "non-finite fibre coordinate");
  }
  for (const AssetAttachment& a : m.attachments) {
    if (a.vertex >= nv) return Fail(error, where + "attachment vertex out of range");
    if (a.bone >= num_bones) return Fail(error, where + "attachment bone out of range");
    if (!Finite(a.local_position)) return Fail(error, where + "non-finite attachment position");
  }
  return true;
}

}  // namespace

bool ValidateAsset(const Asset& asset, std::string* error) {
  std::set<std::string> names;
  for (const AssetBone& b : asset.bones) {
    if (b.name.empty()) return Fail(error, "bone with an empty name");
    if (!names.insert(b.name).second) return Fail(error, "duplicate bone name '" + b.name + "'");
    if (!Finite(b.rest_rotation) || !Finite(b.rest_translation))
      return Fail(error, "bone '" + b.name + "': non-finite rest transform");
  }
  for (const AssetMuscle& m : asset.muscles)
    if (!ValidateMuscle(m, asset.bones.size(), error)) return false;
  for (size_t i = 0; i < asset.colliders.size(); ++i) {
    const AssetCollider& c = asset.colliders[i];
    const Collider& s = c.shape;
    const std::string where = "collider " + std::to_string(i) + ": ";
    if (c.bone < -1 || c.bone >= static_cast<int64_t>(asset.bones.size()))
      return Fail(error, where + "bone out of range");
    if (!Finite(s.rotation) || !Finite(s.translation) || !Finite(s.radius) || !Finite(s.half_length) ||
        !Finite(s.half_extents) || !Finite(s.friction))
      return Fail(error, where + "non-finite value");
    if (s.radius < 0.0 || s.half_length < 0.0 || (s.half_extents.array() < 0.0).any() || s.friction < 0.0)
      return Fail(error, where + "negative size or friction");
  }
  for (const AssetEmbedding& e : asset.embeddings) {
    const std::string where = "embedding '" + e.mesh + "': ";
    if (e.muscle >= asset.muscles.size()) return Fail(error, where + "muscle out of range");
    for (const AssetTetPoint& p : e.vertices)
      if (!ValidTetPoint(p, asset.muscles[e.muscle])) return Fail(error, where + "invalid tet point");
  }
  for (const AssetSkinBinding& s : asset.skin_bindings) {
    const std::string where = "skin binding '" + s.mesh + "': ";
    for (const AssetSkinVertex& v : s.vertices) {
      if (!Finite(v.falloff)) return Fail(error, where + "non-finite falloff");
      if (v.muscle == -1) continue;
      if (v.muscle < 0 || v.muscle >= static_cast<int64_t>(asset.muscles.size()))
        return Fail(error, where + "muscle out of range");
      if (!ValidTetPoint(v.point, asset.muscles[v.muscle])) return Fail(error, where + "invalid tet point");
    }
  }
  return true;
}

std::vector<uint8_t> SerializeAsset(const Asset& asset) {
  Writer w;
  const uint32_t chunk_count = 3 + static_cast<uint32_t>(asset.muscles.size() + asset.embeddings.size() +
                                                         asset.skin_bindings.size());
  w.bytes.insert(w.bytes.end(), std::begin(kMagic), std::end(kMagic));
  w.U32(kAssetVersion);
  w.U32(0);    // Flags.
  w.F64(1.0);  // Meters per unit.
  w.U8(kUpAxisY);
  w.U8(kRightHanded);
  w.U16(0);
  w.U32(chunk_count);

  size_t chunk = w.BeginChunk(kTagInfo);
  w.Str(asset.generator);
  w.Str(asset.source);
  w.EndChunk(chunk);

  chunk = w.BeginChunk(kTagBones);
  w.Count(asset.bones.size());
  for (const AssetBone& b : asset.bones) {
    w.Str(b.name);
    w.Mat3(b.rest_rotation);
    w.Vec3(b.rest_translation);
  }
  w.EndChunk(chunk);

  for (const AssetMuscle& m : asset.muscles) {
    chunk = w.BeginChunk(kTagMuscle);
    WriteMuscle(m, &w);
    w.EndChunk(chunk);
  }

  chunk = w.BeginChunk(kTagColliders);
  w.Count(asset.colliders.size());
  for (const AssetCollider& c : asset.colliders) {
    w.I32(c.bone);
    w.U8(ShapeCode(c.shape.shape));
    w.Mat3(c.shape.rotation);
    w.Vec3(c.shape.translation);
    w.F64(c.shape.radius);
    w.F64(c.shape.half_length);
    w.Vec3(c.shape.half_extents);
    w.F64(c.shape.friction);
  }
  w.EndChunk(chunk);

  for (const AssetEmbedding& e : asset.embeddings) {
    chunk = w.BeginChunk(kTagEmbedding);
    w.Str(e.mesh);
    w.U32(e.muscle);
    w.Count(e.vertices.size());
    for (const AssetTetPoint& p : e.vertices) w.TetPoint(p);
    w.EndChunk(chunk);
  }

  for (const AssetSkinBinding& s : asset.skin_bindings) {
    chunk = w.BeginChunk(kTagSkin);
    w.Str(s.mesh);
    w.Count(s.vertices.size());
    for (const AssetSkinVertex& v : s.vertices) {
      w.I32(v.muscle);
      w.TetPoint(v.point);
      w.F64(v.falloff);
    }
    w.EndChunk(chunk);
  }
  return std::move(w.bytes);
}

bool DeserializeAsset(const uint8_t* data, size_t size, Asset* asset, std::string* error) {
  if (size < kHeaderSize || std::memcmp(data, kMagic, sizeof(kMagic)) != 0)
    return Fail(error, "not a .pdasset file (bad magic)");
  Reader r(data + sizeof(kMagic), size - sizeof(kMagic));
  const uint32_t version = r.U32();
  if (version == 0 || version > kAssetVersion)
    return Fail(error, "unsupported .pdasset version " + std::to_string(version) + " (this build reads up to " +
                           std::to_string(kAssetVersion) + "); re-bake the asset or update the runtime");
  if (r.U32() != 0) return Fail(error, "unknown .pdasset flags");
  const double meters_per_unit = r.F64();
  const uint8_t up_axis = r.U8();
  const uint8_t handedness = r.U8();
  r.U16();
  if (meters_per_unit != 1.0 || up_axis != kUpAxisY || handedness != kRightHanded)
    return Fail(error, ".pdasset is not in meters, +Y up, right-handed");
  const uint32_t chunk_count = r.U32();

  Asset result;
  bool seen_info = false, seen_bones = false, seen_colliders = false;
  for (uint32_t c = 0; c < chunk_count && r.ok(); ++c) {
    const uint32_t tag = r.U32();
    const uint64_t payload = r.U64();
    if (!r.ok()) break;
    if (payload > r.remaining()) {
      r.Fail("chunk extends past the end of the file (truncated file)");
      break;
    }
    Reader chunk = r.Sub(static_cast<size_t>(payload));
    ReadChunk(tag, &chunk, &result, &seen_info, &seen_bones, &seen_colliders);
    if (!chunk.ok()) return Fail(error, chunk.error());
  }
  if (!r.ok()) return Fail(error, r.error());
  if (r.remaining() != 0) return Fail(error, "trailing bytes after the last chunk");
  if (!ValidateAsset(result, error)) return false;
  *asset = std::move(result);
  return true;
}

bool SaveAsset(const std::string& path, const Asset& asset, std::string* error) {
  if (!ValidateAsset(asset, error)) return false;
  const std::vector<uint8_t> bytes = SerializeAsset(asset);
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) return Fail(error, "cannot open '" + path + "' for writing");
  out.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
  if (!out) return Fail(error, "failed writing '" + path + "'");
  return true;
}

bool LoadAsset(const std::string& path, Asset* asset, std::string* error) {
  std::ifstream in(path, std::ios::binary);
  if (!in) return Fail(error, "cannot open '" + path + "'");
  const std::vector<uint8_t> bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  if (in.bad()) return Fail(error, "failed reading '" + path + "'");
  if (!DeserializeAsset(bytes.data(), bytes.size(), asset, error)) {
    if (error) *error = path + ": " + *error;
    return false;
  }
  return true;
}

}  // namespace pd
