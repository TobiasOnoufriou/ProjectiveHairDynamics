#include "TetMesh.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

// Faces of a positively oriented tet (a, b, c, d), each wound so that its normal
// points away from the opposite vertex.
constexpr int kTetFaces[4][3] = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}};

// Order-independent key of a face or edge (vertex ids sorted, packed in 64 bits;
// 21 bits per id, so up to 2M vertices).
uint64_t FaceKey(unsigned int a, unsigned int b, unsigned int c) {
  unsigned int v[3] = {a, b, c};
  std::sort(v, v + 3);
  return (static_cast<uint64_t>(v[0]) << 42) | (static_cast<uint64_t>(v[1]) << 21) | v[2];
}

// Reads the next non-empty line with '#' comments removed into *line.
bool NextDataLine(std::istream& in, std::istringstream* line) {
  std::string text;
  while (std::getline(in, text)) {
    const size_t hash = text.find('#');
    if (hash != std::string::npos) text.erase(hash);
    if (text.find_first_not_of(" \t\r") == std::string::npos) continue;
    line->clear();
    line->str(text);
    return true;
  }
  return false;
}

}  // namespace

ScalarType TetMesh::SignedVolume(size_t t, const EigenMatrix3x& x) const {
  const Tet& tet = tets[t];
  const EigenVector3 a = x.row(tet[0]).transpose();
  EigenMatrix3 d;
  d.col(0) = x.row(tet[1]).transpose() - a;
  d.col(1) = x.row(tet[2]).transpose() - a;
  d.col(2) = x.row(tet[3]).transpose() - a;
  return d.determinant() / 6.0;
}

ScalarType TetMesh::TotalVolume() const {
  ScalarType v = 0.0;
  for (size_t t = 0; t < tets.size(); ++t) v += SignedVolume(t);
  return v;
}

size_t TetMesh::FixOrientation() {
  size_t flipped = 0;
  for (size_t t = 0; t < tets.size(); ++t) {
    if (SignedVolume(t) < 0.0) {
      std::swap(tets[t][2], tets[t][3]);
      ++flipped;
    }
  }
  return flipped;
}

VectorX TetMesh::LumpedMasses(ScalarType density) const {
  VectorX m = VectorX::Zero(vertices.rows());
  for (size_t t = 0; t < tets.size(); ++t) {
    const ScalarType quarter = 0.25 * density * std::abs(SignedVolume(t));
    for (unsigned int v : tets[t]) m[v] += quarter;
  }
  return m;
}

std::vector<TetMesh::Triangle> TetMesh::BoundaryTriangles() const {
  // Interior faces are shared by exactly two tets; boundary faces by one.
  std::unordered_map<uint64_t, int> count;
  count.reserve(tets.size() * 4);
  for (const Tet& tet : tets) {
    for (const auto& f : kTetFaces) ++count[FaceKey(tet[f[0]], tet[f[1]], tet[f[2]])];
  }
  std::vector<Triangle> boundary;
  for (const Tet& tet : tets) {
    for (const auto& f : kTetFaces) {
      if (count[FaceKey(tet[f[0]], tet[f[1]], tet[f[2]])] == 1) {
        boundary.push_back({tet[f[0]], tet[f[1]], tet[f[2]]});
      }
    }
  }
  return boundary;
}

std::vector<std::array<unsigned int, 2>> TetMesh::Edges() const {
  std::vector<std::array<unsigned int, 2>> edges;
  edges.reserve(tets.size() * 6);
  for (const Tet& tet : tets) {
    for (int i = 0; i < 4; ++i) {
      for (int j = i + 1; j < 4; ++j) {
        edges.push_back({std::min(tet[i], tet[j]), std::max(tet[i], tet[j])});
      }
    }
  }
  std::sort(edges.begin(), edges.end());
  edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
  return edges;
}

TetMesh MakeTetBeam(unsigned int nx, unsigned int ny, unsigned int nz, const EigenVector3& size,
                    const EigenVector3& origin) {
  TetMesh mesh;
  nx = std::max(nx, 1u);
  ny = std::max(ny, 1u);
  nz = std::max(nz, 1u);
  const unsigned int vx = nx + 1, vy = ny + 1, vz = nz + 1;
  auto index = [&](unsigned int i, unsigned int j, unsigned int k) { return i + vx * (j + vy * k); };

  mesh.vertices.resize(vx * vy * vz, 3);
  for (unsigned int k = 0; k < vz; ++k)
    for (unsigned int j = 0; j < vy; ++j)
      for (unsigned int i = 0; i < vx; ++i) {
        mesh.vertices.row(index(i, j, k)) << origin.x() + size.x() * i / nx,
            origin.y() + size.y() * j / ny, origin.z() + size.z() * k / nz;
      }

  // Kuhn subdivision: one tet per permutation (a, b, c) of the axes, walking from
  // corner 000 to corner 111 by adding axis a, then b, then c. Corner bits: x = 1, y = 2, z = 4.
  constexpr int kPerms[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
  mesh.tets.reserve(6 * nx * ny * nz);
  for (unsigned int k = 0; k < nz; ++k)
    for (unsigned int j = 0; j < ny; ++j)
      for (unsigned int i = 0; i < nx; ++i) {
        auto corner = [&](int bits) { return index(i + (bits & 1), j + ((bits >> 1) & 1), k + ((bits >> 2) & 1)); };
        for (const auto& p : kPerms) {
          const int c1 = 1 << p[0];
          const int c2 = c1 | (1 << p[1]);
          mesh.tets.push_back({corner(0), corner(c1), corner(c2), corner(7)});
        }
      }
  mesh.FixOrientation();  // Odd permutations come out negatively oriented.
  return mesh;
}

bool LoadTetGen(const std::string& path, TetMesh* mesh, std::string* error) {
  auto fail = [&](const std::string& message) {
    if (error) *error = message;
    return false;
  };

  std::ifstream node_file(path + ".node");
  if (!node_file) return fail("cannot open " + path + ".node");
  std::istringstream line;
  int num_points = 0, dim = 0, num_attributes = 0, has_marker = 0;
  if (!NextDataLine(node_file, &line) || !(line >> num_points >> dim) || num_points <= 0 || dim != 3) {
    return fail(path + ".node: bad header (expected '<#points> 3 ...')");
  }
  line >> num_attributes >> has_marker;  // Optional.

  TetMesh result;
  result.vertices.resize(num_points, 3);
  int base = 0;
  for (int i = 0; i < num_points; ++i) {
    int id = 0;
    double x = 0, y = 0, z = 0;
    if (!NextDataLine(node_file, &line) || !(line >> id >> x >> y >> z)) {
      return fail(path + ".node: point " + std::to_string(i) + " is missing or malformed");
    }
    if (i == 0) base = id;  // TetGen numbers everything from the first point's id.
    if (id - base != i) return fail(path + ".node: point ids must be consecutive");
    result.vertices.row(i) << x, y, z;
  }

  std::ifstream ele_file(path + ".ele");
  if (!ele_file) return fail("cannot open " + path + ".ele");
  int num_tets = 0, nodes_per_tet = 0;
  if (!NextDataLine(ele_file, &line) || !(line >> num_tets >> nodes_per_tet) || num_tets <= 0 ||
      (nodes_per_tet != 4 && nodes_per_tet != 10)) {
    return fail(path + ".ele: bad header (expected '<#tets> 4|10 ...')");
  }
  result.tets.resize(num_tets);
  for (int t = 0; t < num_tets; ++t) {
    int id = 0;
    int n[4];
    if (!NextDataLine(ele_file, &line) || !(line >> id >> n[0] >> n[1] >> n[2] >> n[3])) {
      return fail(path + ".ele: tet " + std::to_string(t) + " is missing or malformed");
    }
    for (int c = 0; c < 4; ++c) {
      const int v = n[c] - base;
      if (v < 0 || v >= num_points) return fail(path + ".ele: tet " + std::to_string(t) + " has a bad vertex id");
      result.tets[t][c] = static_cast<unsigned int>(v);
    }
  }

  result.FixOrientation();
  *mesh = std::move(result);
  return true;
}

bool SaveTetGen(const std::string& path, const TetMesh& mesh) {
  std::ofstream node_file(path + ".node");
  std::ofstream ele_file(path + ".ele");
  if (!node_file || !ele_file) return false;
  node_file << mesh.NumVertices() << " 3 0 0\n" << std::setprecision(17);
  for (unsigned int i = 0; i < mesh.NumVertices(); ++i) {
    node_file << i << ' ' << mesh.vertices(i, 0) << ' ' << mesh.vertices(i, 1) << ' ' << mesh.vertices(i, 2)
              << '\n';
  }
  ele_file << mesh.NumTets() << " 4 0\n";
  for (size_t t = 0; t < mesh.NumTets(); ++t) {
    const TetMesh::Tet& tet = mesh.tets[t];
    ele_file << t << ' ' << tet[0] << ' ' << tet[1] << ' ' << tet[2] << ' ' << tet[3] << '\n';
  }
  return static_cast<bool>(node_file) && static_cast<bool>(ele_file);
}

EigenMatrix3x VertexNormals(const EigenMatrix3x& x, const std::vector<TetMesh::Triangle>& triangles) {
  EigenMatrix3x n = EigenMatrix3x::Zero(x.rows(), 3);
  for (const TetMesh::Triangle& tri : triangles) {
    const EigenVector3 a = x.row(tri[0]).transpose();
    // Unnormalised face normal: its length is twice the area, so this is area-weighted.
    const EigenVector3 face = (x.row(tri[1]).transpose() - a).cross(x.row(tri[2]).transpose() - a);
    for (unsigned int v : tri) n.row(v) += face.transpose();
  }
  for (Eigen::Index i = 0; i < n.rows(); ++i) {
    const ScalarType len = n.row(i).norm();
    if (len > 0.0) n.row(i) /= len;
  }
  return n;
}
