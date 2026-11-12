#pragma once
#include <array>
#include <functional>
#include <string>
#include <vector>

#include "math_headers.h"

// A tetrahedral volume mesh: vertex positions plus 4 vertex indices per tet.
// No OpenGL here, so it can be built and tested headless.
//
// Convention: every tet (a, b, c, d) is positively oriented, i.e.
// det[b - a, c - a, d - a] > 0. FixOrientation() enforces this after loading,
// and the boundary faces rely on it to point outwards.
struct TetMesh {
  using Tet = std::array<unsigned int, 4>;
  using Triangle = std::array<unsigned int, 3>;

  EigenMatrix3x vertices;  // N x 3, row i = vertex i.
  std::vector<Tet> tets;
  // Optional muscle fibre direction per tet (unit vector in rest space), same
  // order as tets. Empty if the mesh has no fibres (e.g. loaded from TetGen).
  std::vector<EigenVector3> fibres;

  unsigned int NumVertices() const { return static_cast<unsigned int>(vertices.rows()); }
  size_t NumTets() const { return tets.size(); }

  // Signed volume det[b - a, c - a, d - a] / 6 of tet t at positions x.
  ScalarType SignedVolume(size_t t, const EigenMatrix3x& x) const;
  ScalarType SignedVolume(size_t t) const { return SignedVolume(t, vertices); }
  ScalarType TotalVolume() const;
  // Swaps two vertices of every negatively oriented tet. Returns how many were flipped.
  size_t FixOrientation();

  // Mean of the 4 vertices of tet t at positions x.
  EigenVector3 Centroid(size_t t, const EigenMatrix3x& x) const;
  EigenVector3 Centroid(size_t t) const { return Centroid(t, vertices); }
  // Samples a direction field at every tet centroid (rest pose) and stores it,
  // normalised, as the tet's fibre direction. A zero direction stays zero.
  void SetFibresFromField(const std::function<EigenVector3(const EigenVector3&)>& field);

  // Lumped masses: each tet's mass density * volume is split equally over its 4 vertices.
  VectorX LumpedMasses(ScalarType density) const;

  // Boundary surface: faces that belong to exactly one tet, wound counter-clockwise
  // when seen from outside (normal = (b - a) x (c - a) points out of the volume).
  std::vector<Triangle> BoundaryTriangles() const;
  // Unique tet edges as index pairs (for a wireframe and for edge springs).
  std::vector<std::array<unsigned int, 2>> Edges() const;
};

// Box [origin, origin + size] split into nx * ny * nz cubes, each cut into 6 tets
// sharing the cube's main diagonal (Kuhn / Freudenthal subdivision). Every cube is
// cut the same way, so neighbouring cubes share their face diagonals and the mesh
// is conforming. Vertex (i, j, k) has index i + (nx + 1) * (j + (ny + 1) * k).
// Fibres run along the longest axis of the box (x on ties).
TetMesh MakeTetBeam(unsigned int nx, unsigned int ny, unsigned int nz, const EigenVector3& size,
                    const EigenVector3& origin = EigenVector3::Zero());

// Reads a TetGen mesh (https://wias-berlin.de/software/tetgen/fformats.html):
// `<path>.node` holds "<#points> 3 <#attributes> <has marker>" then
// "<id> x y z [attributes] [marker]", and `<path>.ele` holds "<#tets> <4|10> <#attributes>"
// then "<id> n0 n1 n2 n3 ...". '#' starts a comment. Ids may be 0- or 1-based (TetGen
// uses the first point id as the base); quadratic (10-node) tets keep their corners.
// Tets are re-oriented to positive volume. Returns false and sets *error on failure.
bool LoadTetGen(const std::string& path_without_extension, TetMesh* mesh, std::string* error);
// Writes `<path>.node` / `<path>.ele` (0-based), readable by LoadTetGen and TetGen tools.
bool SaveTetGen(const std::string& path_without_extension, const TetMesh& mesh);

// Area-weighted vertex normals of a triangle surface at positions x (N x 3).
// Vertices not on the surface get a zero normal.
EigenMatrix3x VertexNormals(const EigenMatrix3x& x, const std::vector<TetMesh::Triangle>& triangles);
