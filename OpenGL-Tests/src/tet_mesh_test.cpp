// Copyright 2022 Tobias Onoufriou
// Tests for tetrahedral meshes: the procedural beam, boundary extraction,
// TetGen I/O and the tet beam scene running on the PD solver.
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "gtest/gtest.h"
#include "PDSolver.h"
#include "TetMesh.h"
#include "TetScene.h"

namespace {

// Closed-surface check: every undirected edge of a closed, consistently wound
// surface is used exactly twice, once in each direction.
bool IsClosedAndConsistent(const std::vector<TetMesh::Triangle>& tris) {
  std::map<std::pair<unsigned int, unsigned int>, int> directed;
  for (const auto& t : tris) {
    for (int i = 0; i < 3; ++i) ++directed[{t[i], t[(i + 1) % 3]}];
  }
  for (const auto& [edge, count] : directed) {
    if (count != 1) return false;
    auto it = directed.find({edge.second, edge.first});
    if (it == directed.end() || it->second != 1) return false;
  }
  return true;
}

// Volume enclosed by a closed surface (divergence theorem); positive when the
// triangles are wound counter-clockwise seen from outside.
ScalarType EnclosedVolume(const EigenMatrix3x& x, const std::vector<TetMesh::Triangle>& tris) {
  ScalarType v = 0.0;
  for (const auto& t : tris) {
    v += x.row(t[0]).dot(x.row(t[1]).cross(x.row(t[2]))) / 6.0;
  }
  return v;
}

std::string TempPath(const char* name) {
  return ::testing::TempDir() + name;
}

}  // namespace

TEST(TetMeshTest, BeamCountsAndVolume) {
  const TetMesh mesh = MakeTetBeam(4, 2, 3, EigenVector3(2.0, 1.0, 1.5));
  EXPECT_EQ(mesh.NumVertices(), 5u * 3u * 4u);
  EXPECT_EQ(mesh.NumTets(), 6u * 4u * 2u * 3u);
  EXPECT_NEAR(mesh.TotalVolume(), 2.0 * 1.0 * 1.5, 1e-12);
  for (size_t t = 0; t < mesh.NumTets(); ++t) EXPECT_GT(mesh.SignedVolume(t), 0.0) << "tet " << t;
}

// The 6 Kuhn tets of a unit cube all have volume 1/6.
TEST(TetMeshTest, KuhnTetsAreEqualVolume) {
  const TetMesh mesh = MakeTetBeam(1, 1, 1, EigenVector3(1.0, 1.0, 1.0));
  ASSERT_EQ(mesh.NumTets(), 6u);
  for (size_t t = 0; t < 6; ++t) EXPECT_NEAR(mesh.SignedVolume(t), 1.0 / 6.0, 1e-15);
}

TEST(TetMeshTest, LumpedMassesSumToTotalMass) {
  const TetMesh mesh = MakeTetBeam(5, 2, 2, EigenVector3(1.0, 0.2, 0.2));
  const VectorX m = mesh.LumpedMasses(1000.0);
  EXPECT_NEAR(m.sum(), 1000.0 * 0.04, 1e-9);
  EXPECT_GT(m.minCoeff(), 0.0);
}

// A box surface split into triangles: 2 per cell face on each of the 6 sides,
// closed, outward-wound, and enclosing the box volume.
TEST(TetMeshTest, BoundaryIsClosedOutwardSurface) {
  const unsigned int nx = 4, ny = 2, nz = 3;
  const TetMesh mesh = MakeTetBeam(nx, ny, nz, EigenVector3(2.0, 1.0, 1.5), EigenVector3(0.3, -0.2, 0.1));
  const auto tris = mesh.BoundaryTriangles();
  EXPECT_EQ(tris.size(), 2u * 2u * (nx * ny + ny * nz + nx * nz));
  EXPECT_TRUE(IsClosedAndConsistent(tris));
  EXPECT_NEAR(EnclosedVolume(mesh.vertices, tris), mesh.TotalVolume(), 1e-12);

  // Interior vertices (the beam has (nx-1)(ny-1)(nz-1) of them) are not on the surface.
  std::set<unsigned int> on_surface;
  for (const auto& t : tris) on_surface.insert(t.begin(), t.end());
  EXPECT_EQ(on_surface.size(), mesh.NumVertices() - (nx - 1) * (ny - 1) * (nz - 1));
}

// Euler's formula for a tet mesh of a ball: V - E + F - T = 1.
TEST(TetMeshTest, EdgesSatisfyEulerCharacteristic) {
  const TetMesh mesh = MakeTetBeam(3, 2, 2, EigenVector3(1.0, 1.0, 1.0));
  const size_t v = mesh.NumVertices(), e = mesh.Edges().size(), t = mesh.NumTets();
  const size_t boundary = mesh.BoundaryTriangles().size();
  const size_t faces = (4 * t + boundary) / 2;  // Interior faces are shared by 2 tets.
  EXPECT_EQ(static_cast<long long>(v) - e + faces - t, 1);
}

TEST(TetMeshTest, FixOrientationFlipsInvertedTets) {
  TetMesh mesh = MakeTetBeam(2, 1, 1, EigenVector3(1.0, 1.0, 1.0));
  std::swap(mesh.tets[0][0], mesh.tets[0][1]);
  std::swap(mesh.tets[3][2], mesh.tets[3][3]);
  EXPECT_LT(mesh.SignedVolume(0), 0.0);
  EXPECT_EQ(mesh.FixOrientation(), 2u);
  for (size_t t = 0; t < mesh.NumTets(); ++t) EXPECT_GT(mesh.SignedVolume(t), 0.0);
}

TEST(TetMeshTest, VertexNormalsPointOutwards) {
  const TetMesh mesh = MakeTetBeam(2, 2, 2, EigenVector3(1.0, 1.0, 1.0));
  const EigenMatrix3x n = VertexNormals(mesh.vertices, mesh.BoundaryTriangles());
  const EigenVector3 centre(0.5, 0.5, 0.5);
  for (unsigned int i = 0; i < mesh.NumVertices(); ++i) {
    const EigenVector3 out = mesh.vertices.row(i).transpose() - centre;
    if (out.norm() < 1e-9) {
      EXPECT_EQ(n.row(i).norm(), 0.0) << "interior vertex has a normal";
    } else {
      EXPECT_NEAR(n.row(i).norm(), 1.0, 1e-12);
      EXPECT_GT(n.row(i).dot(out.transpose()), 0.0) << "vertex " << i;
    }
  }
}

TEST(TetMeshTest, TetGenRoundTrip) {
  const TetMesh mesh = MakeTetBeam(3, 2, 1, EigenVector3(1.0, 0.5, 0.25), EigenVector3(-1.0, 2.0, 0.5));
  const std::string path = TempPath("pd_roundtrip");
  ASSERT_TRUE(SaveTetGen(path, mesh));

  TetMesh loaded;
  std::string error;
  ASSERT_TRUE(LoadTetGen(path, &loaded, &error)) << error;
  EXPECT_EQ(loaded.NumVertices(), mesh.NumVertices());
  ASSERT_EQ(loaded.NumTets(), mesh.NumTets());
  EXPECT_EQ((loaded.vertices - mesh.vertices).cwiseAbs().maxCoeff(), 0.0);
  for (size_t t = 0; t < mesh.NumTets(); ++t) EXPECT_EQ(loaded.tets[t], mesh.tets[t]);
  std::remove((path + ".node").c_str());
  std::remove((path + ".ele").c_str());
}

// TetGen output as written by the tool: 1-based ids, comments, attributes,
// boundary markers, and an inverted tet that must be re-oriented.
TEST(TetMeshTest, LoadsOneBasedTetGenWithComments) {
  const std::string path = TempPath("pd_onebased");
  {
    std::ofstream node(path + ".node");
    node << "# Generated by tetgen\n"
            "5  3  1  1\n"
            "   1    0 0 0   7.0  1\n"
            "   2    1 0 0   7.0  1\n"
            "   3    0 1 0   7.0  1   # trailing comment\n"
            "\n"
            "   4    0 0 1   7.0  1\n"
            "   5    1 1 1   7.0  0\n";
    std::ofstream ele(path + ".ele");
    ele << "2  4  0\n"
           "  1   1 2 3 4\n"
           "  2   2 3 4 5\n";  // Negatively oriented as written.
  }
  TetMesh mesh;
  std::string error;
  ASSERT_TRUE(LoadTetGen(path, &mesh, &error)) << error;
  ASSERT_EQ(mesh.NumVertices(), 5u);
  ASSERT_EQ(mesh.NumTets(), 2u);
  EXPECT_EQ(mesh.tets[0][0], 0u);
  EXPECT_NEAR(mesh.SignedVolume(0), 1.0 / 6.0, 1e-15);
  EXPECT_NEAR(mesh.SignedVolume(1), 1.0 / 3.0, 1e-15);
  EXPECT_TRUE(IsClosedAndConsistent(mesh.BoundaryTriangles()));
  std::remove((path + ".node").c_str());
  std::remove((path + ".ele").c_str());
}

TEST(TetMeshTest, LoadReportsErrors) {
  TetMesh mesh;
  std::string error;
  EXPECT_FALSE(LoadTetGen(TempPath("pd_does_not_exist"), &mesh, &error));
  EXPECT_NE(error.find(".node"), std::string::npos);

  const std::string path = TempPath("pd_bad_ids");
  {
    std::ofstream node(path + ".node");
    node << "4 3 0 0\n0 0 0 0\n1 1 0 0\n2 0 1 0\n3 0 0 1\n";
    std::ofstream ele(path + ".ele");
    ele << "1 4 0\n0 0 1 2 9\n";
  }
  EXPECT_FALSE(LoadTetGen(path, &mesh, &error));
  EXPECT_NE(error.find("bad vertex id"), std::string::npos) << error;
  std::remove((path + ".node").c_str());
  std::remove((path + ".ele").c_str());
}

// --- TetScene on the solver ------------------------------------------------------

TEST(TetSceneTest, BuildsClampedBeam) {
  TetScene scene;
  PDSolver sim;
  scene.Build(&sim);
  ASSERT_TRUE(sim.IsReady());
  EXPECT_EQ(sim.NumVertices(), scene.Mesh().NumVertices());
  EXPECT_EQ(sim.NumConstraints(), scene.Mesh().Edges().size());
  // One 3x3 grid of vertices at the clamped end (10 x 2 x 2 cells).
  EXPECT_EQ(sim.NumFixed(), 9u);
  EXPECT_NEAR(sim.GetMasses().sum(), 1000.0 * 0.04, 1e-9);
  EXPECT_EQ(scene.SurfaceIndices().size() % 3, 0u);
  EXPECT_TRUE(scene.LineIndices().empty());
  scene.m_show_wireframe = true;
  EXPECT_EQ(scene.LineIndices().size(), 2 * sim.NumConstraints());
}

TEST(TetSceneTest, RestPoseIsStableWithoutGravity) {
  TetScene scene;
  PDSolver sim;
  sim.m_gravity_constant = 0.0;
  scene.Build(&sim);
  const EigenMatrix3x initial = sim.GetPositions();
  EXPECT_NEAR(sim.ConstraintEnergy(initial), 0.0, 1e-18);
  for (int i = 0; i < 30; ++i) sim.Step();
  EXPECT_LT((sim.GetPositions() - initial).cwiseAbs().maxCoeff(), 1e-9);
}

// The cantilever sags under gravity: the free end drops, the clamped end stays,
// and the tet mesh does not invert.
TEST(TetSceneTest, CantileverSagsWithoutInverting) {
  TetScene scene;
  PDSolver sim;
  scene.Build(&sim);
  const EigenMatrix3x initial = sim.GetPositions();
  for (int i = 0; i < 300; ++i) sim.Step();
  const EigenMatrix3x& x = sim.GetPositions();

  ScalarType tip_drop = 0.0;
  for (unsigned int i = 0; i < sim.NumVertices(); ++i) {
    if (sim.IsFixed(i)) {
      EXPECT_EQ((x.row(i) - initial.row(i)).norm(), 0.0);
    } else if (initial(i, 0) > 0.49) {
      tip_drop = std::max(tip_drop, initial(i, 1) - x(i, 1));
    }
  }
  EXPECT_GT(tip_drop, 1e-3) << "beam does not sag";
  EXPECT_LT(tip_drop, 1.0) << "beam collapsed";
  for (size_t t = 0; t < scene.Mesh().NumTets(); ++t) {
    EXPECT_GT(scene.Mesh().SignedVolume(t, x), 0.0) << "tet " << t << " inverted";
  }
}
