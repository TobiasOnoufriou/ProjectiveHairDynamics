#pragma once
#include <vector>

#include "GLCore.h"
#include <GLCoreUtils.h>

#include "GlmEigen.h"

// Draws a deforming triangle surface (e.g. the boundary of a tet mesh) with
// per-vertex normals recomputed from the current positions every upload.
class SurfaceRenderer {
 public:
  SurfaceRenderer();
  ~SurfaceRenderer();

  // positions: N x 3 simulation vertices; indices: triangle triples into them.
  void Upload(const EigenMatrix3x& positions, const std::vector<unsigned int>& indices);
  void Draw(const GLCore::Utils::PerspectiveCameraController& camera);
  void SetColor(const glm::vec3& color) { m_color = color; }

 private:
  GLuint m_vao = 0, m_vbo = 0, m_ebo = 0;
  std::vector<unsigned int> m_indices;  // Last uploaded, to skip re-uploading the EBO.
  std::vector<float> m_vertex_data;     // Interleaved position, normal.
  GLCore::Utils::Shader* m_shader;
  glm::vec3 m_color = glm::vec3(0.80f, 0.33f, 0.30f);
};
