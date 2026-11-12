#pragma once
#include <vector>
#include "GLCore.h"
#include <GLCoreUtils.h>
#include "math_headers.h"

// Draws simulated hair strands as GL_LINES. The buffers are created once and
// the vertex data is re-uploaded every frame.
class StrandRenderer {
 public:
  StrandRenderer();
  ~StrandRenderer();

  // positions: N x 3 particle positions, indices: pairs of particle indices.
  void Upload(const EigenMatrix3x& positions, const std::vector<unsigned int>& indices);
  void Draw(const GLCore::Utils::PerspectiveCameraController& camera);

 private:
  GLuint vao_ = 0, vbo_ = 0, ebo_ = 0;
  GLsizei index_count_ = 0;
  std::vector<float> vertex_data_;
  GLCore::Utils::Shader* shader_;
};
