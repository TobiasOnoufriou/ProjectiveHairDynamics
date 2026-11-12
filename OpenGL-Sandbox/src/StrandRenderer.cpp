#include "StrandRenderer.h"

StrandRenderer::StrandRenderer() {
  shader_ = GLCore::Utils::Shader::FromGLSLTextFiles(
    "src/shaders/debug_vertex.vs",
    "src/shaders/debug_fragment.fs"
  );
  glGenVertexArrays(1, &vao_);
  glGenBuffers(1, &vbo_);
  glGenBuffers(1, &ebo_);

  glBindVertexArray(vao_);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo_);
  glBindVertexArray(0);
}

StrandRenderer::~StrandRenderer() {
  glDeleteBuffers(1, &vbo_);
  glDeleteBuffers(1, &ebo_);
  glDeleteVertexArrays(1, &vao_);
  delete shader_;
}

void StrandRenderer::Upload(const EigenMatrix3x& positions, const std::vector<unsigned int>& indices) {
  vertex_data_.resize(positions.rows() * 3);
  for (Eigen::Index i = 0; i < positions.rows(); ++i) {
    vertex_data_[3 * i + 0] = static_cast<float>(positions(i, 0));
    vertex_data_[3 * i + 1] = static_cast<float>(positions(i, 1));
    vertex_data_[3 * i + 2] = static_cast<float>(positions(i, 2));
  }

  glBindVertexArray(vao_);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_);
  glBufferData(GL_ARRAY_BUFFER, vertex_data_.size() * sizeof(float), vertex_data_.data(), GL_DYNAMIC_DRAW);
  if (indices != indices_) {
    indices_ = indices;
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);
    index_count_ = static_cast<GLsizei>(indices.size());
  }
  glBindVertexArray(0);
}

void StrandRenderer::Draw(const GLCore::Utils::PerspectiveCameraController& camera) {
  if (index_count_ == 0) return;
  glUseProgram(shader_->GetRendererID());

  glm::mat4 model(1.0f);
  GLint location = glGetUniformLocation(shader_->GetRendererID(), "model");
  glUniformMatrix4fv(location, 1, GL_FALSE, glm::value_ptr(model));
  location = glGetUniformLocation(shader_->GetRendererID(), "u_ViewProjection");
  glUniformMatrix4fv(location, 1, GL_FALSE, glm::value_ptr(camera.GetCamera().GetViewProjectionMatrix()));

  glBindVertexArray(vao_);
  glDrawElements(GL_LINES, index_count_, GL_UNSIGNED_INT, (GLvoid*)0);
  glBindVertexArray(0);
}
