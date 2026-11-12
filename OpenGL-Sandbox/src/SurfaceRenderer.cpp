#include "SurfaceRenderer.h"

#include <vector>

#include "TetMesh.h"

SurfaceRenderer::SurfaceRenderer() {
  m_shader = GLCore::Utils::Shader::FromGLSLTextFiles("src/shaders/surface.vs", "src/shaders/surface.fs");
  glGenVertexArrays(1, &m_vao);
  glGenBuffers(1, &m_vbo);
  glGenBuffers(1, &m_ebo);

  glBindVertexArray(m_vao);
  glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
  const GLsizei stride = 6 * sizeof(GLfloat);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<GLvoid*>(0));
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<GLvoid*>(3 * sizeof(GLfloat)));
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ebo);
  glBindVertexArray(0);
}

SurfaceRenderer::~SurfaceRenderer() {
  glDeleteBuffers(1, &m_vbo);
  glDeleteBuffers(1, &m_ebo);
  glDeleteVertexArrays(1, &m_vao);
  delete m_shader;
}

void SurfaceRenderer::Upload(const EigenMatrix3x& positions, const std::vector<unsigned int>& indices) {
  std::vector<TetMesh::Triangle> triangles(indices.size() / 3);
  for (size_t t = 0; t < triangles.size(); ++t) {
    triangles[t] = {indices[3 * t], indices[3 * t + 1], indices[3 * t + 2]};
  }
  const EigenMatrix3x normals = VertexNormals(positions, triangles);

  m_vertex_data.resize(positions.rows() * 6);
  for (Eigen::Index i = 0; i < positions.rows(); ++i) {
    for (int c = 0; c < 3; ++c) {
      m_vertex_data[6 * i + c] = static_cast<float>(positions(i, c));
      m_vertex_data[6 * i + 3 + c] = static_cast<float>(normals(i, c));
    }
  }

  glBindVertexArray(m_vao);
  glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
  glBufferData(GL_ARRAY_BUFFER, m_vertex_data.size() * sizeof(float), m_vertex_data.data(), GL_DYNAMIC_DRAW);
  if (indices != m_indices) {
    m_indices = indices;
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);
  }
  glBindVertexArray(0);
}

void SurfaceRenderer::Draw(const GLCore::Utils::PerspectiveCameraController& camera) {
  if (m_indices.empty()) return;
  const GLuint program = m_shader->GetRendererID();
  glUseProgram(program);
  glUniformMatrix4fv(glGetUniformLocation(program, "u_ViewProjection"), 1, GL_FALSE,
                     glm::value_ptr(camera.GetCamera().GetViewProjectionMatrix()));
  const glm::vec3 eye = camera.GetCamera().GetPosition();
  glUniform3f(glGetUniformLocation(program, "u_CameraPosition"), eye.x, eye.y, eye.z);
  glUniform3f(glGetUniformLocation(program, "u_Color"), m_color.x, m_color.y, m_color.z);

  // Push the surface back slightly so a wireframe drawn on top stays visible.
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1.0f, 1.0f);
  glBindVertexArray(m_vao);
  glDrawElements(GL_TRIANGLES, static_cast<GLsizei>(m_indices.size()), GL_UNSIGNED_INT, nullptr);
  glBindVertexArray(0);
  glDisable(GL_POLYGON_OFFSET_FILL);
}
