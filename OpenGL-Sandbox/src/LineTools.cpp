#include "LineTools.h"


LineTools::LineTools() {
  std::cout << "Debug tool init" << std::endl;
  shader_ = GLCore::Utils::Shader::FromGLSLTextFiles(
    "src/shaders/debug_vertex.vs",
    "src/shaders/debug_fragment.fs"
  );
  model_ = glm::mat4(1.0f);
  this->InitLine(glm::vec3(0.0f, 0.0f, 0.0f),glm::vec3(0.0f, 0.0f, 10.0f), glm::mat4(1.0f));

}

LineTools::LineTools(glm::vec3 cam_pos) {
  std::cout << "Debug tool init" << std::endl;
  shader_ = GLCore::Utils::Shader::FromGLSLTextFiles(
    "src/shaders/debug_vertex.vs",
    "src/shaders/debug_fragment.fs"
  );
  model_ = glm::mat4(1.0f);
  glm::vec3 cam_pos_offset = cam_pos;
  cam_pos_offset.x += 1.0f;
  this->InitLine(cam_pos, cam_pos_offset, glm::mat4(1.0f));
}

void LineTools::InitLine(glm::vec3 min_vec,glm::vec3 max_vec, glm::mat4 model) {
  //Min line vector.
  line_vertices[0] = min_vec.x;
  line_vertices[1] = min_vec.y;
  line_vertices[2] = min_vec.z;

  //Max line vector. Direction need to do something with the camera to change this below.
  line_vertices[3] = max_vec.x;
  line_vertices[4] = max_vec.y;
  line_vertices[5] = max_vec.z;

  model_ = model;
}


LineTools::~LineTools() {

}



void LineTools::DrawDebugLine(const GLCore::Utils::PerspectiveCameraController& camera) {
  glUseProgram(shader_->GetRendererID());

  glGenVertexArrays(1, &vao_);
  glGenBuffers(1, &vbo_);


  glBindVertexArray(vao_);

  glBindBuffer(GL_ARRAY_BUFFER, vbo_);
  glBufferData(GL_ARRAY_BUFFER, sizeof(line_vertices), line_vertices, GL_STATIC_DRAW);


  GLuint location = glGetUniformLocation(shader_->GetRendererID(), "model");
  glad_glUniformMatrix4fv(location, 1, GL_FALSE, glm::value_ptr(model_));

  location = glGetUniformLocation(shader_->GetRendererID(), "u_ViewProjection");
  glUniformMatrix4fv(location, 1, GL_FALSE, glm::value_ptr(camera.GetCamera().GetViewProjectionMatrix()));

  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
  

  //Pass in shader to this function.

  glDrawArrays(GL_LINES, 0, 2);

  glDisableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
}