#pragma once
#include "GLCore.h"
#include <GLCoreUtils.h>

class LineTools {
public:
  LineTools();
  LineTools(glm::vec3 cam_pos);
  ~LineTools();
  void InitLine(glm::vec3 min_vec,glm::vec3 max_vec, glm::mat4 model);
  void DrawDebugLine(const GLCore::Utils::PerspectiveCameraController& camera);
private:
  GLuint vbo_, vao_;
  glm::mat4 model_;
  GLCore::Utils::Shader* shader_;

  float line_vertices[6];
};