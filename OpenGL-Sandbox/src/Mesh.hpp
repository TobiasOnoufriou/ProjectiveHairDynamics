// Copyright 2022 Tobias Onoufriou

#pragma once
#include "GLCore.h"
#include <GLCoreUtils.h>


#include "Hittable.h"
#include "BoxCollider.h"

#include <glm/gtc/type_ptr.hpp>

#include <string>
#include <vector>
#include <cmath>


struct Vertex {
  glm::vec3 position;
  glm::vec3 normal;
  glm::vec2 tex_coords;
  glm::vec3 tangent;
  glm::vec3 bitangent;
};

struct Texture {
  unsigned int id;
  std::string type;
  std::string path;
};

// Mesh would inherit a collider and material

class Mesh: public BoxCollider {
public:
  Mesh(
    std::vector<Vertex> vertices,
    std::vector<GLuint> indices,
    std::vector<Texture> textures
  );
  ~Mesh();

  void Draw(const GLCore::Utils::PerspectiveCameraController& camera);
  void ChangeVertexPosition(glm::vec3 pos, int index);
  // Ray vs axis-aligned bounding box (slab test). A zero direction component
  // gives +-inf slab distances, which fmin/fmax handle correctly.
  bool intersectionTest(glm::vec3 origin, glm::vec3 dir) {
    float tmin = 0.0f, tmax = INFINITY;

    for (int i = 0; i < 3; ++i) {
      const float inv_dir = 1.0f / dir[i];
      const float t1 = (pMin[i] - origin[i]) * inv_dir;
      const float t2 = (pMax[i] - origin[i]) * inv_dir;

      tmin = std::fmax(tmin, std::fmin(t1, t2));
      tmax = std::fmin(tmax, std::fmax(t1, t2));
    }

    return tmin <= tmax;
  }
  void TranslateShape(glm::vec3 translation);
  void SetupMesh();
  void ComputeNormal();
  virtual void CalculateBoundingVolume() override;
  virtual void DrawBoundingBox(const GLCore::Utils::PerspectiveCameraController& cam) override;


  bool GetCollider() { return collider; }

 public:
  std::vector<Vertex> vertices_;
  std::vector<GLuint> indices_;
  std::vector<Texture> textures_;

  //Scale and Position
  glm::vec3 Scale_, Position_;
  glm::mat4 model;
  
  // Mesh/Material Shader.
  GLCore::Utils::Shader* shader_;
  GLCore::Utils::Shader* debug_shader_;

  GLuint vbo_, ebo_, vao_;
};
