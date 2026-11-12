#pragma once
#include "GlmEigen.h"

class BoxCollider {
public:
  BoxCollider() {
    collider = true;
  }
  BoxCollider(glm::vec3 pMin, glm::vec3 pMax, glm::vec3 c) {
    this->c = c;
    this->pMin = pMin;
    this->pMax = pMax;
  }
  ~BoxCollider() {

  }

  virtual void DrawBoundingBox(const GLCore::Utils::PerspectiveCameraController& cam) = 0;

  virtual void CalculateBoundingVolume() = 0;
  

  //for now public eventually move to private.
protected:
  glm::vec3 pMin, pMax;
  glm::vec3 c;
  bool collider; // turns on and off the collider
};

// TODO: Add OBB Bounding Box;