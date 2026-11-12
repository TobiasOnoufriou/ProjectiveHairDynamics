// Copyright 2022 Tobias Onoufriou
// Conversions between the solver's Eigen types (PD-Core, no glm) and the glm
// types used by the renderer and camera code of the sandbox.
#pragma once
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "math_headers.h"

glm::vec3 Eigen2GLM(const EigenVector3& eigen_vector);
EigenVector3 GLM2Eigen(const glm::vec3& glm_vector);
