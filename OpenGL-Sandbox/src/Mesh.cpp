// Copyright 2022 Tobias Onoufriou

#include "Mesh.hpp"
#include "math_headers.h"


Mesh::Mesh(std::vector<Vertex> vertices,
  std::vector<GLuint> indices,
  std::vector<Texture> textures) {
  vertices_ = vertices;
  indices_ = indices;
  textures_ = textures;

  shader_ = GLCore::Utils::Shader::FromGLSLTextFiles(
    "src/shaders/vertex.vs",
    "src/shaders/fragment.fs"
  );

  debug_shader_ = GLCore::Utils::Shader::FromGLSLTextFiles(
    "src/shaders/debug_vertex.vs",
    "src/shaders/debug_fragment.fs"
  );


  model = glm::mat4(1.0f);
  Scale_ = glm::vec3(1.0f, 1.0f, 1.0f);
  Position_ = glm::vec3(0.0f, 0.0f, 1.0f);
  collider = true;
}

Mesh::~Mesh() {

}

void Mesh::ComputeNormal() {
  glm::vec3 zero(0.0);
  for (std::vector<Vertex>::iterator n = vertices_.begin(); n != vertices_.end(); ++n) {

  }
}


void Mesh::ChangeVertexPosition(glm::vec3 pos, int index) {
  this->vertices_.at(0).position *= pos;
}

void Mesh::TranslateShape(glm::vec3 translation) {
  Position_ += translation;
  model = glm::translate(model, Position_);
  pMax = pMax + Position_;
  pMin = pMin + Position_;


  std::cout << "Max: " << std::endl;
  std::cout << pMax.x << std::endl;
  std::cout << pMax.y << std::endl;
  std::cout << pMax.z << std::endl;

  std::cout << "Min: " << std::endl;
  std::cout << pMin.x << std::endl;
  std::cout << pMin.y << std::endl;
  std::cout << pMin.z << std::endl;
}

void Mesh::CalculateBoundingVolume() {
  pMin = this->vertices_[0].position;
  pMax = this->vertices_[0].position;

  for (int i = 1; i < this->vertices_.size(); ++i) {
    pMin = glm::min(pMin, this->vertices_[i].position);
    pMax = glm::max(pMax, this->vertices_[i].position);
  }
}


void Mesh::DrawBoundingBox(const GLCore::Utils::PerspectiveCameraController& cam) {
  if (vertices_.size() == 0)
    return;
  glUseProgram(debug_shader_->GetRendererID());

  GLfloat vertices[] = {
  -0.5, -0.5, -0.5, 1.0,
   0.5, -0.5, -0.5, 1.0,
   0.5,  0.5, -0.5, 1.0,
  -0.5,  0.5, -0.5, 1.0,
  -0.5, -0.5,  0.5, 1.0,
   0.5, -0.5,  0.5, 1.0,
   0.5,  0.5,  0.5, 1.0,
  -0.5,  0.5,  0.5, 1.0,
  };

  glm::mat4 ve;

  GLuint vbo_vertices;
  glGenBuffers(1, &vbo_vertices);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_vertices);
  glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  GLushort elements[] = {
  0, 1, 2, 3,
  4, 5, 6, 7,
  0, 4, 1, 5, 2, 6, 3, 7
  };

  int location = glGetUniformLocation(debug_shader_->GetRendererID(), "u_ViewProjection");
  glUniformMatrix4fv(location, 1, GL_FALSE, glm::value_ptr(cam.GetCamera().GetViewProjectionMatrix()));

   //Add a move and rotate function within the mesh class.

  glm::vec3 size = glm::vec3(pMax.x - pMin.x, pMax.y - pMin.y, pMax.z - pMin.z);
  glm::vec3 center = glm::vec3((pMin.x + pMax.x), (pMin.y + pMax.y) , (pMin.z + pMax.z));
  glm::mat4 transform = glm::scale(glm::mat4(1), size);;

  ve = model * transform;


  location = glGetUniformLocation(debug_shader_->GetRendererID(), "model");
  glUniformMatrix4fv(location, 1, GL_FALSE, glm::value_ptr(ve));

  GLuint ibo_elements;
  glGenBuffers(1, &ibo_elements);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_elements);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(elements), elements, GL_STATIC_DRAW);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

  glBindBuffer(GL_ARRAY_BUFFER, vbo_vertices);


  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(GLfloat)*4, (void*)0);


  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_elements);
  glDrawElements(GL_LINE_LOOP, 4, GL_UNSIGNED_SHORT, 0);
  glDrawElements(GL_LINE_LOOP, 4, GL_UNSIGNED_SHORT, (GLvoid*)(4 * sizeof(GLushort)));
  glDrawElements(GL_LINES, 8, GL_UNSIGNED_SHORT, (GLvoid*)(8 * sizeof(GLushort)));
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

  glBindBuffer(GL_ARRAY_BUFFER, 0);

  glDeleteBuffers(1, &vbo_vertices);
  glDeleteBuffers(1, &ibo_elements);
}

void Mesh::SetupMesh() {
  if(!glIsBuffer(vao_))
    glGenVertexArrays(1, &vao_);
  if (!glIsBuffer(vbo_))
    glGenBuffers(1, &vbo_);
  if (!glIsBuffer(ebo_))
    glGenBuffers(1, &ebo_);

  glBindVertexArray(vao_);


  //std::cout << vertices_.at(1).position.x << std::endl;
  //for (unsigned int i = 0; i < vertices_.size(); ++i) {
  //}

  glBindBuffer(GL_ARRAY_BUFFER, vbo_);
  glBufferData(GL_ARRAY_BUFFER, vertices_.size() * sizeof(Vertex), &vertices_[0], GL_DYNAMIC_DRAW);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo_);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices_.size() * sizeof(unsigned int), &indices_[0], GL_STATIC_DRAW);

  // vertex position
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
  // vertex normals
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));
  // vertex texture coords
  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, tex_coords));
  // vertex tangent
  glEnableVertexAttribArray(3);
  glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, tangent));
  // vertex bitangent
  glEnableVertexAttribArray(4);
  glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, bitangent));
  
  glBindVertexArray(0);
}

// Draw will be added to the update function.
// Implementation of Shader is needed.
// Could be defined as a Template
void Mesh::Draw(const GLCore::Utils::PerspectiveCameraController& camera) {


  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glUseProgram(shader_->GetRendererID());

  int location = glGetUniformLocation(shader_->GetRendererID(), "u_ViewProjection");
  glUniformMatrix4fv(location, 1, GL_FALSE, glm::value_ptr(camera.GetCamera().GetViewProjectionMatrix()));

  // Add a move and rotate function within the mesh class.
  //model = glm::scale(model, Scale_);

  location = glGetUniformLocation(shader_->GetRendererID(), "model");
  glUniformMatrix4fv(location, 1, GL_FALSE, glm::value_ptr(model));
  

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, 0);
  location = glGetUniformLocation(shader_->GetRendererID(), "t0");
  glUniform1i(location, 0);



  // get uniform location.
  // bind appropriate textures
  unsigned int diffuseNr = 1;
  unsigned int specularNr = 1;
  unsigned int normalNr = 1;
  unsigned int heightNr = 1;

  SetupMesh();

  for (unsigned int i = 0; i < textures_.size(); i++)
  {
    glActiveTexture(GL_TEXTURE0 + i);
    std::string number;
    std::string name = textures_[i].type;
    if (name == "texture_diffuse")
      number = std::to_string(diffuseNr++);
    else if (name == "texture_specular")
      number = std::to_string(specularNr++);
    else if (name == "texture_normal")
      number = std::to_string(normalNr++);
    else if (name == "texture_height")
      number = std::to_string(heightNr++);

    glUniform1i(glGetUniformLocation(shader_->GetRendererID(), (name + number).c_str()), i);
    glBindTexture(GL_TEXTURE_2D, textures_[i].id);
  }

  // draw mesh
  glBindVertexArray(vao_);
  glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, 0);
  
  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);
  glDisableVertexAttribArray(3);

  glBindBuffer(GL_ARRAY_BUFFER, 0);

  // always good practice to set everything back to defaults once configured.
  glActiveTexture(GL_TEXTURE0);

}


