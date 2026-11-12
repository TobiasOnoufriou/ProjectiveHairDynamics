// Copyright 2022 Tobias Onoufriou
#include "ModelLoader.hpp"

#include "GlmEigen.h"


ModelLoader::ModelLoader() {

}

ModelLoader::ModelLoader(
  std::string const& path,
  bool gamma = false):
  gamme_correction_(gamma) {
    LoadModel(path);
}

ModelLoader::~ModelLoader() {

}

std::unique_ptr<Mesh>& ModelLoader::LoadModel(std::string const& path) {
  kScene = aiImportFile(path.c_str(),
    aiProcess_Triangulate | aiProcess_FlipUVs | aiProcess_GenSmoothNormals | aiProcess_CalcTangentSpace
    );

  if (!kScene ||
    kScene->mFlags &
    AI_SCENE_FLAGS_INCOMPLETE ||
    !kScene->mRootNode) {
    std::cout << "ERROR::ASSIMP:: " << std::endl;
    return std::unique_ptr<Mesh>(nullptr);
  }
  
  directory_ = path.substr(0, path.find_last_of('/'));

  this->ProcessNode(kScene->mRootNode, kScene);


  if (meshes_->GetCollider()) {
    //CenterAndScale(&vertices[0].position, sizeof(Vertex), vertices.size(), 0);
    meshes_->CalculateBoundingVolume();
  }

  return meshes_;
}

void ModelLoader::ProcessNode(aiNode* node, const aiScene* kScene) {
  if (node->mMeshes == nullptr) {
    aiMesh* mesh = kScene->mMeshes[0];
    ProcessMesh(mesh, kScene);
  }
  
  
  else {
    for (unsigned int i = 0; i < node->mNumChildren; i++) {
      //unsigned int index = node->mMeshes[i];
      aiMesh* mesh = kScene->mMeshes[node->mMeshes[i]];
      ProcessMesh(mesh, kScene);
    }
  } 
  
  if (node->mNumChildren > 0) {
    for (unsigned int i = 0; i < node->mNumChildren; i++) {
      this->ProcessNode(node->mChildren[i], kScene);
    }
  }
}



void ModelLoader::ProcessMesh(aiMesh* mesh, const aiScene* kScene) {
  std::vector<Vertex> vertices;
  std::vector<unsigned int> indices;
  std::vector<Texture> textures;


  for (unsigned int i = 0; i < mesh->mNumVertices; i++) {
    Vertex vertex;
    glm::vec3 vector;

    //Use Eigen to assign current positions.
    // positions
    vector.x = mesh->mVertices[i].x;
    vector.y = mesh->mVertices[i].y;
    vector.z = mesh->mVertices[i].z;
    vertex.position = vector;
    // normals
    if (mesh->HasNormals()) {
      vector.x = mesh->mNormals[i].x;
      vector.y = mesh->mNormals[i].y;
      vector.z = mesh->mNormals[i].z;
      vertex.normal = vector;
    }
    if (mesh->mTextureCoords[0]) {
      glm::vec2 vec;
      // a vertex can contain up to 8 different texture coordinates. We thus make the assumption that we won't 
      // use models where a vertex can have multiple texture coordinates so we always take the first set (0).
      vec.x = mesh->mTextureCoords[0][i].x;
      vec.y = mesh->mTextureCoords[0][i].y;
      vertex.tex_coords = vec;
      // tangent
      if (mesh->HasTangentsAndBitangents()) {
        vector.x = mesh->mTangents[i].x;
        vector.y = mesh->mTangents[i].y;
        vector.z = mesh->mTangents[i].z;
        vertex.tangent = vector;
        
        // bitangent
        vector.x = mesh->mBitangents[i].x;
        vector.y = mesh->mBitangents[i].y;
        vector.z = mesh->mBitangents[i].z;
        vertex.bitangent = vector;
      }
    } else {
      vertex.tex_coords = glm::vec2(0.0f, 0.0f);
    }
    vertices.push_back(vertex);
  }
  for (unsigned int i = 0; i < mesh->mNumFaces; i++) {
    aiFace face = mesh->mFaces[i];
    // retrieve all indices of the face and store them in the indices vector
    for (unsigned int j = 0; j < face.mNumIndices; j++)
      indices.push_back(face.mIndices[j]);
  }
  
  
  aiMaterial* material = kScene->mMaterials[mesh->mMaterialIndex];


  std::vector<Texture> diffuseMaps = LoadMaterialTextures(material, aiTextureType_DIFFUSE, "texture_diffuse");
  textures.insert(textures.end(), diffuseMaps.begin(), diffuseMaps.end());
  std::vector<Texture> specularMaps = LoadMaterialTextures(material, aiTextureType_SPECULAR, "texture_specular");
  textures.insert(textures.end(), specularMaps.begin(), specularMaps.end());
  std::vector<Texture> normalMaps = LoadMaterialTextures(material, aiTextureType_HEIGHT, "texture_normal");
  textures.insert(textures.end(), normalMaps.begin(), normalMaps.end());
  std::vector<Texture> heightMaps = LoadMaterialTextures(material, aiTextureType_AMBIENT, "texture_height");
  textures.insert(textures.end(), heightMaps.begin(), heightMaps.end());

  meshes_ = std::make_unique<Mesh>(vertices, indices, textures);
}

template<typename Vec>
std::pair<Vec, Vec> ModelLoader::GetExtents(const Vec* pts, size_t stride, size_t count){
  unsigned char* base = (unsigned char*)pts;
  Vec pmin(*(Vec*)base);
  Vec pmax(*(Vec*)base);

  for (size_t i = 0; i < count; ++i, base += stride) {
    const Vec& pt = *(Vec*)base;
    pmin = glm::min(pmin, pt);
    pmax = glm::max(pmax, pt);
  }
  return std::make_paiz(pmin, pmax);
}

template<typename Vec>
void ModelLoader::CenterAndScale(Vec* pts, size_t stride, size_t count, const typename Vec::value_type& size) {
  typedef typename Vec::value_type Scalar;

  std::pair<Vec, Vec> exts = GetExtents(pts, stride, count);

  const Vec center = (exts.first * Scalar(0.5f)) + (exts.second * Scalar(0.5f));

  meshes_->c = center;
  meshes_->pMin = exts.first;
  meshes_->pMax = exts.second;
}



unsigned int ModelLoader::TextureFromFile(std::string const& path, const std::string& directory, bool gamma) {
  std::string filename = std::string(path);
  filename = directory + '/' + filename;

  unsigned int textureID;
  glGenTextures(1, &textureID);

  int width, height, nrComponents;
  unsigned char* data = stbi_load(filename.c_str(), &width, &height, &nrComponents, 0);
  if (data) {
    GLenum format;
    if (nrComponents == 1)
      format = GL_RED;
    else if (nrComponents == 3)
      format = GL_RGB;
    else if (nrComponents == 4)
      format = GL_RGBA;

    glBindTexture(GL_TEXTURE_2D, textureID);
    glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
    glGenerateMipmap(GL_TEXTURE_2D);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    stbi_image_free(data);
  }
  else {
    std::cout << "Texture failed to load at path: " << path << std::endl;
    stbi_image_free(data);
  }

  return textureID;
}

std::vector<Texture> ModelLoader::LoadMaterialTextures(aiMaterial* mat, aiTextureType type, std::string typeName) {
  std::vector<Texture> textures;
  for (unsigned int i = 0; i < mat->GetTextureCount(type); i++) {
    aiString str;
    mat->GetTexture(type, i, &str);
    bool skip = false;
    for (unsigned int j = 0; j < textures_loaded_.size(); j++) {
      if (std::strcmp(textures_loaded_[j].path.data(), str.C_Str()) == 0) {
        textures.push_back(textures_loaded_[j]);
        skip = true;
        break;
      }
    }
    if (!skip) {
      Texture texture;
      texture.id = TextureFromFile(str.C_Str(), this->directory_, false);
      texture.type = typeName;
      texture.path = str.C_Str();
      textures.push_back(texture);
      textures_loaded_.push_back(texture);
    }
  }
  return textures;
}


