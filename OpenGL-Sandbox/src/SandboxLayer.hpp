#pragma once

#include <GLCore.h>
#include <GLCoreUtils.h>

#include "ModelLoader.hpp"
#include "LineTools.h"
#include "Simulation.h"
#include "StrandRenderer.h"


class SandboxLayer : public GLCore::Layer
{
public:
	SandboxLayer(uint32_t w, uint32_t h);
	virtual ~SandboxLayer();

	virtual void OnAttach() override;
	virtual void OnDetach() override;
	virtual void OnEvent(GLCore::Event& event) override;
	virtual void OnUpdate(GLCore::Timestep ts) override;
	virtual void OnImGuiRender() override;
private:
  std::unique_ptr<ModelLoader> model_loader_;
  std::unique_ptr<Mesh> mesh_;
  std::unique_ptr<Mesh> plane_;
  std::unique_ptr<LineTools> line_;

  std::unique_ptr<Simulation> m_sim;
  std::unique_ptr<StrandRenderer> m_strand_renderer;
  bool m_run_sim = true;

  const aiScene* kScene; // Current Scene.
  
  glm::vec3 velocity = { 0.1f, 0.0f, -0.1f };

  GLCore::Utils::Shader* shader_;
  GLCore::Utils::PerspectiveCameraController camera_controller_;

  uint32_t width, height;
};