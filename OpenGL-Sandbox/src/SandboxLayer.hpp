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
  // Returns the particle whose screen projection is closest to the mouse
  // (within kPickRadiusPx), or -1. Writes its NDC depth to out_ndc_depth.
  int PickParticle(float mouse_x, float mouse_y, float* out_ndc_depth) const;
  // Unprojects a mouse position at a fixed NDC depth back into world space.
  glm::vec3 MouseToWorld(float mouse_x, float mouse_y, float ndc_depth) const;

  void DrawSimulationTab();
  void DrawCameraTab();

  static constexpr float kPickRadiusPx = 15.0f;
  float m_drag_ndc_depth = 0.0f;  // Keeps the dragged particle at a constant depth from the camera.

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