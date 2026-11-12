#include "SandboxLayer.hpp"
#include "Mesh.hpp"
#include "glm/gtx/rotate_vector.hpp"
#include "GLCore/Core/Input.h"
#include "GLCore/Core/MouseButtonCodes.h"

#include <cmath>

using namespace GLCore;
using namespace GLCore::Utils;

SandboxLayer::SandboxLayer(uint32_t w, uint32_t h)
  : camera_controller_(h > 0 ? (float)w / (float)h : 16.0f / 9.0f), width(w), height(h)
{
  camera_controller_.SetViewportSize(w, h);
}

SandboxLayer::~SandboxLayer()
{
}

float angleBetween(
  glm::vec3 a,
  glm::vec3 b,
  glm::vec3 origin
) {
  glm::vec3 da = glm::normalize(a - origin);
  glm::vec3 db = glm::normalize(b - origin);
  return glm::acos(glm::dot(da, db));
}

void SandboxLayer::OnAttach()
{
	EnableGLDebugging();
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  //glLineWidth(1);
  //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  model_loader_ = std::make_unique<ModelLoader>();// std::make_unique<ModelLoader>("src/assets/backpack.obj", false);
  // Scene should also be returned to access the the right node, maybe storing each mesh into a hash map.

  line_ = std::make_unique<LineTools>(this->camera_controller_.GetCamera().GetPosition());
  mesh_ = std::move(model_loader_->LoadModel("src/assets/monkey.obj"));
  plane_ = std::move(model_loader_->LoadModel("src/assets/Plane.obj"));
  m_sim = std::make_unique<Simulation>();
  {
    // Test scene: a 5x5 patch of strands, starting horizontal so they swing
    // down under gravity and come to rest hanging.
    std::vector<EigenVector3> roots;
    for (int i = 0; i < 5; ++i)
      for (int j = 0; j < 5; ++j)
        roots.emplace_back(-0.4 + 0.2 * i, 1.5, -0.4 + 0.2 * j);
    m_sim->InitHairStrands(roots, EigenVector3(1.0, 0.0, 0.0), 20, 0.05);
  }
  m_strand_renderer = std::make_unique<StrandRenderer>();
  m_strand_renderer->Upload(m_sim->GetPositions(), m_sim->GetSegmentIndices());

  this->kScene = model_loader_->getScene();

  this->mesh_->TranslateShape(glm::vec3(0.0f, 0.0f, 0.0f));
  plane_->TranslateShape(glm::vec3(0.0f, -1.0f, 0.0f));
}

void SandboxLayer::OnDetach()
{
	// Shutdown here
}

int SandboxLayer::PickParticle(float mouse_x, float mouse_y, float* out_ndc_depth) const
{
  const glm::mat4& vp = camera_controller_.GetCamera().GetViewProjectionMatrix();
  const EigenMatrix3x& pos = m_sim->GetPositions();

  int best = -1;
  float best_dist2 = kPickRadiusPx * kPickRadiusPx;
  float best_depth = 0.0f;
  for (Eigen::Index i = 0; i < pos.rows(); ++i) {
    glm::vec4 clip = vp * glm::vec4(pos(i, 0), pos(i, 1), pos(i, 2), 1.0);
    if (clip.w <= 1e-6f) continue;  // Behind the camera.
    glm::vec3 ndc = glm::vec3(clip) / clip.w;
    // NDC -> window pixels (origin top-left, like the mouse coordinates).
    float sx = (ndc.x * 0.5f + 0.5f) * width;
    float sy = (1.0f - (ndc.y * 0.5f + 0.5f)) * height;
    float d2 = (sx - mouse_x) * (sx - mouse_x) + (sy - mouse_y) * (sy - mouse_y);
    // Prefer the closest on screen; break near-ties with the one nearer the camera.
    if (d2 < best_dist2 - 1.0f || (best >= 0 && std::abs(d2 - best_dist2) <= 1.0f && ndc.z < best_depth)) {
      best = static_cast<int>(i);
      best_dist2 = d2;
      best_depth = ndc.z;
    }
  }
  if (out_ndc_depth) *out_ndc_depth = best_depth;
  return best;
}

glm::vec3 SandboxLayer::MouseToWorld(float mouse_x, float mouse_y, float ndc_depth) const
{
  glm::vec4 ndc(2.0f * mouse_x / width - 1.0f, 1.0f - 2.0f * mouse_y / height, ndc_depth, 1.0f);
  glm::vec4 world = glm::inverse(camera_controller_.GetCamera().GetViewProjectionMatrix()) * ndc;
  return glm::vec3(world) / world.w;
}

void SandboxLayer::OnEvent(Event& event)
{
  camera_controller_.OnEvent(event);
  EventDispatcher dispatcher(event);
  dispatcher.Dispatch<WindowResizeEvent>(
    [&](WindowResizeEvent& e)
    {
      if (e.GetWidth() > 0 && e.GetHeight() > 0) {
        width = e.GetWidth();
        height = e.GetHeight();
      }
      return false;
    });
  // Left click picks the nearest hair particle and drags it while held.
  // (ImGui windows swallow the press first when the cursor is over them.)
  dispatcher.Dispatch<MouseButtonPressedEvent>(
    [&](MouseButtonPressedEvent& e)
    {
      if (e.GetMouseButton() != HZ_MOUSE_BUTTON_LEFT) return false;
      float depth = 0.0f;
      int particle = PickParticle(e.GetMouseX(), e.GetMouseY(), &depth);
      if (particle < 0) return false;
      m_drag_ndc_depth = depth;
      m_sim->BeginDrag(particle, GLM2Eigen(MouseToWorld(e.GetMouseX(), e.GetMouseY(), depth)));
      return true;
    });
  dispatcher.Dispatch<MouseButtonReleasedEvent>(
    [&](MouseButtonReleasedEvent& e)
    {
      if (e.GetMouseButton() != HZ_MOUSE_BUTTON_LEFT || !m_sim->IsDragging()) return false;
      m_sim->EndDrag();
      return true;
    });
}

void SandboxLayer::OnUpdate(Timestep ts)
{
  // Don't turn the camera while the cursor is over ImGui or a strand is being dragged.
  camera_controller_.SetMouseLookBlocked(ImGui::GetIO().WantCaptureMouse || m_sim->IsDragging());
  camera_controller_.OnUpdate(ts);
  // Follow the mouse every frame (also keeps the target right while the camera moves).
  if (m_sim->IsDragging()) {
    m_sim->UpdateDrag(GLM2Eigen(MouseToWorld(Input::GetMouseX(), Input::GetMouseY(), m_drag_ndc_depth)));
  }
	// Render here
  glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
  glStencilFunc(GL_ALWAYS, 1, 0xFF);
  glStencilMask(0xFF);
  glLineWidth(2);


  this->mesh_->DrawBoundingBox(this->camera_controller_);
  this->plane_->DrawBoundingBox(this->camera_controller_);
  this->mesh_->Draw(this->camera_controller_);
  this->plane_->Draw(this->camera_controller_);
  this->line_->DrawDebugLine(this->camera_controller_);

  if (m_run_sim) {
    m_sim->Update(ts.GetSeconds());
    m_strand_renderer->Upload(m_sim->GetPositions(), m_sim->GetSegmentIndices());
  }
  m_strand_renderer->Draw(camera_controller_);
}

void SandboxLayer::OnImGuiRender() {
  ImGui::Begin("Hair Simulation");
  if (ImGui::BeginTabBar("Tabs")) {
    if (ImGui::BeginTabItem("Simulation")) {
      DrawSimulationTab();
      ImGui::EndTabItem();
    }
    if (ImGui::BeginTabItem("Camera")) {
      DrawCameraTab();
      ImGui::EndTabItem();
    }
    ImGui::EndTabBar();
  }
  ImGui::End();
}

void SandboxLayer::DrawCameraTab() {
  float sensitivity = camera_controller_.GetMouseSensitivity();
  if (ImGui::SliderFloat("Mouse sensitivity", &sensitivity, 0.01f, 1.0f, "%.2f deg/px"))
    camera_controller_.SetMouseSensitivity(sensitivity);

  int mode = static_cast<int>(camera_controller_.GetMouseLookMode());
  ImGui::Text("Mouse look");
  bool changed = ImGui::RadioButton("Hold Left Alt", &mode, static_cast<int>(MouseLookMode::HoldAlt));
  ImGui::SameLine();
  changed |= ImGui::RadioButton("Always (no key)", &mode, static_cast<int>(MouseLookMode::Always));
  if (changed) camera_controller_.SetMouseLookMode(static_cast<MouseLookMode>(mode));

  bool invert = camera_controller_.GetInvertY();
  if (ImGui::Checkbox("Invert Y", &invert)) camera_controller_.SetInvertY(invert);

  ImGui::TextDisabled("WASD to move. Look pauses over this window and while dragging hair.");
}

void SandboxLayer::DrawSimulationTab() {
  ImGui::Checkbox("Run", &m_run_sim);
  if (ImGui::Button("Step") && !m_run_sim) {
    m_sim->Step();
    m_strand_renderer->Upload(m_sim->GetPositions(), m_sim->GetSegmentIndices());
  }
  ImGui::SameLine();
  if (ImGui::Button("Reset")) {
    m_sim->Reset();
    m_strand_renderer->Upload(m_sim->GetPositions(), m_sim->GetSegmentIndices());
  }

  ImGui::SliderInt("Iterations", &m_sim->m_iterations_per_step, 1, 50);
  float damping = static_cast<float>(m_sim->m_damping_coefficient);
  if (ImGui::SliderFloat("Damping", &damping, 0.0f, 0.2f)) m_sim->m_damping_coefficient = damping;
  float gravity = static_cast<float>(m_sim->m_gravity_constant);
  if (ImGui::SliderFloat("Gravity", &gravity, 0.0f, 20.0f)) m_sim->m_gravity_constant = gravity;

  // Stiffness changes the system matrix, so it has to be refactored.
  bool rebuild = false;
  float stretch = static_cast<float>(m_sim->m_stiffness_stretch);
  if (ImGui::SliderFloat("Stretch stiffness", &stretch, 1e3f, 1e7f, "%.0f", 4.0f)) {
    m_sim->m_stiffness_stretch = stretch;
    rebuild = true;
  }
  float bending = static_cast<float>(m_sim->m_stiffness_bending);
  if (ImGui::SliderFloat("Bending stiffness", &bending, 1.0f, 1e5f, "%.0f", 4.0f)) {
    m_sim->m_stiffness_bending = bending;
    rebuild = true;
  }
  float drag = static_cast<float>(m_sim->m_stiffness_drag);
  if (ImGui::SliderFloat("Mouse drag stiffness", &drag, 1e2f, 1e6f, "%.0f", 4.0f)) {
    m_sim->m_stiffness_drag = drag;
    rebuild = true;
  }
  if (rebuild) m_sim->RebuildSystemMatrix();
  ImGui::Text(m_sim->IsDragging() ? "Dragging particle %d" : "Left-click a strand to drag it", m_sim->DraggedParticle());
}