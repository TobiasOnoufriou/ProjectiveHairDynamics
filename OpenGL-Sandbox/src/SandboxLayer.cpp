#include "SandboxLayer.hpp"
#include "Mesh.hpp"
#include "glm/gtx/rotate_vector.hpp"
#include "GLCore/Core/Input.h"
#include "GLCore/Core/MouseButtonCodes.h"
#include "HairScene.h"
#include "MuscleScene.h"
#include "TetScene.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <memory>

using namespace GLCore;
using namespace GLCore::Utils;

namespace {

// Selectable scenes.
struct SceneEntry {
  const char* name;
  std::unique_ptr<Scene> (*create)();
};
const SceneEntry kScenes[] = {
  {"Hair", []() -> std::unique_ptr<Scene> { return std::make_unique<HairScene>(); }},
  {"Tet beam", []() -> std::unique_ptr<Scene> { return std::make_unique<TetScene>(); }},
  {"Muscle", []() -> std::unique_ptr<Scene> { return std::make_unique<MuscleScene>(); }},
};
constexpr int kNumScenes = static_cast<int>(sizeof(kScenes) / sizeof(kScenes[0]));

}  // namespace

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
  m_sim = std::make_unique<PDSolver>();
  m_sim->m_track_energy = true;
  // Anderson needs ~4x fewer iterations than plain PD for the same error on the
  // tet scenes at no measurable per-iteration cost (see the ConvergenceVsTime
  // benchmark), so the app starts with it; the solver default stays plain.
  m_sim->m_acceleration = PDAcceleration::kAnderson;
  m_strand_renderer = std::make_unique<StrandRenderer>();
  m_surface_renderer = std::make_unique<SurfaceRenderer>();
  m_collider_renderer = std::make_unique<SurfaceRenderer>();
  m_collider_renderer->SetColor(glm::vec3(0.45f, 0.50f, 0.58f));
  LoadScene(0);

  this->kScene = model_loader_->getScene();

  this->mesh_->TranslateShape(glm::vec3(0.0f, 0.0f, 0.0f));
  plane_->TranslateShape(glm::vec3(0.0f, -1.0f, 0.0f));
}

void SandboxLayer::LoadScene(int index)
{
  if (index < 0 || index >= kNumScenes) return;
  m_sim->EndDrag();
  m_scene_index = index;
  m_scene = kScenes[index].create();
  m_scene->Build(m_sim.get());
  UploadStrands();
}

void SandboxLayer::UploadStrands()
{
  m_strand_renderer->Upload(m_sim->GetPositions(), m_scene->LineIndices());
  m_surface_renderer->Upload(m_sim->GetPositions(), m_scene->SurfaceIndices());
  UploadColliders();
}

void SandboxLayer::UploadColliders()
{
  std::vector<EigenVector3> points;
  std::vector<unsigned int> indices;
  if (m_show_colliders) {
    for (const Collider& c : m_sim->GetColliders()) {
      if (c.visible) TriangulateCollider(c, 4.0, &points, &indices);
    }
  }
  EigenMatrix3x positions(points.size(), 3);
  for (size_t k = 0; k < points.size(); ++k) positions.row(k) = points[k].transpose();
  m_collider_renderer->Upload(positions, indices);
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
    if (m_sim->IsFixed(static_cast<unsigned int>(i))) continue;  // Pinned/bone vertices can't be dragged.
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
    m_scene->Update(m_sim.get(), ts.GetSeconds());  // Kinematic targets (bones) first.
    m_sim->Update(ts.GetSeconds());
    UploadStrands();
  }
  // Exponential moving average so the ms/step readout is legible.
  m_avg_step_ms += 0.05 * (m_sim->GetStats().step_ms - m_avg_step_ms);
  m_surface_renderer->Draw(camera_controller_);
  m_collider_renderer->Draw(camera_controller_);
  m_strand_renderer->Draw(camera_controller_);
}

void SandboxLayer::OnImGuiRender() {
  ImGui::Begin("Projective Dynamics");
  if (ImGui::BeginTabBar("Tabs")) {
    if (ImGui::BeginTabItem("Simulation")) {
      DrawSimulationTab();
      ImGui::EndTabItem();
    }
    if (ImGui::BeginTabItem("Metrics")) {
      DrawMetricsTab();
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
  if (ImGui::BeginCombo("Scene", kScenes[m_scene_index].name)) {
    for (int i = 0; i < kNumScenes; ++i) {
      if (ImGui::Selectable(kScenes[i].name, i == m_scene_index)) LoadScene(i);
    }
    ImGui::EndCombo();
  }

#ifdef _DEBUG
  ImGui::TextColored(ImVec4(1, 0.6f, 0.2f, 1), "Debug build: the solver runs 10-100x slower. Use Release.");
#endif
  ImGui::Checkbox("Run", &m_run_sim);
  if (ImGui::Button("Step") && !m_run_sim) {
    m_scene->Update(m_sim.get(), m_sim->m_h);
    m_sim->Step();
    UploadStrands();
  }
  ImGui::SameLine();
  if (ImGui::Button("Reset")) {
    m_sim->Reset();
    UploadStrands();
  }
  ImGui::SameLine();
  if (ImGui::Checkbox("Show colliders", &m_show_colliders)) UploadColliders();

  ImGui::SliderInt("Iterations", &m_sim->m_iterations_per_step, 1, 50);
  // Iterate extrapolation (RHS-side only, the system matrix is unchanged).
  static const char* kAccelerations[] = {"None", "Chebyshev (Wang 2015)", "Anderson (Peng et al. 2018)"};
  int acceleration = static_cast<int>(m_sim->m_acceleration);
  if (ImGui::Combo("Acceleration", &acceleration, kAccelerations, IM_ARRAYSIZE(kAccelerations))) {
    m_sim->m_acceleration = static_cast<PDAcceleration>(acceleration);
  }
  if (m_sim->m_acceleration == PDAcceleration::kChebyshev) {
    float rho = static_cast<float>(m_sim->m_chebyshev_rho);
    if (ImGui::SliderFloat("Spectral radius rho", &rho, 0.5f, 0.9999f, "%.4f")) m_sim->m_chebyshev_rho = rho;
    float gamma = static_cast<float>(m_sim->m_chebyshev_gamma);
    if (ImGui::SliderFloat("Under-relaxation gamma", &gamma, 0.3f, 1.0f)) m_sim->m_chebyshev_gamma = gamma;
  } else if (m_sim->m_acceleration == PDAcceleration::kAnderson) {
    ImGui::SliderInt("Anderson window", &m_sim->m_anderson_window, 1, PDSolver::kMaxAndersonWindow);
  }
  float damping = static_cast<float>(m_sim->m_damping_coefficient);
  if (ImGui::SliderFloat("Damping", &damping, 0.0f, 0.2f)) m_sim->m_damping_coefficient = damping;
  float gravity = static_cast<float>(m_sim->m_gravity_constant);
  if (ImGui::SliderFloat("Gravity", &gravity, 0.0f, 20.0f)) m_sim->m_gravity_constant = gravity;

  // Stiffness changes the system matrix, so it has to be refactored.
  bool rebuild_matrix = false;
  bool rebuild_scene = false;
  float drag = static_cast<float>(m_sim->m_stiffness_drag);
  if (ImGui::SliderFloat("Mouse drag stiffness", &drag, 1e2f, 1e6f, "%.0f", 4.0f)) {
    m_sim->m_stiffness_drag = drag;
    rebuild_matrix = true;
  }

  ImGui::Separator();
  ImGui::Text("%s parameters", m_scene->Name());
  for (const SceneParam& p : m_scene->Params()) {
    float v = static_cast<float>(*p.value);
    const bool changed = ImGui::SliderFloat(p.name, &v, p.min, p.max, p.logarithmic ? "%.3g" : "%.3f",
                                            p.logarithmic ? 4.0f : 1.0f);
    if (changed) *p.value = v;
    if (changed && p.effect == SceneParam::kSystemMatrix) rebuild_matrix = true;
    // Rebuilding the scene resets the pose, so wait until the slider is released.
    if (p.effect == SceneParam::kRebuildScene && ImGui::IsItemDeactivatedAfterEdit()) rebuild_scene = true;
  }
  if (auto* tet = dynamic_cast<TetScene*>(m_scene.get())) {
    if (ImGui::Checkbox("Show tet edges", &tet->m_show_wireframe)) UploadStrands();
    ImGui::Text("%zu tets, %zu surface triangles, volume %.4f m^3", tet->Mesh().NumTets(),
                tet->SurfaceIndices().size() / 3, tet->Mesh().TotalVolume());
  }
  if (auto* muscle = dynamic_cast<MuscleScene*>(m_scene.get())) {
    ImGui::Checkbox("Animate elbow + activation", &muscle->m_animate);
    ImGui::Checkbox("Muscle-bone contact", &muscle->m_bone_contact);
    ImGui::SameLine();
    // Which vertices have contact slots is part of the system: rebuild.
    if (ImGui::Checkbox("Interior vertices too", &muscle->m_interior_contact)) rebuild_scene = true;
    if (ImGui::Checkbox("Show tet edges", &muscle->m_show_wireframe)) UploadStrands();
    // Current muscle volume over rest volume.
    const TetMesh& mesh = muscle->Mesh();
    ScalarType volume = 0.0;
    for (size_t t = 0; t < mesh.NumTets(); ++t) volume += mesh.SignedVolume(t, m_sim->GetPositions());
    ImGui::Text("%zu tets, muscle volume %.1f%% of rest", mesh.NumTets(), 100.0 * volume / mesh.TotalVolume());
  }
  if (rebuild_scene) {
    m_sim->EndDrag();
    m_scene->Build(m_sim.get());
    UploadStrands();
  } else if (rebuild_matrix) {
    m_sim->RebuildSystemMatrix();
  }
  ImGui::Text(m_sim->IsDragging() ? "Dragging vertex %d" : "Left-click a vertex to drag it", m_sim->DraggedVertex());
}

void SandboxLayer::DrawMetricsTab() {
  const PDStepStats& stats = m_sim->GetStats();
  ImGui::Text("Vertices: %u (%u fixed)   Constraints: %zu", m_sim->NumVertices(), m_sim->NumFixed(),
              m_sim->NumConstraints());
  ImGui::Text("Step: %.3f ms (avg %.3f ms)", stats.step_ms, m_avg_step_ms);
  ImGui::Text("  local %.3f ms, global %.3f ms", stats.local_ms, stats.global_ms);
  ImGui::Text("Max strain: %.3e   Max volume error: %.3e", stats.max_strain, stats.max_volume_error);
  if (m_sim->m_acceleration != PDAcceleration::kNone) {
    ImGui::Text("Rejected accelerated iterates: %d", stats.rejected_iterations);
  }
  if (!m_sim->GetColliders().empty()) {
    ImGui::Text("Contacts: %d of %zu slots, max penetration %.2f mm", stats.num_contacts, m_sim->NumContactSlots(),
                1e3 * stats.max_penetration);
  }

  // Race-free local step: every constraint projects into its own slot, so the
  // OpenMP loop gives the same result as the serial one. Compare the ms/step.
  ImGui::Checkbox("Parallel local step (OpenMP)", &m_sim->m_parallel_local_step);
  ImGui::SameLine();
  ImGui::TextDisabled("%d threads", PDSolver::MaxThreads());

  ImGui::Checkbox("Track energy per iteration", &m_sim->m_track_energy);
  if (m_sim->m_track_energy && !stats.energy.empty()) {
    // Plot E(x_k) - E_final on a log scale so convergence is visible across
    // many orders of magnitude.
    const ScalarType e_final = stats.energy.back();
    m_energy_plot.clear();
    for (ScalarType e : stats.energy) {
      m_energy_plot.push_back(static_cast<float>(std::log10(std::max(e - e_final, 1e-12))));
    }
    ImGui::Text("Energy: %.6e -> %.6e", stats.energy.front(), e_final);
    ImGui::PlotLines("log10(E_k - E_final)", m_energy_plot.data(), static_cast<int>(m_energy_plot.size()),
                     0, nullptr, -12.0f, FLT_MAX, ImVec2(0, 80));
    bool monotone = true;
    for (size_t k = 1; k < stats.energy.size(); ++k) {
      monotone &= stats.energy[k] <= stats.energy[k - 1] * (1.0 + 1e-9) + 1e-12;
    }
    if (!monotone) ImGui::TextColored(ImVec4(1, 0.4f, 0.4f, 1), "Energy increased during the last step");
  }
}
