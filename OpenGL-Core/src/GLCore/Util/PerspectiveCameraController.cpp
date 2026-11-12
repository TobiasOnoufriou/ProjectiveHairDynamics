#include "glpch.h"
#include "PerspectiveCameraController.h"

#include "GLCore/Core/Input.h"
#include "GLCore/Core/KeyCodes.h"

#include <glm/gtc/constants.hpp>

namespace GLCore::Utils {

	PerspectiveCameraController::PerspectiveCameraController(float aspectRatio, bool rotation)
		: m_AspectRatio(aspectRatio), m_Camera(kFovDegrees, m_AspectRatio, kNearClip, kFarClip), m_Rotation(rotation)
	{
    UpdateOrientation();
	}

	void PerspectiveCameraController::UpdateOrientation()
	{
    m_Camera.SetDirection(glm::vec3{
      cos(vertical_angle_) * sin(horizontal_angle_),
      sin(vertical_angle_),
      cos(vertical_angle_) * cos(horizontal_angle_) }
    );

    m_Camera.SetRight(glm::vec3{
     sin(horizontal_angle_ - glm::half_pi<float>()),
     0,
     cos(horizontal_angle_ - glm::half_pi<float>())
    });

    m_Camera.SetUp(glm::cross(m_Camera.GetRight(), m_Camera.GetDirection()));
	}

	void PerspectiveCameraController::OnUpdate(Timestep ts)
	{
		if (Input::IsKeyPressed(HZ_KEY_A))
		{
      m_CameraPosition -= m_Camera.GetRight() * m_CameraTranslationSpeed * float(ts);
		}
		else if (Input::IsKeyPressed(HZ_KEY_D))
		{
      m_CameraPosition += m_Camera.GetRight() * m_CameraTranslationSpeed * float(ts);
		}

		if (Input::IsKeyPressed(HZ_KEY_W))
		{
      m_CameraPosition += m_Camera.GetDirection() * m_CameraTranslationSpeed * float(ts);
		}
		else if (Input::IsKeyPressed(HZ_KEY_S))
		{
      m_CameraPosition -= m_Camera.GetDirection() * m_CameraTranslationSpeed * float(ts);
		}

    // Mouse look: rotate by how far the cursor moved since last frame.
    {
      const bool wants_look = !m_MouseLookBlocked &&
        (m_MouseLookMode == MouseLookMode::Always || Input::IsKeyPressed(HZ_KEY_LEFT_ALT));
      const float mx = Input::GetMouseX(), my = Input::GetMouseY();
      if (wants_look) {
        if (m_MouseLookActive) {
          const float rad_per_px = glm::radians(m_MouseSensitivity);
          const float dy = (my - m_LastMouseY) * (m_InvertY ? -1.0f : 1.0f);
          horizontal_angle_ -= rad_per_px * (mx - m_LastMouseX);
          vertical_angle_ -= rad_per_px * dy;
          // Stop just short of straight up/down so the view never flips.
          const float limit = glm::radians(89.0f);
          vertical_angle_ = glm::clamp(vertical_angle_, -limit, limit);
          UpdateOrientation();
        }
      }
      m_MouseLookActive = wants_look;
      m_LastMouseX = mx;
      m_LastMouseY = my;
    }

		if (m_Rotation)
		{
			if (Input::IsKeyPressed(HZ_KEY_Q))
				m_CameraRotation += m_CameraRotationSpeed * ts;
			if (Input::IsKeyPressed(HZ_KEY_E))
				m_CameraRotation -= m_CameraRotationSpeed * ts;

			if (m_CameraRotation > 180.0f)
				m_CameraRotation -= 360.0f;
			else if (m_CameraRotation <= -180.0f) 
				m_CameraRotation += 360.0f;

			m_Camera.SetRotation(m_CameraRotation);
		}

    m_Camera.SetPosition(m_CameraPosition);
		m_CameraTranslationSpeed = m_ZoomLevel;
	}

	void PerspectiveCameraController::OnEvent(Event& e)
	{
		EventDispatcher dispatcher(e);
		dispatcher.Dispatch<MouseScrolledEvent>(GLCORE_BIND_EVENT_FN(PerspectiveCameraController::OnMouseScrolled));
		dispatcher.Dispatch<WindowResizeEvent>(GLCORE_BIND_EVENT_FN(PerspectiveCameraController::OnWindowResized));
	}

	bool PerspectiveCameraController::OnMouseScrolled(MouseScrolledEvent& e)
	{
    return false;
	}

	void PerspectiveCameraController::SetViewportSize(uint32_t width, uint32_t height)
	{
		if (width == 0 || height == 0) return;  // Minimised.
		m_ViewportWidth = width;
		m_ViewportHeight = height;
		m_AspectRatio = (float)width / (float)height;
		m_Camera.SetProjection(kFovDegrees, m_AspectRatio, kNearClip, kFarClip);
	}

	bool PerspectiveCameraController::OnWindowResized(WindowResizeEvent& e)
  {
    SetViewportSize(e.GetWidth(), e.GetHeight());
    return false;
	}

}