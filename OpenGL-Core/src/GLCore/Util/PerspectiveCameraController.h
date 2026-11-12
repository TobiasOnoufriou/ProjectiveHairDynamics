#pragma once

#include "PerspectiveCamera.h"
#include "GLCore/Core/Timestep.h"

#include "GLCore/Core/Constants.h"

#include "GLCore/Events/ApplicationEvent.h"
#include "GLCore/Events/MouseEvent.h"

namespace GLCore::Utils {

	enum class MouseLookMode
	{
		HoldAlt = 0,  // Rotate only while Left Alt is held.
		Always = 1,   // Rotate whenever the mouse moves.
	};

	class PerspectiveCameraController
	{
	public:
		PerspectiveCameraController(float aspectRatio, bool rotation = false);

		void OnUpdate(Timestep ts);
		void OnEvent(Event& e);

		PerspecitiveCamera& GetCamera() { return m_Camera; }
		const PerspecitiveCamera& GetCamera() const { return m_Camera; }

		float GetZoomLevel() const { return m_ZoomLevel; }
		void SetZoomLevel(float level) { m_ZoomLevel = level; }

		// Window size in pixels: sets the aspect ratio.
		void SetViewportSize(uint32_t width, uint32_t height);

		// Mouse look. Sensitivity is in degrees of rotation per pixel of mouse movement.
		float GetMouseSensitivity() const { return m_MouseSensitivity; }
		void SetMouseSensitivity(float degrees_per_pixel) { m_MouseSensitivity = degrees_per_pixel; }
		MouseLookMode GetMouseLookMode() const { return m_MouseLookMode; }
		void SetMouseLookMode(MouseLookMode mode) { m_MouseLookMode = mode; }
		bool GetInvertY() const { return m_InvertY; }
		void SetInvertY(bool invert) { m_InvertY = invert; }
		// Lets the owner pause mouse look for a frame (e.g. cursor over UI, dragging an object).
		void SetMouseLookBlocked(bool blocked) { m_MouseLookBlocked = blocked; }

    glm::vec3 GetRaycastDirection(int mouse_x, int mouse_y);
	private:
		void UpdateOrientation();
		bool OnMouseScrolled(MouseScrolledEvent& e);
		bool OnWindowResized(WindowResizeEvent& e);
	private:
		float m_AspectRatio;
		uint32_t m_ViewportWidth = 1280, m_ViewportHeight = 720;
		static constexpr float kNearClip = 0.01f;
		static constexpr float kFarClip = 1000.0f;
		static constexpr float kFovDegrees = 60.0f;
		float m_ZoomLevel = 1.0f;
    float m_FovY = 90.0f;
		PerspecitiveCamera m_Camera;

		bool m_Rotation;

		glm::vec3 m_CameraPosition = { 0.0f, 0.0f, 0.0f };
    glm::vec3 direction = {};
    float horizontal_angle_ = 0.0f, vertical_angle_ = 0.0f;
    float m_MouseSensitivity = 0.15f;  // Degrees per pixel.
    MouseLookMode m_MouseLookMode = MouseLookMode::HoldAlt;
    bool m_InvertY = false;
    bool m_MouseLookBlocked = false;
    bool m_MouseLookActive = false;    // Was looking last frame (false => re-seed last position).
    float m_LastMouseX = 0.0f, m_LastMouseY = 0.0f;
		float m_CameraRotation = 0.0f; //In degrees, in the anti-clockwise direction
		float m_CameraTranslationSpeed = 100.0f, m_CameraRotationSpeed = 180.0f;
	};

}