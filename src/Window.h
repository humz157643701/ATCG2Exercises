/** \addtogroup application
Application logic
*  @{
*/

/*!
\file Window.h
*/

#include <GameWindow.h>
#include <memory>
#include <Scene.h>
#include <string>
#include <sstream>
#include <BSpline.h>
#include <PerfDisplay.h>


//! Main application class
class Window : public GameWindow
{
public:
	Window();
	~Window();

private:
	//! Do simulation tasks such as triggering simulations steps for the particle systems or move the camera
	GLvoid update(GLdouble dtime) override;
	//! Render the current scene with the currently selected renderer
	GLvoid render(GLdouble dtime) override;
	//! Initialize all resources, load the scene list...
	GLvoid init() override;
	//! Clean up the mess
	GLvoid shutdown() override;

	//! Callback for window resize events
	void onWindowResize(int width, int height) override;
	//! Callback for framebuffer resize events
	void onFrameBufferResize(int width, int height) override;

	//! Callback for key events
	void onKey(Key key, Action action, Modifier modifier) override;
	//! Callback for mouse move events
	void onMouseMove(MousePosition mouseposition) override;
	//! Callback for mouse button events
	void onMouseButton(MouseButton button, Action action, Modifier modifier) override;
	//! Callback for mouse scroll events
	void onMouseScroll(double xscroll, double yscroll) override;

	//! Loads the list of scenes from the scene file
	void loadSceneList();
	//! Loads the scene with index index
	void loadScene(size_t index);

	//data
	//! Current scene
	std::unique_ptr<Scene> m_scene;
	//! Current renderer
	std::unique_ptr<IRenderer> m_renderer;
	//! scene simulation stops when this flag is set false
	bool m_updateScene;

	//scene list
	//! Paths to all available scene files
	std::vector<std::string> m_scenes;
	//! Index of the current scene
	size_t m_currentScene;
	//! BSpline for camera travel
	BSpline m_bspline;
	glm::vec3 m_camera_lookat;
	bool m_cam_spline_move;
	float m_current_spline_pos;
	float m_spline_speed;
	//! Lock the camera view
	bool lockview;

	// For FPS display
	TextRenderer textrenderer;
	bool displayFPS;
	float fpsAccum;
	float fpspos;
};

/** @}*/