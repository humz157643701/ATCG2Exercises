#include "Window.h"
#include <fstream>
#include <SceneLoader.h>
#include <WardRenderer.h>

Window::Window() :
	GameWindow(
		1280,
		720,
		false,
		true,
		4,
		5,
		"Ward Renderer",
		4
	),
	textrenderer("assets/fonts/arial.ttf", 128, 500)
{
}

Window::~Window()
{
}

GLvoid Window::update(GLdouble dtime)
{
	if (m_updateScene)
	{
		m_scene->update(static_cast<float>(dtime));
	}

	if (displayFPS)
	{
		fpsAccum += dtime;
		if (fpsAccum >= static_cast<float>(PERF_INTERVAL))
		{
			textrenderer.setString("Avg. FPS: " + std::to_string(m_currentTM.avgfps) + "\nMin. FPS: " + std::to_string(m_currentTM.minfps),
				{ fpspos, 0.1f }, 0.05f, getAspectRatio(), { 1.0f, 1.0f, 1.0f, 1.0f }, 0.2f);
			fpsAccum = 0.0f;
		}

		if (input->getKeyState(Key::PageUp) == KeyState::Pressed)
			fpspos += 0.2f * static_cast<float>(dtime);
		if (input->getKeyState(Key::PageDown) == KeyState::Pressed)
			fpspos -= 0.2f * static_cast<float>(dtime);
	}

	if (!m_cam_spline_move)
	{
		float mod = (input->getKeyState(Key::LeftShift) == KeyState::Pressed ? 4.0f : 1.0f);
		if (input->getKeyState(Key::W) == KeyState::Pressed)
		{
			m_scene->m_camera.forward(dtime * 5.0f * mod);
		}
		if (input->getKeyState(Key::S) == KeyState::Pressed)
		{
			m_scene->m_camera.backward(dtime * 5.0f * mod);
		}
		if (input->getKeyState(Key::A) == KeyState::Pressed)
		{
			m_scene->m_camera.left(dtime * 5.0f * mod);
		}
		if (input->getKeyState(Key::D) == KeyState::Pressed)
		{
			m_scene->m_camera.right(dtime * 5.0f * mod);
		}
	}
	else
	{
		m_scene->m_camera.transform.setWorldPosition(m_bspline.at(m_current_spline_pos));
		m_scene->m_camera.transform.lookinto(glm::normalize(m_camera_lookat - m_scene->m_camera.transform.getWorldPosition()));
		m_scene->m_camera.setViewDirty();
		m_current_spline_pos += dtime * m_spline_speed;
		if (m_current_spline_pos > 1.0f)
			m_current_spline_pos = 0.0f;

		if (input->getKeyState(Key::W) == KeyState::Pressed)
		{
			m_spline_speed += dtime * 0.1f;
		}
		if (input->getKeyState(Key::S) == KeyState::Pressed)
		{
			m_spline_speed = glm::max(m_spline_speed - static_cast<float>(dtime) * 0.1f, 0.0f);
		}
	}
}

GLvoid Window::render(GLdouble dtime)
{
	GLint newvp[] = { GLint(getFrameBufferWidth() / 2), getFrameBufferHeight() };
	m_scene->m_camera.setViewport(newvp[0], newvp[1]);

	for (auto x : m_scene->m_opaque_models.front().meshes)
	{
		x->setMaterial(m_scene->m_materials.front().get());
	}

	glViewport(0, 0, newvp[0], newvp[1]);
	m_renderer->render(m_scene.get(), dtime);

	//Set material of slate box to NN stuff
	for (auto x : m_scene->m_opaque_models.front().meshes)
	{
		x->setMaterial(m_scene->m_materials.back().get());
	}
	glViewport(newvp[0], 0, newvp[0], getFrameBufferHeight());
	//Render once via ward renderer foor the things we don't have a NN material of
	m_renderer->render(m_scene.get(), dtime, false, false);
	//Temp copy of our models. Should copy fine, I hope?
	auto tempmodels = std::vector<Model>(m_scene->m_opaque_models.begin()+1, m_scene->m_opaque_models.end());
	//Remove non axf material models
	m_scene->m_opaque_models.erase(m_scene->m_opaque_models.begin() + 1, m_scene->m_opaque_models.end());
	// TODO Render using ggx renderer


	//Put models back in
	m_scene->m_opaque_models.insert(m_scene->m_opaque_models.end(), tempmodels.begin(), tempmodels.end());

	if (displayFPS)
	{
		textrenderer.draw();
	}
}

GLvoid Window::init() 
{
	setCursorVisible(false);
	loadSceneList();
	if (m_scenes.size() == 0)
		throw std::logic_error("Scene list is empty");
	loadScene(0);
	m_currentScene = 0;

	m_cam_spline_move = false;
	m_current_spline_pos = 0.0f;
	m_spline_speed = 0.1f;

	// create renderer
	m_renderer = std::unique_ptr<WardRenderer>(new WardRenderer(m_scene.get()));

	m_updateScene = true;

	lockview = false;

	displayFPS = false;
	fpspos = 0.625f;
	fpsAccum = 0.0f;
	textrenderer.setString("", { fpspos, 0.06 }, 0.05, getAspectRatio(), { 1.0f, 1.0f, 1.0f, 1.0f }, 0.2f);
}

GLvoid Window::shutdown()
{
	// TODO: cleanup
}

void Window::onWindowResize(int width, int height)
{
}

void Window::onFrameBufferResize(int width, int height)
{
	m_scene->m_camera.setViewport(width, height);
	// recreate renderer
	m_renderer.reset(new WardRenderer(m_scene.get()));
}

void Window::onKey(Key key, Action action, Modifier modifier)
{
	if (key == Key::Escape && action == Action::Down)
		quit();

	if (key == Key::U && action == Action::Down)
	{
		m_updateScene = !m_updateScene;
	}

	if (key == Key::L && action == Action::Down)
	{
		lockview = !lockview;
	}

	// set lookat point for camera spline
	if (key == Key::H && action == Action::Down && !m_cam_spline_move)
	{
		m_camera_lookat = m_scene->m_camera.transform.getWorldPosition();
	}
	
	// add control point to camera spline
	if (key == Key::B && action == Action::Down && !m_cam_spline_move)
	{
		m_bspline.appendControlPoint(m_scene->m_camera.transform.getWorldPosition());
	}

	// clear camera spline
	if (key == Key::N && action == Action::Down && !m_cam_spline_move)
	{
		m_bspline.clear();
	}

	// start camera travel along spline
	if (key == Key::M && action == Action::Down)
	{
		m_cam_spline_move = !m_cam_spline_move;
		m_current_spline_pos = 0.0f;
	}

	if (key == Key::X && action == Action::Down)
		displayFPS = !displayFPS;

}

void Window::onMouseMove(MousePosition mouseposition)
{
	
	float dx = mouseposition.X - mouseposition.oldX;
	float dy = mouseposition.Y - mouseposition.oldY;
	if (!lockview && !m_cam_spline_move)
		m_scene->m_camera.rotateView(-dy * 0.005f, -dx * 0.005f);
}

void Window::onMouseButton(MouseButton button, Action action, Modifier modifier)
{
}

void Window::onMouseScroll(double xscroll, double yscroll)
{
}

void Window::loadSceneList()
{
	std::ifstream sl("assets/scenes/scenelist.txt");
	if (!sl.is_open())
	{
		throw std::logic_error("Scene list couldn't be loaded");
	}
	while (!sl.eof())
	{
		std::string sentry;
		std::getline(sl, sentry);
		if (sentry.length() > 0)
		{
			m_scenes.push_back("assets/scenes/" + sentry);
		}
	}
	m_currentScene = 0;
}

void Window::loadScene(size_t index)
{
	glFinish();
	m_scene.reset();
	m_scene = std::move(SceneLoader::loadScene(m_scenes[index]));
	m_scene->m_camera.setViewport(getFrameBufferWidth(), getFrameBufferHeight());
}