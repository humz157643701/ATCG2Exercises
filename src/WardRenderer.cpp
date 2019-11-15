#include "WardRenderer.h"
#include <libheaders.h>
#include <Scene.h>
#include <Shader.h>

WardRenderer::WardRenderer() :
	m_opaque_shader(ShaderProgram::createShaderProgram("assets/shaders/ward/ward_opaque.vert", "assets/shaders/ward/ward_opaque.frag"))
{
	
}

WardRenderer::~WardRenderer()
{
}

void WardRenderer::render(Scene * scene, double dt, bool measure)
{
	glClearColor(scene->clearColor.x, scene->clearColor.y, scene->clearColor.z, scene->clearColor.w);
	glClearDepth(1.0f);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glDepthMask(GL_TRUE);
	glDisable(GL_BLEND);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	m_opaque_shader.use();
	scene->m_camera.bind(&m_opaque_shader, "camera");

	for (size_t i = 0; i < scene->m_dirlights.size(); ++i)
		scene->m_dirlights[i].bind(&m_opaque_shader, scene->m_camera.getViewMatrix(), ("dirlights[" + std::to_string(i) + "]").c_str());
	m_opaque_shader.setUniform("dirlightcount", static_cast<GLint>(scene->m_dirlights.size()));

	for (size_t i = 0; i < scene->m_pointlights.size(); ++i)
		scene->m_pointlights[i].bind(&m_opaque_shader, scene->m_camera.getViewMatrix(), ("pointlights[" + std::to_string(i) + "]").c_str());
	m_opaque_shader.setUniform("pointlightcount", static_cast<GLint>(scene->m_pointlights.size()));

	for (size_t i = 0; i < scene->m_ambientlights.size(); ++i)
		scene->m_ambientlights[i].bind(&m_opaque_shader, ("ambientlights[" + std::to_string(i) + "]").c_str());
	m_opaque_shader.setUniform("ambientlightcount", static_cast<GLint>(scene->m_ambientlights.size()));

	m_opaque_shader.setUniform("exposure", scene->tm_exposure);

	scene->drawOpaqueWithMaterials(&m_opaque_shader);
}

double WardRenderer::getLastTransparentRenderTime()
{
	return 0.0;
}
