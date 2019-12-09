#include "OpaqueRenderer.h"
#include <libheaders.h>
#include <Scene.h>
#include <Shader.h>

OpaqueRenderer::OpaqueRenderer() :
	shader(ShaderProgram::createShaderProgram("assets/shaders/opaque/opaque.vert", "assets/shaders/opaque/opaque.frag"))
{
	
}

OpaqueRenderer::~OpaqueRenderer()
{
}

void OpaqueRenderer::render(Scene * scene, double dt, bool measure)
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

	shader.use();
	scene->m_camera.bind(&shader, "camera");

	for (size_t i = 0; i < scene->m_dirlights.size(); ++i)
		scene->m_dirlights[i].bind(&shader, scene->m_camera.getViewMatrix(), ("dirlights[" + std::to_string(i) + "]").c_str());
	shader.setUniform("dirlightcount", static_cast<GLint>(scene->m_dirlights.size()));

	for (size_t i = 0; i < scene->m_pointlights.size(); ++i)
		scene->m_pointlights[i].bind(&shader, scene->m_camera.getViewMatrix(), ("pointlights[" + std::to_string(i) + "]").c_str());
	shader.setUniform("pointlightcount", static_cast<GLint>(scene->m_pointlights.size()));

	for (size_t i = 0; i < scene->m_ambientlights.size(); ++i)
		scene->m_ambientlights[i].bind(&shader, ("ambientlights[" + std::to_string(i) + "]").c_str());
	shader.setUniform("ambientlightcount", static_cast<GLint>(scene->m_ambientlights.size()));

	shader.setUniform("tmwhite", scene->maxlum);
	//shader.setUniform("normalsfromtex", scene->enableNormalMapping);

	scene->drawOpaqueWithMaterials(&shader);
	scene->drawTransparentWithMaterials(&shader);
}

double OpaqueRenderer::getLastTransparentRenderTime()
{
	return 0.0;
}
