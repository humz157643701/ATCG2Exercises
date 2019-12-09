#include "WeightedBlendedOITRenderer.h"
#include <Scene.h>
#include <Primitives.h>
#include <chrono>

WeightedBlendedOITRenderer::WeightedBlendedOITRenderer(GLsizei width, GLsizei height, GLsizei samples) :
	opaqueshader(ShaderProgram::createShaderProgram("assets/shaders/wboit/opaquelinear.vert", "assets/shaders/wboit/opaquelinear.frag")),
	combineshader(ShaderProgram::createShaderProgram("assets/shaders/wboit/combineshader.vert", "assets/shaders/wboit/combineshader.frag")),
	accumulationshader(ShaderProgram::createShaderProgram("assets/shaders/wboit/accumulationshader.vert", "assets/shaders/wboit/accumulationshader.frag")),
	m_width(width),
	m_height(height)
{
	if (samples == 0)
	{
		//create opaque rendertargets and framebuffer
		glGenTextures(1, &rt_opaque_color);
		glBindTexture(GL_TEXTURE_2D, rt_opaque_color);
		glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA32F, width, height);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

		glGenRenderbuffers(1, &rb_opaque_depth);
		glBindRenderbuffer(GL_RENDERBUFFER, rb_opaque_depth);
		glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT32F, width, height);

		glGenFramebuffers(1, &fb_o);
		glBindFramebuffer(GL_FRAMEBUFFER, fb_o);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, rt_opaque_color, 0);
		glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rb_opaque_depth);
		GLenum db[]{ GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1 };
		glDrawBuffers(1, db);

		//create transparent framebuffers and rendertargets

		glGenTextures(1, &rt_t_coloracc);
		glBindTexture(GL_TEXTURE_2D, rt_t_coloracc);
		glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA32F, width, height);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

		glGenTextures(1, &rt_t_revealage);
		glBindTexture(GL_TEXTURE_2D, rt_t_revealage);
		glTexStorage2D(GL_TEXTURE_2D, 1, GL_R32F, width, height);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

		glGenFramebuffers(1, &fb_t);
		glBindFramebuffer(GL_FRAMEBUFFER, fb_t);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, rt_t_coloracc, 0);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, rt_t_revealage, 0);
		glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rb_opaque_depth);
		glDrawBuffers(2, db);
	}
	else
	{

	}
}

WeightedBlendedOITRenderer::~WeightedBlendedOITRenderer()
{
	if (fb_o)
		glDeleteFramebuffers(1, &fb_o);
	if (fb_t)
		glDeleteFramebuffers(1, &fb_t);

	if (rb_opaque_depth)
		glDeleteRenderbuffers(1, &rb_opaque_depth);

	if (rt_opaque_color)
		glDeleteTextures(1, &rt_opaque_color);
	if (rt_t_coloracc)
		glDeleteTextures(1, &rt_t_coloracc);
	if (rt_t_revealage)
		glDeleteTextures(1, &rt_t_revealage);
}

void WeightedBlendedOITRenderer::render(Scene * scene, double dt, bool measure)
{
	renderOpaque(scene, dt);
	std::chrono::time_point<std::chrono::high_resolution_clock> t1; if(measure) t1 = std::chrono::high_resolution_clock::now();
	renderTransparent(scene, dt);
	combineAndDisplay(scene, dt);
	if(measure) {glFinish(); ltrt = 1e-3 * static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - t1).count());}
}

void WeightedBlendedOITRenderer::renderOpaque(Scene * scene, double dt)
{
	glBindFramebuffer(GL_FRAMEBUFFER, fb_o);
	glViewport(0, 0, m_width, m_height);
	glClearColor(scene->clearColor.x, scene->clearColor.y, scene->clearColor.z, scene->clearColor.w);
	glClearDepth(1.0f);
	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glDepthMask(GL_TRUE);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	opaqueshader.use();
	scene->m_camera.bind(&opaqueshader, "camera");

	for (size_t i = 0; i < scene->m_dirlights.size(); ++i)
		scene->m_dirlights[i].bind(&opaqueshader, scene->m_camera.getViewMatrix(), ("dirlights[" + std::to_string(i) + "]").c_str());
	opaqueshader.setUniform("dirlightcount", static_cast<GLint>(scene->m_dirlights.size()));

	for (size_t i = 0; i < scene->m_pointlights.size(); ++i)
		scene->m_pointlights[i].bind(&opaqueshader, scene->m_camera.getViewMatrix(), ("pointlights[" + std::to_string(i) + "]").c_str());
	opaqueshader.setUniform("pointlightcount", static_cast<GLint>(scene->m_pointlights.size()));

	for (size_t i = 0; i < scene->m_ambientlights.size(); ++i)
		scene->m_ambientlights[i].bind(&opaqueshader, ("ambientlights[" + std::to_string(i) + "]").c_str());
	opaqueshader.setUniform("ambientlightcount", static_cast<GLint>(scene->m_ambientlights.size()));

	//opaqueshader.setUniform("normalsfromtex", scene->enableNormalMapping);

	scene->drawOpaqueWithMaterials(&opaqueshader);
}

void WeightedBlendedOITRenderer::renderTransparent(Scene * scene, double dt)
{
	glBindFramebuffer(GL_FRAMEBUFFER, fb_t);
	glViewport(0, 0, m_width, m_height);

	//clear coloracc to (0, 0, 0, 0)
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	GLenum dbc[] { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1 };
	glDrawBuffers(1, &dbc[0]);
	glClear(GL_COLOR_BUFFER_BIT);
	//clear revealage to (1)
	glClearColor(1.0f, 0.0f, 0.0f, 0.0f);
	glDrawBuffers(1, &dbc[1]);
	glClear(GL_COLOR_BUFFER_BIT);
	//reset drawbuffers
	glDrawBuffers(2, dbc);	
	
	glEnable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glDepthMask(GL_FALSE);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);

	glBlendEquation(GL_FUNC_ADD);
	glBlendFuncSeparatei(0, GL_ONE, GL_ONE, GL_ONE, GL_ONE);
	glBlendFuncSeparatei(1, GL_ZERO, GL_ONE_MINUS_SRC_COLOR, GL_ZERO, GL_ZERO);

	//render transparent geometry
	accumulationshader.use();

	scene->m_camera.bind(&accumulationshader, "camera");

	for (size_t i = 0; i < scene->m_dirlights.size(); ++i)
		scene->m_dirlights[i].bind(&accumulationshader, scene->m_camera.getViewMatrix(), ("dirlights[" + std::to_string(i) + "]").c_str());
	accumulationshader.setUniform("dirlightcount", static_cast<GLint>(scene->m_dirlights.size()));

	for (size_t i = 0; i < scene->m_pointlights.size(); ++i)
		scene->m_pointlights[i].bind(&accumulationshader, scene->m_camera.getViewMatrix(), ("pointlights[" + std::to_string(i) + "]").c_str());
	accumulationshader.setUniform("pointlightcount", static_cast<GLint>(scene->m_pointlights.size()));

	for (size_t i = 0; i < scene->m_ambientlights.size(); ++i)
		scene->m_ambientlights[i].bind(&accumulationshader, ("ambientlights[" + std::to_string(i) + "]").c_str());
	accumulationshader.setUniform("ambientlightcount", static_cast<GLint>(scene->m_ambientlights.size()));

	//accumulationshader.setUniform("normalsfromtex", scene->enableNormalMapping);
	accumulationshader.setUniform("usedepthfunc", scene->wboitusedepthfunc);

	scene->drawTransparentWithMaterials(&accumulationshader);
}

void WeightedBlendedOITRenderer::combineAndDisplay(Scene * scene, double dt)
{
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glViewport(0, 0, m_width, m_height);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glDisable(GL_BLEND);
	glDisable(GL_DEPTH_TEST);
	glDepthMask(GL_FALSE);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glClear(GL_COLOR_BUFFER_BIT);

	combineshader.use();
	combineshader.setUniform("tmwhite", scene->maxlum);

	//bind fg and bg textures
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, rt_opaque_color);
	combineshader.setUniform("opqbuf", 0);

	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, rt_t_coloracc);
	combineshader.setUniform("coloracc", 1);

	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, rt_t_revealage);
	combineshader.setUniform("revealage", 2);

	Primitives::drawNDCQuad();
}

double WeightedBlendedOITRenderer::getLastTransparentRenderTime()
{
	return ltrt;
}
