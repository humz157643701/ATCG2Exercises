#include "DepthPeelRenderer.h"
#include <Scene.h>
#include <Primitives.h>
#include <chrono>

DepthPeelRenderer::DepthPeelRenderer(GLsizei width, GLsizei height, GLsizei samples) :
	peelshader(ShaderProgram::createShaderProgram("assets/shaders/depthpeeling/peelshader.vert", "assets/shaders/depthpeeling/peelshader.frag")),
	combineshader(ShaderProgram::createShaderProgram("assets/shaders/depthpeeling/combineshader.vert", "assets/shaders/depthpeeling/combineshader.frag")),
	opaqueshader(ShaderProgram::createShaderProgram("assets/shaders/depthpeeling/opaquelinear.vert", "assets/shaders/depthpeeling/opaquelinear.frag")),
	blendshader(ShaderProgram::createShaderProgram("assets/shaders/depthpeeling/blendshader.vert", "assets/shaders/depthpeeling/blendshader.frag")),
	m_width(width),
	m_height(height),
	fb_o(0),
	rt_opaque_color(0),
	rb_opaque_depth(0),	
	fb_t1(0),
	fb_t2(0),
	fb_blend(0),	
	rt_t_color(0),
	rt_t_acc(0),
	rt_t1_depth(0),
	rt_t2_depth(0)
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
		GLenum db[]{ GL_COLOR_ATTACHMENT0 };
		glDrawBuffers(1, db);

		//create transparent framebuffers and rendertargets
		glGenTextures(1, &rt_t_color);
		glBindTexture(GL_TEXTURE_2D, rt_t_color);
		glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA32F, width, height);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

		glGenTextures(1, &rt_t_acc);
		glBindTexture(GL_TEXTURE_2D, rt_t_acc);
		glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA32F, width, height);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

		glGenTextures(1, &rt_t1_depth);
		glBindTexture(GL_TEXTURE_2D, rt_t1_depth);
		glTexStorage2D(GL_TEXTURE_2D, 1, GL_DEPTH_COMPONENT32F, width, height);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

		glGenTextures(1, &rt_t2_depth);
		glBindTexture(GL_TEXTURE_2D, rt_t2_depth);
		glTexStorage2D(GL_TEXTURE_2D, 1, GL_DEPTH_COMPONENT32F, width, height);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

		glGenFramebuffers(1, &fb_t1);
		glBindFramebuffer(GL_FRAMEBUFFER, fb_t1);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, rt_t_color, 0);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, rt_t1_depth, 0);
		glDrawBuffers(1, db);

		glGenFramebuffers(1, &fb_t2);
		glBindFramebuffer(GL_FRAMEBUFFER, fb_t2);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, rt_t_color, 0);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, rt_t2_depth, 0);
		glDrawBuffers(1, db);

		glGenFramebuffers(1, &fb_blend);
		glBindFramebuffer(GL_FRAMEBUFFER, fb_blend);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, rt_t_acc, 0);
		glDrawBuffers(1, db);

		//occlusion query
		glGenQueries(1, &m_occq);

	}
	else
	{

	}
}

DepthPeelRenderer::~DepthPeelRenderer()
{
	if (fb_o)
		glDeleteFramebuffers(1, &fb_o);
	if (fb_t1)
		glDeleteFramebuffers(1, &fb_t1);
	if (fb_t2)
		glDeleteFramebuffers(1, &fb_t2);
	if (fb_blend)
		glDeleteFramebuffers(1, &fb_blend);

	if (rt_opaque_color)
		glDeleteTextures(1, &rt_opaque_color);
	if (rt_t_color)
		glDeleteTextures(1, &rt_t_color);
	if (rt_t_acc)
		glDeleteTextures(1, &rt_t_acc);
	if (rt_t1_depth)
		glDeleteTextures(1, &rt_t1_depth);
	if (rt_t2_depth)
		glDeleteTextures(1, &rt_t2_depth);

	if (rb_opaque_depth)
		glDeleteRenderbuffers(1, &rb_opaque_depth);

	if (m_occq)
		glDeleteQueries(1, &m_occq);
}

void DepthPeelRenderer::render(Scene * scene, double dt, bool measure)
{
	renderOpaque(scene, dt);
	std::chrono::time_point<std::chrono::high_resolution_clock> t1; if(measure) t1 = std::chrono::high_resolution_clock::now();
	renderTransparent(scene, dt);
	combineAndDisplay(scene, dt);
	if(measure) {glFinish(); ltrt = 1e-3 * static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - t1).count());}
}

void DepthPeelRenderer::renderOpaque(Scene * scene, double dt)
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

void DepthPeelRenderer::renderTransparent(Scene * scene, double dt)
{
	glBindFramebuffer(GL_FRAMEBUFFER, fb_t2);
	glViewport(0, 0, m_width, m_height);
	glClearDepth(0.0f);
	glDepthMask(GL_TRUE);
	glClear(GL_DEPTH_BUFFER_BIT);

	glBindFramebuffer(GL_FRAMEBUFFER, fb_blend);
	glViewport(0, 0, m_width, m_height);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	glDepthFunc(GL_LESS);
	glEnable(GL_CULL_FACE); //assume transparent models to be double-sided
	glCullFace(GL_BACK);

	//blend stuff
	glBlendEquationSeparate(GL_FUNC_ADD, GL_FUNC_ADD);
	glBlendFuncSeparate(GL_DST_ALPHA, GL_ONE, GL_ZERO, GL_ONE_MINUS_SRC_ALPHA);

	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	GLint anyfragmentsdrawn = 1;
	bool enableOcclusionQueries = scene->depthPeelMaxPasses == 0;
	//int layercount = 0;
	for (int i = 0; (i < scene->depthPeelMaxPasses) || (enableOcclusionQueries && anyfragmentsdrawn); ++i)
	{
		//++layercount;
		//peel pass
		peelshader.use();
		glDisable(GL_BLEND);
		glEnable(GL_DEPTH_TEST);
		glDepthMask(GL_TRUE);

		//coose fb and depth texture
		GLuint cwb = (i % 2 == 0 ? fb_t1 : fb_t2);
		GLuint pdepth = (i % 2 == 0 ? rt_t2_depth : rt_t1_depth);

		//blit opaque depth onto current write buffer and initialize color to 0000
		glBindFramebuffer(GL_READ_FRAMEBUFFER, fb_o);
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, cwb);		
		glBlitFramebuffer(0, 0, m_width, m_width, 0, 0, m_width, m_width, GL_DEPTH_BUFFER_BIT, GL_NEAREST);
		glClear(GL_COLOR_BUFFER_BIT);

		//bind depth texture from previous pass
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, pdepth);
		peelshader.setUniform("pdepth", 0);
		peelshader.resetTU(1);

		scene->m_camera.bind(&peelshader, "camera");

		for (size_t i = 0; i < scene->m_dirlights.size(); ++i)
			scene->m_dirlights[i].bind(&peelshader, scene->m_camera.getViewMatrix(), ("dirlights[" + std::to_string(i) + "]").c_str());
		peelshader.setUniform("dirlightcount", static_cast<GLint>(scene->m_dirlights.size()));

		for (size_t i = 0; i < scene->m_pointlights.size(); ++i)
			scene->m_pointlights[i].bind(&peelshader, scene->m_camera.getViewMatrix(), ("pointlights[" + std::to_string(i) + "]").c_str());
		peelshader.setUniform("pointlightcount", static_cast<GLint>(scene->m_pointlights.size()));

		for (size_t i = 0; i < scene->m_ambientlights.size(); ++i)
			scene->m_ambientlights[i].bind(&peelshader, ("ambientlights[" + std::to_string(i) + "]").c_str());
		peelshader.setUniform("ambientlightcount", static_cast<GLint>(scene->m_ambientlights.size()));

		//peelshader.setUniform("normalsfromtex", scene->enableNormalMapping);

		if (enableOcclusionQueries)
			glBeginQuery(GL_ANY_SAMPLES_PASSED, m_occq);

		scene->drawTransparentWithMaterials(&peelshader);
		
		if (enableOcclusionQueries)
		{
			glEndQuery(GL_ANY_SAMPLES_PASSED);
			glGetQueryObjectiv(m_occq, GL_QUERY_RESULT, &anyfragmentsdrawn);
		}

		//blend pass
		blendshader.use();
		glBindFramebuffer(GL_FRAMEBUFFER, fb_blend);
		glEnable(GL_BLEND);
		glDisable(GL_DEPTH_TEST);
		glDepthMask(GL_FALSE);
		
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, rt_t_color);
		blendshader.setUniform("layer", 0);

		Primitives::drawNDCQuad();
	}
	//int dummy = 0;
}

void DepthPeelRenderer::combineAndDisplay(Scene * scene, double dt)
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
	glBindTexture(GL_TEXTURE_2D, rt_t_acc);
	combineshader.setUniform("accbuf", 1);

	Primitives::drawNDCQuad();
}

double DepthPeelRenderer::getLastTransparentRenderTime()
{
	return ltrt;
}
