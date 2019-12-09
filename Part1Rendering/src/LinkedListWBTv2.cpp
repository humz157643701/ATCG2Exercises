#include "LinkedListWBTv2.h"
#include <Scene.h>
#include <Primitives.h>
#include <cmath>
#include <chrono>

//size of one fragment element?
//RGBA|D|N
//4444|4|4
/*
struct frag
{
vec4 color;    16 <= alignment
float depth;    4
int next;       4
pad				8
};	=> size is 32
*/
#define FRAG_ELEMENT_SIZE 32

INDLLWBTV2::INDLLWBTV2(GLsizei width, GLsizei height, GLsizei samples, GLsizei max_elements) :
	opaqueshader(ShaderProgram::createShaderProgram("assets/shaders/INDLLWBTV2/opaquelinear.vert", "assets/shaders/INDLLWBTV2/opaquelinear.frag")),
	depthshader(ShaderProgram::createShaderProgram("assets/shaders/INDLLWBTV2/depth.vert", "assets/shaders/INDLLWBTV2/depth.frag")),
	combineshader(ShaderProgram::createShaderProgram("assets/shaders/INDLLWBTV2/combineshader.vert", "assets/shaders/INDLLWBTV2/combineshader.frag")),
	accumulationshader(ShaderProgram::createShaderProgram("assets/shaders/INDLLWBTV2/accumulationshader.vert", "assets/shaders/INDLLWBTV2/accumulationshader.frag")),
	blendshader(ShaderProgram::createShaderProgram("assets/shaders/INDLLWBTV2/blendshader.vert", "assets/shaders/INDLLWBTV2/blendshader.frag")),
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
		GLenum db[]{ GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2, GL_COLOR_ATTACHMENT3, GL_COLOR_ATTACHMENT4, GL_COLOR_ATTACHMENT5, GL_COLOR_ATTACHMENT6, GL_COLOR_ATTACHMENT7 };
		glDrawBuffers(1, db);

		//transparent depth prepass
		// min / max depth
		glGenTextures(1, &rt_transparent_min_max_depth);
		glBindTexture(GL_TEXTURE_2D, rt_transparent_min_max_depth);
		glTexStorage2D(GL_TEXTURE_2D, 1, GL_RG32F, width, height);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		// fragment count
		glGenTextures(1, &rt_frag_count);
		glBindTexture(GL_TEXTURE_2D, rt_frag_count);
		glTexStorage2D(GL_TEXTURE_2D, 1, GL_R32F, width, height);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

		glGenFramebuffers(1, &fb_d);
		glBindFramebuffer(GL_FRAMEBUFFER, fb_d);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, rt_transparent_min_max_depth, 0);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, rt_frag_count, 0);
		glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rb_opaque_depth);
		glDrawBuffers(2, db);

		// depth histogram

		for (size_t i = 0; i < (INDLLWBTV2_HISTOGRAM_BINS / 4); ++i)
		{
			glGenTextures(1, &rt_depth_hist[i]);
			glBindTexture(GL_TEXTURE_2D, rt_depth_hist[i]);
			glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA32F, width, height);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		}		

		//create framebuffer without color attachments, only reuse opaque depth buffer
		glGenFramebuffers(1, &fb_acc);
		glBindFramebuffer(GL_FRAMEBUFFER, fb_acc);
		glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rb_opaque_depth);
		for (size_t i = 0; i < (INDLLWBTV2_HISTOGRAM_BINS / 4); ++i)
		{
			glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0 + GLenum(i), GL_TEXTURE_2D, rt_depth_hist[i], 0);
		}
		glDrawBuffers(8, &db[0]);

		//create atomic counter
		glGenBuffers(1, &b_llcounter);
		glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, b_llcounter);
		GLuint cinit{ 0 };
		glBufferStorage(GL_ATOMIC_COUNTER_BUFFER, sizeof(GLuint), &cinit, GL_MAP_WRITE_BIT);
		glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, 0);

		//create offset image
		glGenTextures(1, &im_offsetbuffer);
		glBindTexture(GL_TEXTURE_2D, im_offsetbuffer);
		glTexStorage2D(GL_TEXTURE_2D, 1, GL_R32I, width, height);

		//create storage buffer for fragment data
		glGenBuffers(1, &b_fragmentdata);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, b_fragmentdata);
		glBufferStorage(GL_SHADER_STORAGE_BUFFER, max_elements * FRAG_ELEMENT_SIZE, reinterpret_cast<void*>(0), 0);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, b_fragmentdata);
		fragbufsize = max_elements * FRAG_ELEMENT_SIZE;
		m_max_elements = max_elements;

		//buffer for holding bucket definitions
		// glGenBuffers(1, &b_bucketbuffer);
		// glBindBuffer(GL_SHADER_STORAGE_BUFFER, b_bucketbuffer);
		// glBufferStorage(GL_SHADER_STORAGE_BUFFER, INDLLWBTV2_MAX_BUCKETS * sizeof(GLfloat) * 2, reinterpret_cast<void*>(0), GL_MAP_WRITE_BIT);
		// glBindBuffer(GL_SHADER_STORAGE_BUFFER, b_bucketbuffer);
		numbins = 0;

		//create blend framebuffer and render texture
		glGenTextures(1, &rt_blend);
		glBindTexture(GL_TEXTURE_2D, rt_blend);
		glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA32F, width, height);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

		glGenFramebuffers(1, &fb_blend);
		glBindFramebuffer(GL_FRAMEBUFFER, fb_blend);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, rt_blend, 0);
		glDrawBuffers(1, db);
	}
	else
	{

	}
}

INDLLWBTV2::~INDLLWBTV2()
{
	if (fb_o)
		glDeleteFramebuffers(1, &fb_o);
	if (rb_opaque_depth)
		glDeleteRenderbuffers(1, &rb_opaque_depth);
	if (rt_opaque_color)
		glDeleteTextures(1, &rt_opaque_color);

	if (fb_d)
		glDeleteFramebuffers(1, &fb_d);
	if (rt_transparent_min_max_depth)
		glDeleteTextures(1, &rt_transparent_min_max_depth);

	if (fb_acc)
		glDeleteFramebuffers(1, &fb_acc);
	if (im_offsetbuffer)
		glDeleteTextures(1, &im_offsetbuffer);
	if (b_fragmentdata)
		glDeleteBuffers(1, &b_fragmentdata);
	if (b_llcounter)
		glDeleteBuffers(1, &b_llcounter);

	// if (b_bucketbuffer)
	// 	glDeleteBuffers(1, &b_bucketbuffer);

	if (fb_blend)
		glDeleteFramebuffers(1, &fb_blend);
	if (rt_blend)
		glDeleteTextures(1, &rt_blend);
}

void INDLLWBTV2::render(Scene * scene, double dt, bool measure)
{
	renderOpaque(scene, dt);
	std::chrono::time_point<std::chrono::high_resolution_clock> t1; if(measure) t1 = std::chrono::high_resolution_clock::now();
	renderTransparent(scene, dt);
	combineAndDisplay(scene, dt);
	if(measure) {glFinish(); ltrt = 1e-3 * static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - t1).count());}
}

void INDLLWBTV2::renderOpaque(Scene * scene, double dt)
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

void INDLLWBTV2::renderTransparent(Scene * scene, double dt)
{
	GLenum db[]{ GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1 };
	// --------------------------------- prepass: min, max depth and fragment count ---------------------------------------

	//render dual depth map
	glBindFramebuffer(GL_FRAMEBUFFER, fb_d);
	glViewport(0, 0, m_width, m_height);
	// clear minmax depth buffer
	glDrawBuffers(1, &db[0]);
	glClearColor(-scene->m_camera.getFar(), 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT);
	// clear frag count
	glDrawBuffers(1, &db[1]);
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT);
	// use depth map of opaque pass to discard useless fragments
	glEnable(GL_DEPTH_TEST);
	glDepthMask(GL_FALSE);
	glDepthFunc(GL_LESS);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	// reset draw buffers
	glDrawBuffers(2, &db[0]);
	// set max blend mode for depth maps
	glEnable(GL_BLEND);
	glBlendEquationi(0, GL_MAX);
	glBlendFunci(0, GL_ZERO, GL_ZERO);
	// additive blending for frag count
	glBlendEquationi(1, GL_FUNC_ADD);
	glBlendFunci(1, GL_ONE, GL_ONE);

	depthshader.use();
	scene->m_camera.bind(&depthshader, "camera");
	scene->drawTransparent(&depthshader);

	// -------------------------------------------- fragment lists and depth histogram ---------------------------------------------

	//render accumulation buffer
	glBindFramebuffer(GL_FRAMEBUFFER, fb_acc);
	glViewport(0, 0, m_width, m_height);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glDepthMask(GL_FALSE);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);

	accumulationshader.use();

	// clear depth histogram targets
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT);

	//clear fragment counter and bind to binding point 0
	glMemoryBarrier(GL_BUFFER_UPDATE_BARRIER_BIT);
	glBindBufferBase(GL_ATOMIC_COUNTER_BUFFER, 0, b_llcounter);
	GLuint* buf = static_cast<GLuint*>(glMapBufferRange(GL_ATOMIC_COUNTER_BUFFER, 0, 4, GL_MAP_WRITE_BIT));
	*buf = 0u;
	glUnmapBuffer(GL_ATOMIC_COUNTER_BUFFER);

	//clear offset buffer
	GLint clearval = -1;
	glClearTexSubImage(im_offsetbuffer, 0, 0, 0, 0, m_width, m_height, 1, GL_RED_INTEGER, GL_INT, &clearval);

	//bind offset buffer to image unit
	glBindImageTexture(0, im_offsetbuffer, 0, GL_FALSE, 0, GL_READ_WRITE, GL_R32I);

	//bind fragment buffer to indexed binding target
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, b_fragmentdata);
	accumulationshader.setUniform("max_fragments", static_cast<GLuint>(m_max_elements));

	//render transparent geometry
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

	// for depth histogram
	// bind dual depth map
	auto tu = accumulationshader.getCurrentTU();
	glActiveTexture(GL_TEXTURE0 + tu);
	glBindTexture(GL_TEXTURE_2D, rt_transparent_min_max_depth);
	accumulationshader.setUniform("dual_depth_map", tu);
	accumulationshader.resetTU(tu + 1);
	// activate additive blending on all histogram targets
	glEnable(GL_BLEND);
	for (size_t i = 0; i < (INDLLWBTV2_HISTOGRAM_BINS / 4); ++i)
	{
		glBlendEquationi(i, GL_FUNC_ADD);
		glBlendFunci(i, GL_ONE, GL_ONE);
	}

	scene->drawTransparentWithMaterials(&accumulationshader);

	//guarantee visibility of all the above side effects
	glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT | GL_ATOMIC_COUNTER_BARRIER_BIT | GL_SHADER_STORAGE_BARRIER_BIT);

	//update bucketbuffer if neccessary
	GLuint newnumbins = glm::clamp(static_cast<GLuint>(scene->m_camera.buckets.size()), 1u, INDLLWBTV2_MAX_BINS);
	if (numbins == 0 || numbins != newnumbins)
	{
		// glBindBuffer(GL_SHADER_STORAGE_BUFFER, b_bucketbuffer);
		// GLfloat* buf = static_cast<GLfloat*>(glMapBufferRange(GL_SHADER_STORAGE_BUFFER, 0, newnumbuckets * 2 * sizeof(GLfloat), GL_MAP_WRITE_BIT));
		// for (size_t i = 0; i < newnumbuckets; ++i)
		// {
		// 	buf[i * 2] = scene->m_camera.buckets[i].x;
		// 	buf[i * 2 + 1] = scene->m_camera.buckets[i].y;
		// }
		// glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);
		// glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
		numbins = newnumbins;
	}

	// ------------------------------------- processing the fragment buffers ---------------------------------------
	glBindFramebuffer(GL_FRAMEBUFFER, fb_blend);
	glViewport(0, 0, m_width, m_height);
	glDisable(GL_BLEND);
	glDisable(GL_DEPTH_TEST);
	glDepthMask(GL_FALSE);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);

	blendshader.use();
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, im_offsetbuffer);
	blendshader.setUniform("im_offsetbuffer", 0);
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, b_fragmentdata);
	//glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, b_bucketbuffer);
	blendshader.setUniform("numbins", numbins);

	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, rt_transparent_min_max_depth);
	blendshader.setUniform("depth_map", 1);
	blendshader.setUniform("dwoffset", scene->dwoffset);
	//blendshader.setUniform("proj", glm::vec2(scene->m_camera.getProjectionMatrix()[3][2], scene->m_camera.getProjectionMatrix()[2][2]));

	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, rt_frag_count);
	blendshader.setUniform("frag_count", 2);

	// bind depth histogram textures
	for (std::size_t i = 0; i < INDLLWBTV2_HISTOGRAM_BINS / 4; ++i)
	{
		glActiveTexture(GL_TEXTURE3 + static_cast<GLenum>(i));
		glBindTexture(GL_TEXTURE_2D, rt_depth_hist[i]);
		blendshader.setUniform(("depth_hist" + std::to_string(i)).c_str(), static_cast<GLint>(i + 3));
	}

	Primitives::drawNDCQuad();
}

void INDLLWBTV2::combineAndDisplay(Scene * scene, double dt)
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
	glBindTexture(GL_TEXTURE_2D, rt_blend);
	combineshader.setUniform("tbuf", 1);

	Primitives::drawNDCQuad();
}

double INDLLWBTV2::getLastTransparentRenderTime()
{
	return ltrt;
}
