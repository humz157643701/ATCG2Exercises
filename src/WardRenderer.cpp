#include "WardRenderer.h"
#include <Primitives.h>

WardRenderer::WardRenderer(Scene* scn) :
	m_opaque_shader(ShaderProgram::createShaderProgram("assets/shaders/ward/ward_opaque.vert", "assets/shaders/ward/ward_opaque.frag")),
	m_ermap_to_cubemap(ShaderProgram::createShaderProgram("assets/shaders/ward/envconv.vert", "assets/shaders/ward/envconv.frag", "assets/shaders/ward/envconv.geom")),
	m_skybox_shader(ShaderProgram::createShaderProgram("assets/shaders/ward/skybox.vert", "assets/shaders/ward/skybox.frag"))
{
	glEnable(GL_TEXTURE_CUBE_MAP_SEAMLESS);
	// convert er environment map to cube map
	m_skybox = Texture::TCB(GL_RGB32F, scn->m_skyboxres, scn->m_skyboxres, Texture::calculateMipMapLevels(scn->m_skyboxres, scn->m_skyboxres));
	
	glGenFramebuffers(1, &m_fb_erconv);
	glBindFramebuffer(GL_FRAMEBUFFER, m_fb_erconv);
	glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, m_skybox->tex, 0);

	// projection and view matrices for env map converison
	glm::mat4 pmat = glm::perspective(glm::radians(90.0f), 1.0f, 0.1f, 10.0f);
	glm::mat4 layermats[] = {  //px, nx, py, ny, pz, nz
		// due to weird opengl cube map sampling, we're rendering x and z faces rotated
		pmat * glm::lookAt(glm::vec3(0.0f), glm::vec3(1.0f, 0.0f, 0.0f), glm::vec3(0.0f, -1.0f, 0.0f)),
		pmat * glm::lookAt(glm::vec3(0.0f), glm::vec3(-1.0f, 0.0f, 0.0f), glm::vec3(0.0f, -1.0f, 0.0f)),
		pmat * glm::lookAt(glm::vec3(0.0f), glm::vec3(0.0f, 1.0f, 0.0f), glm::vec3(0.0f, 0.0f, 1.0f)),
		pmat * glm::lookAt(glm::vec3(0.0f), glm::vec3(0.0f, -1.0f, 0.0f), glm::vec3(0.0f, 0.0f, -1.0f)),
		pmat * glm::lookAt(glm::vec3(0.0f), glm::vec3(0.0f, 0.0f, 1.0f), glm::vec3(0.0f, -1.0f, 0.0f)),
		pmat * glm::lookAt(glm::vec3(0.0f), glm::vec3(0.0f, 0.0f, -1.0f), glm::vec3(0.0f, -1.0f, 0.0f))
	};

	glDisable(GL_DEPTH_TEST);
	glDepthMask(GL_FALSE);
	glViewport(0, 0, scn->m_skyboxres, scn->m_skyboxres);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	m_ermap_to_cubemap.use();

	for (size_t i = 0; i < 6; i++)
		m_ermap_to_cubemap.setUniform(("u_layer_matrices[" + std::to_string(i) + "]").c_str(), layermats[i], false);

	m_ermap_to_cubemap.bindTex("u_enver", scn->m_er_env_map.get());
	Primitives::drawNDCCube();
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	m_skybox->bind(0);
	m_skybox->generateMipMap();
	
	float max_aniso;
	glGetFloatv(GL_MAX_TEXTURE_MAX_ANISOTROPY, &max_aniso);
	m_skybox->setTexParam(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	m_skybox->setTexParam(GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	m_skybox->setTexParam(GL_TEXTURE_WRAP_S, GL_REPEAT);
	m_skybox->setTexParam(GL_TEXTURE_WRAP_T, GL_REPEAT);
	m_skybox->setTexParamArr(GL_TEXTURE_MAX_ANISOTROPY, &max_aniso, 1);

	m_skybox->unbind();
	
	glDeleteFramebuffers(1, &m_fb_erconv);
}

WardRenderer::~WardRenderer()
{
}

void WardRenderer::render(Scene * scene, double dt, bool measure, bool clear, GLuint fbo)
{
	glBindFramebuffer(GL_FRAMEBUFFER, fbo);
	glEnable(GL_SCISSOR_TEST);
	scene->m_camera.updateGLViewport();
	scene->m_camera.updateGLScissor();
	glClearColor(scene->clearColor.x, scene->clearColor.y, scene->clearColor.z, scene->clearColor.w);
	glClearDepth(1.0f);
	glEnable(GL_CULL_FACE);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// for rendering the environment map disable depth writes
	glDepthMask(GL_FALSE);
	glDisable(GL_DEPTH_TEST);
	// we're rendering the inside of a cube, so cull front faces
	glCullFace(GL_BACK); // excuse me, WTF?!

	// render environment map
	m_skybox_shader.use();
	m_skybox_shader.bindTex("u_skybox", m_skybox.get());
	m_skybox_shader.setUniform("u_view_matrix", scene->m_camera.getViewRotationMatrix(), false);
	m_skybox_shader.setUniform("u_projection_matrix", scene->m_camera.getProjectionMatrix(), false);
	m_skybox_shader.setUniform("exposure", scene->tm_exposure);
	Primitives::drawNDCCube();

	// render opaque geometry
	// enable depth writes
	glEnable(GL_DEPTH_TEST);
	glDepthMask(GL_TRUE);
	// cull back faces
	glCullFace(GL_BACK);

	m_opaque_shader.use();
	m_opaque_shader.bindTex("skybox", m_skybox.get());
	m_opaque_shader.setUniform("skybox_res", static_cast<float>(scene->m_skyboxres));
	m_opaque_shader.setUniform("skybox_lodlevels", static_cast<float>(m_skybox->miplevels));
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

	scene->drawOpaqueWithMaterials(&m_opaque_shader, rendererid());
}

double WardRenderer::getLastTransparentRenderTime()
{
	return 0.0;
}
