/** \addtogroup renderers
Ward renderer
*  @{
* /

/*!
\file WardRenderer.h
*/

#ifndef _WARD_RENDERER_H_
#define _WARD_RENDERER_H_
#include <IRenderer.h>
#include <Shader.h>
#include <Scene.h>

/**
\brief Implements a ward brdf renderer
*/
class WardRenderer : public IRenderer
{
public:
	/*!
	\brief Constructs the renderer.
	*/
	WardRenderer(Scene* scn);
	//! Destructor
	~WardRenderer();
	/*!
	\brief Renders the scene.

	Renders the scene. Measurements are not implemented for this renderer, regardless of the last parameter's value.

	\param scene Pointer to a scene to be rendered
	\param dt Last frame time
	\param measure Unused in this implementation
	*/
	void render(Scene * scene, double dt, bool measure = false, bool clear= true, GLuint fbo = 0) override;
	//! Returns 0.0
	virtual double getLastTransparentRenderTime() override;

	virtual size_t rendererid() const  override
	{
		return 1;
	}

private:
	ShaderProgram m_opaque_shader;
	ShaderProgram m_ermap_to_cubemap;
	ShaderProgram m_skybox_shader;
	//ShaderProgram transparent_shader;
	std::unique_ptr<Texture> m_skybox;
	std::unique_ptr<Texture> m_irradiance_map;

	// framebuffers
	GLuint m_fb_erconv;
};

#endif