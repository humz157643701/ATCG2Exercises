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

/**
\brief Implements a ward brdf renderer
*/
class WardRenderer : public IRenderer
{
public:
	/*!
	\brief Constructs the renderer.
	*/
	WardRenderer();
	//! Destructor
	~WardRenderer();
	/*!
	\brief Renders the scene.

	Renders the scene. Measurements are not implemented for this renderer, regardless of the last parameter's value.

	\param scene Pointer to a scene to be rendered
	\param dt Last frame time
	\param measure Unused in this implementation
	*/
	void render(Scene * scene, double dt, bool measure = false) override;
	//! Returns 0.0
	virtual double getLastTransparentRenderTime() override;

private:
	ShaderProgram m_opaque_shader;
	//ShaderProgram transparent_shader;
};

#endif