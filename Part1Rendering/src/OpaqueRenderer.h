/** \addtogroup renderers
Opaque renderer
*  @{
* /

/*!
\file OpaqueRenderer.h
*/

#ifndef _OPAQUE_RENDERER_H_
#define _OPAQUE_RENDERER_H_
#include <IRenderer.h>
#include <Shader.h>

/**
\brief Implements a standard opaque renderer with no transparency but alpha testing
*/
class OpaqueRenderer : public IRenderer
{
public:
	/*!
	\brief Constructs the renderer.
	*/
	OpaqueRenderer();
	//! Destructor
	~OpaqueRenderer();
	/*!
	\brief Renders the scene.

	Renders the scene. Measurements are not implemented for this renderer, regardless of the last parameter's value.

	\param scene Pointer to a scene to be rendered
	\param dt Last frame time
	\param measure Unused in this implementation
	*/
	void render(Scene * scene, double dt, bool measure) override;
	//! Returns 0.0
	virtual double getLastTransparentRenderTime() override;

private:
	ShaderProgram shader;

	
};

#endif