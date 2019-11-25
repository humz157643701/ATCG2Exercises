/** \addtogroup renderers
Renderer interface
*  @{
* /

/*!
\file IRenderer.h
*/

#ifndef _IRENDERER_H_
#define _IRENDERER_H_
#include <libheaders.h>

class Scene;
class IRenderer
{
public:
	/*!
	\brief Renders the scene.

	Renders the scene. The last parameter, if true, causes a measurement to be triggered.
	The time that it takes to render the transparent part and merge it with the opaque part is recorded. After this function returns, the result can be read via getLastTransparentRenderTime().
	During the measurement glFinish() is called to wait for the completion of all running opengl operations, which is slow.
	Thus, don't set this parameter true unless you want to do an actual measurement.

	\param scene Pointer to a scene to be rendered
	\param dt Last frame time
	\param measure Enables/disables the measurement for this frame
	*/
	virtual void render(Scene* scene, double dt, bool measure = false, bool clear = true, GLuint fbo = 0) = 0;
	/*!
	\brief Returns the last measured time for transparent rendering and merging in milliseconds
	*/
	virtual double getLastTransparentRenderTime() = 0;
	virtual size_t rendererid() const = 0;
	virtual ~IRenderer() {};
};

#endif

/** @}*/