/** \addtogroup renderers
WBOIT renderer
*  @{
* /

/*!
\file WeightedBlendedOITRenderer.h
*/

#ifndef _WBOIT_RENDERER_H_
#define _WBOIT_RENDERER_H_
#include <IRenderer.h>
#include <Shader.h>
#include <Framebuffer.h>

/**
\brief Implements the weighted blended OIT renderer
*/
class WeightedBlendedOITRenderer : public IRenderer
{
public:
	/*!
	\brief Constructs the renderer.
	\param width Horizontal target resolution
	\param height Vertical target resolution
	\param samples Currently unused
	*/
	WeightedBlendedOITRenderer(GLsizei width, GLsizei height, GLsizei samples);
	//! Destructor
	~WeightedBlendedOITRenderer();
	/*!
	\brief Renders the scene.

	Renders the scene. The last parameter, if true, causes a measurement to be triggered.
	The time that it takes to render the transparent part and merge it with the opaque part is recorded. After this function returns, the result can be read via getLastTransparentRenderTime().
	During the measurement glFinish() is called to wait for the completion of all running opengl operations, which is slow.
	Thus, don't set this parameter true unless you wan't to do an actual measurement.

	\param scene Pointer to a scene to be rendered
	\param dt Last frame time
	\param measure Enables/disables the measurement for this frame
	*/
	void render(Scene * scene, double dt, bool measure) override;
	/*!
	\brief Returns the last measured time for transparent rendering and merging in milliseconds
	*/
	virtual double getLastTransparentRenderTime() override;

private:
	void renderOpaque(Scene* scene, double dt);
	void renderTransparent(Scene* scene, double dt);
	void combineAndDisplay(Scene* scene, double dt);


	ShaderProgram opaqueshader;
	ShaderProgram combineshader;
	ShaderProgram accumulationshader;

	GLuint fb_o;
	GLuint rt_opaque_color;
	GLuint rb_opaque_depth;

	GLuint fb_t;	
	GLuint rt_t_coloracc; //weighted color sum in rgb, sum(ai * wi) in a
	GLuint rt_t_revealage; //prod(1 - ai)

	GLsizei m_width;
	GLsizei m_height;
	GLsizei m_samples;

	double ltrt;

	
};

#endif