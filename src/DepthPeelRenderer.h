/** \addtogroup renderers
Depth peeling renderer
*  @{
*/

/*!
\file DepthPeeling.h
*/

#ifndef _DEPTH_PEEL_RENDERER_H_
#define _DEPTH_PEEL_RENDERER_H_
#include <IRenderer.h>
#include <Shader.h>
#include <Framebuffer.h>

/**
	\brief Implements the depth peeling algorithm
*/
class DepthPeelRenderer : public IRenderer
{
public:
	/*!
	\brief Constructs a depth peeling renderer.
	\param width Horizontal target resolution
	\param height Vertical target resolution
	\param samples Currently unused
	*/
	DepthPeelRenderer(GLsizei width, GLsizei height, GLsizei samples);
	//! Destructor
	~DepthPeelRenderer();

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
	ShaderProgram peelshader;
	ShaderProgram combineshader;
	ShaderProgram blendshader;

	GLuint fb_o;
	GLuint rt_opaque_color;
	GLuint rb_opaque_depth;

	GLuint fb_t1;
	GLuint fb_t2;
	GLuint fb_blend;
	
	GLuint rt_t_color;
	GLuint rt_t_acc;
	GLuint rt_t1_depth;
	GLuint rt_t2_depth;

	GLsizei m_width;
	GLsizei m_height;
	GLsizei m_samples;

	GLuint m_occq;

	double ltrt;
	
};

#endif

/** @}*/