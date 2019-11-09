/** \addtogroup renderers
LLI2_DB renderer
*  @{
*/

/*!
\file INDLLWBT.h
*/

#ifndef _LINKED_LIST_WBT_V2_RENDERER_H_
#define _LINKED_LIST_WBT_V2_RENDERER_H_
#include <IRenderer.h>
#include <Shader.h>
#include <Framebuffer.h>
#include <vector>

//! A hard coded maximum of bins, neccessary for buffer creation.
#define INDLLWBTV2_MAX_BINS 64u

#define INDLLWBTV2_HISTOGRAM_BINS 32

/**
\brief Implements the linked list renderer with interpolation, normalization and depth buffer offset
*/
class INDLLWBTV2 : public IRenderer
{
public:
	/*!
	\brief Constructs the renderer.
	\param width Horizontal target resolution
	\param height Vertical target resolution
	\param samples Currently unused
	\param max_elements Maximum size of the fragment buffer
	*/
	INDLLWBTV2(GLsizei width, GLsizei height, GLsizei samples, GLsizei max_elements);
	//! Destructor
	~INDLLWBTV2();
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

	// shaders
	ShaderProgram opaqueshader;
	ShaderProgram depthshader;
	ShaderProgram accumulationshader;
	ShaderProgram blendshader;
	ShaderProgram combineshader;

	// opaque buffers
	GLuint fb_o;
	GLuint rt_opaque_color;
	GLuint rb_opaque_depth;

	// prepass: dual depth, frag count
	GLuint fb_d;
	GLuint rt_transparent_min_max_depth;
	GLuint rt_frag_count;

	// render targets for depth histograms (wednesday)
	
	GLuint fb_depthhist;
	GLuint rt_depth_hist[INDLLWBTV2_HISTOGRAM_BINS / 4];

	// fragment list
	GLuint fb_acc;
	GLuint im_offsetbuffer;
	GLuint b_fragmentdata;
	GLuint b_llcounter;
	GLuint m_max_elements;
	GLsizeiptr fragbufsize;

	// bins
	//GLuint b_bucketbuffer;
	GLuint numbins;

	// final transparent buffer
	GLuint fb_blend;
	GLuint rt_blend;

	GLsizei m_width;
	GLsizei m_height;
	GLsizei m_samples;

	double ltrt;
};

#endif

/** @}*/