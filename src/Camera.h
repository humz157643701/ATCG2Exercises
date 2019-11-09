/** \addtogroup scene_components
Implements a standard, 1st person fly-through camera with perspective projection
*  @{
*/

/*!
\file Camera.h
*/

#ifndef _CAMERA_H_
#define _CAMERA_H_
#include <Transform.h>
#include <libheaders.h>
#include <Shader.h>
#include <vector>

/*!
\brief Implements a standard, 1st person fly-through camera with perspective projection
*/
class Camera
{
	float m_aspectratio;
	float m_width;
	float m_height;
	float m_fovy;
	float m_near;
	float m_far;

	glm::mat4 m_viewMatrix;
	glm::mat4 m_projectionMatrix;

	bool m_viewDirty;
	bool m_projDirty;

public:
	/*!
	\brief Constructs a camera with come default settings
	*/
	Camera();
	/*!
	\brief Constructs a camera
	\param width		Horizontal resolution of the target image
	\param height		Vertical resolution of the target image
	\param fovy			Vertical field of view in radians
	\param near			Position of the near plane. Should be greater than 0.0
	\param far			Position of the far plane. Should be greater than near
	*/
	Camera(int width, int height, float fovy, float near, float far);
	//! Destructor
	~Camera();

	/*!
	\brief Returns the current view matrix	
	\returns Returns the current view matrix
	*/
	const glm::mat4& getViewMatrix();
	/*!
	\brief Returns the current projection matrix
	\returns Returns the current projection matrix
	*/
	const glm::mat4& getProjectionMatrix();

	/*!
	\brief Sets the camera-related uniform parameters for the given shader program
	To use the camera parameters in a shader, define struct like follows and add a uniform variable of this type.
	\code
	struct Camera
	{
		mat4 viewMatrix;
		mat4 projectionMatrix;
		vec3 worldPos;
		float renderWidth;
		float renderHeight;
		float nearPlane;
		float farPlane;
	};
	uniform Camera camera;
	\endcode
	\returns True if all parameters were set successfully
	\param shader	Shader program to set the uniform parameters for
	\param instanceName	Name of the camera struct instance in the shader
	*/
	bool bind(ShaderProgram* shader, const std::string& instanceName = "camera");
	/*!
	\brief Moves the camera forwards into the current view direction
	\param amount		Amount of movement
	*/
	void forward(float amount);
	/*!
	\brief Moves the camera backwards into the current negative view direction
	\param amount		Amount of movement
	*/
	void backward(float amount);
	/*!
	\brief Moves the camera left
	\param amount		Amount of movement
	*/
	void left(float amount);
	/*!
	\brief Moves the camera right
	\param amount		Amount of movement
	*/
	void right(float amount);
	/*!
	\brief Moves the camera up on the global y-axis
	\param amount		Amount of movement
	*/
	void up(float amount);
	/*!
	\brief Moves the camera down on the global y-axis
	\param amount		Amount of movement
	*/
	void down(float amount);
	/*!
	\brief Rotates the camera's view
	\param Pitch Pitch angle in radians
	\param Yaw angle in radians
	*/
	void rotateView(float pitch, float yaw);

	/*!
	\brief Changes the camera's target resolution and aspect ratio
	\param width Desired horizontal resolution
	\param height Desired vertical resolution
	*/
	void setViewport(int width, int height);
	/*!
	\brief Changes the position of the camera's near and far planes
	\param near Desired near plane position
	\param height Desired far plane position
	*/
	void setFrustumBounds(float near, float far);
	/*!
	\brief Changes the camera's field of view
	\param fovy Desired vertical field of view in radians
	*/
	void setFieldOfView(float fovy);
	/*!
	\brief Changes all camera parameters at once
	\param width Desired horizontal resolution
	\param height Desired vertical resolution
	\param height Desired vertical field of view in radians
	\param height Desired near plane position
	\param height Desired far plane position
	*/
	void setParameters(int width, int height, float fovy, float near, float far);

	/*!
	\brief Returns the current horizontal resolution
	\returns Returns the current horizontal resolution
	*/
	float getWidth();
	/*!
	\brief Returns the current vertical resolution
	\returns Returns the current vertical resolution
	*/
	float getHeight();
	/*!
	\brief Returns the current near plane
	\returns Returns the current near plane
	*/
	float getNear();
	/*!
	\brief Returns the current far plane
	\returns Returns the current far plane
	*/
	float getFar();
	/*!
	\brief Returns the current vertical field of view in radians
	\returns Returns the current vertical field of view in radians
	*/
	float getFieldOfView();
	/*!
	\brief Returns the current aspect ratio
	\returns Returns the current aspect ratio
	*/
	float getAspectRatio();


	void setViewDirty();

	//! Current world transformation of the camera
	Transform transform;
};

#endif

/** @}*/