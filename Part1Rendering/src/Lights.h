/** \addtogroup scene_components
Various light classes.
*  @{
* /

/*!
\file Lights.h
*/

#ifndef _LIGHTS_H_
#define _LIGHTS_H_
#include <libheaders.h>
#include <Transform.h>
#include <Shader.h>

/**
\brief Represents a directional light source in the scene
*/
class DirectionalLight
{
public:
	//! World transformation
	Transform transform;
	//! Light color
	glm::vec3 color;
	/*!
	\brief Sets the light-related uniform parameters for the given shader program
	To use the light parameters in a shader, define struct like follows and add a uniform variable of this type.
	\code
	struct DirectionalLight
	{
		float lumint;
		vec3 color;
		vec3 direction;
	};
	uniform DirectionalLight dirlight;
	\endcode
	\returns True if all parameters were set successfully
	\param shader Shader program to set the uniform parameters for
	\param mat	Transformation matrix to transform light parameters into another coordinate space before setting the uniforms
	\param instanceName	Name of the light struct instance in the shader
	*/
	bool bind(ShaderProgram* shader, const glm::mat4& mat, const char* instancename);
	/*!
	\brief Initializes a light instance
	\param _transform	World transformation of the light
	\param _li			Light intensity
	\param _color		Light color
	*/
	DirectionalLight(const Transform& _transform, const glm::vec3& _color);
	/*!
	\brief Initializes a light instance with some default parameters
	*/
	DirectionalLight();
	//! Copy ctor
	DirectionalLight(const DirectionalLight&) = default;
	//! Move ctor
	DirectionalLight(DirectionalLight&&) = default;
	//! Copy assignment operator
	DirectionalLight& operator=(const DirectionalLight&) = default;
	//! Move assignment operator
	DirectionalLight& operator=(DirectionalLight&&) = default;
};

/**
\brief Represents a point light source in the scene
*/
class PointLight
{
public:
	//! World transformation
	Transform transform;
	//! Light color
	glm::vec3 color;
	/*!
	\brief Sets the light-related uniform parameters for the given shader program
	To use the light parameters in a shader, define struct like follows and add a uniform variable of this type.
	\code
	struct PointLight
	{
		float lumint;
		vec3 color;
		vec3 position;
	};
	uniform PointLight pointlight;
	\endcode
	\returns True if all parameters were set successfully
	\param shader Shader program to set the uniform parameters for
	\param mat	Transformation matrix to transform light parameters into another coordinate space before setting the uniforms
	\param instanceName	Name of the light struct instance in the shader
	*/
	bool bind(ShaderProgram* shader, const glm::mat4& mat, const char* instancename);
	/*!
	\brief Initializes a light instance
	\param _transform	World transformation of the light
	\param _li			Light intensity
	\param _color		Light color
	*/
	PointLight(const Transform& _transform, const glm::vec3& _color);
	/*!
	\brief Initializes a light instance with some default parameters
	*/
	PointLight();
	//! Copy ctor
	PointLight(const PointLight&) = default;
	//! Move ctor
	PointLight(PointLight&&) = default;
	//! Copy assignment operator
	PointLight& operator=(const PointLight&) = default;
	//! Move assignment operator
	PointLight& operator=(PointLight&&) = default;
};

/**
\brief Represents an ambient light source in the scene
*/
class AmbientLight
{
public:
	//! Light color
	glm::vec3 color;
	/*!
	\brief Sets the light-related uniform parameters for the given shader program
	To use the light parameters in a shader, define struct like follows and add a uniform variable of this type.
	\code
	struct AmbientLight
	{
		float lumint;
		vec3 color;
	};
	uniform PointLight pointlight;
	\endcode
	\returns True if all parameters were set successfully
	\param shader Shader program to set the uniform parameters for
	\param instanceName	Name of the light struct instance in the shader
	*/
	bool bind(ShaderProgram* shader, const char* instancename);
	/*!
	\brief Initializes a light instance
	\param _li			Light intensity
	\param _color		Light color
	*/
	AmbientLight(const glm::vec3& _color);
	/*!
	\brief Initializes a light instance with some default parameters
	*/
	AmbientLight();
	//! Copy ctor
	AmbientLight(const AmbientLight&) = default;
	//! Move ctor
	AmbientLight(AmbientLight&&) = default;
	//! Copy assignment operator
	AmbientLight& operator=(const AmbientLight&) = default;
	//! Move assignment operator
	AmbientLight& operator=(AmbientLight&&) = default;
};

#endif


/** @}*/