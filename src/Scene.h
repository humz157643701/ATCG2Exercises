/** \addtogroup scene
Scene class
*  @{
*/

/*!
\file Scene.h
*/

#ifndef _SCENE_H_
#define _SCENE_H_
#include <vector>
#include <Mesh.h>
#include <Model.h>
#include <Material.h>
#include <Lights.h>
#include <Camera.h>
#include <IRenderer.h>
#include <ParticleEmitter.h>

/*!
\brief Holds scene description
*/
class Scene
{
public:
	//! Camera
	Camera m_camera;
	//! All loaded textures
	std::vector<std::unique_ptr<Texture>> m_textures;
	//! All materials used in the scene
	std::vector<std::unique_ptr<Material>> m_materials;
	//! All meshes used in the scene
	std::vector<std::unique_ptr<Mesh>> m_meshes;
	//! All opaque models
	std::vector<Model> m_opaque_models;
	//! All transparent models
	std::vector<Model> m_transparent_models;
	//! All directional lights
	std::vector<DirectionalLight> m_dirlights;
	//! ALl point lights
	std::vector<PointLight> m_pointlights;
	//! All ambient lights
	std::vector<AmbientLight> m_ambientlights;
	//! Maximum luminance at which the tone mapping step outputs white
	float maxlum;
	//! Backgroud color / clear color
	glm::vec4 clearColor;
	//! Environment map
	std::unique_ptr<Texture> m_env_map;

	/**
	\brief Draws every opaque model without setting material uniforms
	\param shader Shader being used for drawing
	*/
	void drawOpaque(ShaderProgram* shader);
	/**
	\brief Draws every opaque model with setting material uniforms
	\param shader Shader being used for drawing
	*/
	void drawOpaqueWithMaterials(ShaderProgram* shader);
	/**
	\brief Draws every transparent model without setting material uniforms
	\param shader Shader being used for drawing
	*/
	void drawTransparent(ShaderProgram* shader);
	/**
	\brief Draws every transparent model with setting material uniforms
	\param shader Shader being used for drawing
	*/
	void drawTransparentWithMaterials(ShaderProgram* shader);
	/**
	\brief Triggers a simulation step for every particle system
	\param dt Delta time for the integrator to advance the system
	*/
	void update(float dt);
};

#endif

/** @}*/