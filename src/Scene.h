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
	//! Exposure value for tone mapping
	float tm_exposure;
	//! Backgroud color / clear color
	glm::vec4 clearColor;
	//! Environment map
	std::unique_ptr<Texture> m_cb_env_map;
	std::unique_ptr<Texture> m_er_env_map;

	GLsizei m_skyboxres;

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
	\brief Triggers a simulation step for every entity that needs framerate-independent simulation
	\param dt Delta time for the integrator to advance the system
	*/
	void update(float dt);

	/**
	\brief Returns a pointer to the transformation of scene entity with object id oid, nullptr otherwise
	\param oid Object id to search for
	*/
	Transform* getTransformByOID(std::size_t oid);

	/**
	\brief Returns a pointer to material with given id, nullptr if no material with this id exists
	\param mid Material id to search for
	*/
	Material* getMaterialByID(std::size_t mid);
};

#endif

/** @}*/