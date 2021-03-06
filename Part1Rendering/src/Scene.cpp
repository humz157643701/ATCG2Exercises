#include "Scene.h"

void Scene::drawOpaque(ShaderProgram * shader, size_t rid)
{
	for (auto& m :m_opaque_models)
	{
		m.draw(shader, rid);
	}
}

void Scene::drawOpaqueWithMaterials(ShaderProgram * shader, size_t rid)
{
	for (auto& m : m_opaque_models)
	{
		m.drawWithMaterials(shader, rid);
	}
}

void Scene::drawTransparent(ShaderProgram * shader, size_t rid)
{
	for (auto& m : m_transparent_models)
	{
		m.draw(shader, rid);
	}
}

void Scene::drawTransparentWithMaterials(ShaderProgram * shader, size_t rid)
{
	for (auto& m : m_transparent_models)
	{
		m.drawWithMaterials(shader, rid);
	}
}

void Scene::update(float dt)
{
	
}

Transform* Scene::getTransformByOID(std::size_t oid)
{
	if(m_camera.transform.getOID() == oid)
		return &m_camera.transform;
	for(Model& m : m_opaque_models)
	{
		if(m.transform.getOID() == oid)
			return &m.transform;
	}
	for(Model& m : m_transparent_models)
	{
		if(m.transform.getOID() == oid)
			return &m.transform;
	}
	for(PointLight& p : m_pointlights)
	{
		if(p.transform.getOID() == oid)
			return &p.transform;
	}
	return nullptr;
}

Material* Scene::getMaterialByID(std::size_t mid)
{
	for(std::unique_ptr<Material>& m : m_materials)
		if(m->m_mid == mid)
			return m.get();
	return nullptr;
}
