#include "Scene.h"

void Scene::drawOpaque(ShaderProgram * shader)
{
	for (auto& m :m_opaque_models)
	{
		m.draw(shader);
	}
}

void Scene::drawOpaqueWithMaterials(ShaderProgram * shader)
{
	for (auto& m : m_opaque_models)
	{
		m.drawWithMaterials(shader);
	}
}

void Scene::drawTransparent(ShaderProgram * shader)
{
	for (auto& m : m_transparent_models)
	{
		m.draw(shader);
	}
}

void Scene::drawTransparentWithMaterials(ShaderProgram * shader)
{
	for (auto& m : m_transparent_models)
	{
		m.drawWithMaterials(shader);
	}
}

void Scene::update(float dt)
{
	
}
