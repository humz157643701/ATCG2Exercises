#include "Model.h"

void Model::draw(ShaderProgram* shader, size_t rid)
{
	shader->setUniform("model_matrix", transform.getLocalToWorldMatrix(), false);
	for (auto m : meshes)
	{
		m->draw(rid);
	}
}

void Model::drawWithMaterials(ShaderProgram* shader, size_t rid)
{
	shader->setUniform("model_matrix", transform.getLocalToWorldMatrix(), false);
	for (auto m : meshes)
	{
		m->drawWithMaterial(shader, rid);
	}
}
