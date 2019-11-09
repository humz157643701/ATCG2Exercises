#include "Model.h"

void Model::draw(ShaderProgram* shader)
{
	shader->setUniform("model_matrix", transform.getLocalToWorldMatrix(), false);
	for (auto m : meshes)
	{
		m->draw();
	}
}

void Model::drawWithMaterials(ShaderProgram* shader)
{
	shader->setUniform("model_matrix", transform.getLocalToWorldMatrix(), false);
	for (auto m : meshes)
	{
		m->drawWithMaterial(shader);
	}
}
