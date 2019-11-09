#include "Lights.h"
#include <string>

bool DirectionalLight::bind(ShaderProgram * shader, const glm::mat4& mat, const char* instancename)
{
	bool success = shader->setUniform((std::string(instancename) + ".color").c_str(), color) && success;
	success = shader->setUniform((std::string(instancename) + ".direction").c_str(), glm::vec3(mat * glm::vec4(transform.getDirection(), 0.0f))) && success;
	return success;
}

DirectionalLight::DirectionalLight(const Transform & _transform, const glm::vec3 & _color) :
	transform(_transform),
	color(_color)
{
}

DirectionalLight::DirectionalLight() :
	transform(),
	color(1.0f)
{
}

bool PointLight::bind(ShaderProgram * shader, const glm::mat4& mat, const char* instancename)
{
	bool success = shader->setUniform((std::string(instancename) + ".color").c_str(), color) && success;
	success = shader->setUniform((std::string(instancename) + ".position").c_str(), glm::vec3(mat * glm::vec4(transform.getWorldPosition(), 1.0f))) && success;
	return success;
}

PointLight::PointLight(const Transform & _transform, const glm::vec3 & _color) :
	transform(_transform),
	color(_color)
{
}

PointLight::PointLight() :
	transform(),
	color(1.0f)
{
}

bool AmbientLight::bind(ShaderProgram * shader, const char * instancename)
{
	bool success = shader->setUniform((std::string(instancename) + ".color").c_str(), color) && success;
	return success;
}

AmbientLight::AmbientLight(const glm::vec3 & _color) :
	color(_color)
{
}

AmbientLight::AmbientLight() :
	color(glm::vec3(1.0f))
{
}
