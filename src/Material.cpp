#include "Material.h"

Material::Material() :
	m_diffuse_albedo(nullptr),
	m_normals(nullptr),
	m_specular_albedo(nullptr),
	m_roughness(nullptr),
	m_aniso_rotation(nullptr),
	m_fresnel_f0(nullptr),
	m_displacement(nullptr),
	m_transparency(nullptr),
	m_mid(0)
{
}

Material::Material(
	Texture* diffuse_albedo,
	Texture* normals,
	Texture* specular_albedo,
	Texture* roughness,
	Texture* aniso_rotation,
	Texture* fresnel_f0,
	Texture* displacement,
	Texture* transparency,
	std::size_t mid
	) :
	m_diffuse_albedo(diffuse_albedo),
	m_normals(normals),
	m_specular_albedo(specular_albedo),
	m_roughness(roughness),
	m_aniso_rotation(aniso_rotation),
	m_fresnel_f0(fresnel_f0),
	m_displacement(displacement),
	m_transparency(transparency),
	m_mid(mid)
{
}

void Material::bind(ShaderProgram * _shader)
{
	_shader->bindTex("material.s_diffuse_albedo", m_diffuse_albedo);
	_shader->bindTex("material.s_normals", m_normals);
	_shader->bindTex("material.s_specular_albedo", m_specular_albedo);
	_shader->bindTex("material.s_roughness", m_roughness);
	_shader->bindTex("material.s_aniso_rotation", m_aniso_rotation);
	_shader->bindTex("material.s_fresnel_f0", m_fresnel_f0);
	_shader->bindTex("material.s_displacement", m_displacement);
	_shader->bindTex("material.s_transparency", m_transparency);
}
