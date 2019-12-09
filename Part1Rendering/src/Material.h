/** \addtogroup resources
Material class
*  @{
*/

/*!
\file Material.h
*/

#ifndef _MATERIAL_H_
#define _MATERIAL_H_
#include <Texture.h>
#include <Shader.h>
#include <libheaders.h>

/*!
\brief Implements a material class which stores either textures or constant parameters for diffuse color, reflectivity, shininess and normals
*/
class Material
{
public:
	//! Some default material
	Material();
	
	Material(	Texture* diffuse_albedo,
				Texture* normals,
				Texture* specular_albedo,
				Texture* roughness,
				Texture* aniso_rotation,
				Texture* fresnel_f0,
				Texture* displacement,
				Texture* transparency,
				std::size_t mid = 0
	);

	//! Copy ctor
	Material(const Material& other) = default;
	//! Move ctor
	Material& operator=(const Material& other) = default;

	
	void bind(ShaderProgram* _shader, const glm::vec2& _mscale = glm::vec2(0.0, 0.0));

	Texture* m_diffuse_albedo;
	Texture* m_normals;
	Texture* m_specular_albedo;
	Texture* m_roughness;
	Texture* m_aniso_rotation;
	Texture* m_fresnel_f0;
	Texture* m_displacement;
	Texture* m_transparency;

	std::size_t m_mid;

	//! Renderer id
	std::size_t renderer_id;
};

#endif
/** @}*/