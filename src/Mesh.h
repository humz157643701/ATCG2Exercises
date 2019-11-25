/** \addtogroup resources
Mesh class
*  @{
*/

/*!
\file Mesh.h
*/

#ifndef _MESH_H_
#define _MESH_H_
#include <Material.h>
#include <libheaders.h>
#include <vector>
#include <Shader.h>
#include <memory>
#include <CommonTypes.h>
#include <stdexcept>

/*!
\brief Class that manages geometry data ready to be rendered with OpenGL
*/
class Mesh
{
public:
	//! dtor
	~Mesh();
	//! No copying
	Mesh(const Mesh& other) = delete;
	//! No copying
	Mesh& operator=(const Mesh& other) = delete;
	//! Move ctor
	Mesh(Mesh&& other);
	//! Move assignment operator
	Mesh& operator=(Mesh&& other);
	/**
	\brief Sets the mesh's material uniforms and draws the mesh
	\param shader Shader for which the uniforms shall be set
	*/
	void drawWithMaterial(ShaderProgram* shader, size_t rendererid = 0);
	/**
	\brief Draws the mesh without setting any material uniforms
	*/
	void draw(size_t rendererid = 0);

	void setMaterial(Material* m);

	//! Creates a new mesh from the given vertex and index data and configures the vertex attributes accorting to the last parameter
	static std::unique_ptr<Mesh> createMesh(const std::vector<Vertex>& vd, const std::vector<Index>& id, const std::vector<VertexAttribute>& att);
	//! Creates a new mesh from the given vertex and index data and configures the vertex attributes accorting to the last parameter. Additionally sets a material to be used for renderirng.
	static std::unique_ptr<Mesh> createMesh(const std::vector<Vertex>& vd, const std::vector<Index>& id, const std::vector<VertexAttribute>& att, Material* material, const glm::vec2& material_scale = glm::vec2(1.0f, 1.0f));
private:
	Mesh(GLuint _vbo, GLuint _ibo, GLuint _vao, GLsizei _numIndices, Material* _material, const glm::vec2& material_scale = glm::vec2(0.0f, 0.0f));
	Material* material;
	GLuint vbo;
	GLuint ibo;
	GLuint vao;
	GLsizei numIndices;
	glm::vec2 mscale;
};

#endif
/** @}*/