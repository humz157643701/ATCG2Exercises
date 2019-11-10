#ifndef _COMMON_TYPES_H_
#define _COMMON_TYPES_H_
#include <libheaders.h>

//! Typedef for index data
using Index = GLuint;

//! Used to conveniently specify OpenGL's vertex attributes
struct VertexAttribute
{
	//! Number of components
	GLint n;
	//! OpenGL data type. (GL_FLOAT, GL_INT ...)
	GLenum type;
	//! Stride in bytes
	GLsizei stride;
	//! Attribute offset in bytes
	GLintptr offset;
};

struct Vertex
{
	glm::vec3 position;
	glm::vec2 uv;
	glm::vec3 normal;
	glm::vec3 tangent;
};

#endif