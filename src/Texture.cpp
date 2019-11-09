#include "Texture.h"

Texture::Texture() :
	tex(0),
	type(TexType::T2D),
	target(GL_TEXTURE_2D),
	size(0),
	miplevels(0)
{
}

Texture::~Texture()
{
	if (tex)
		glDeleteTextures(1, &tex);
}

Texture::Texture(Texture && other) :
	tex(other.tex),
	type(other.type),
	target(other.target),
	size(other.size),
	miplevels(other.miplevels),
	path(std::move(other.path))
{
	other.tex = 0;
	other.size = glm::ivec3(0);
	other.miplevels = 0;
}

Texture & Texture::operator=(Texture && other)
{
	if (this == &other)
		return *this;

	if (tex)
		glDeleteTextures(1, &tex);

	tex = other.tex;
	type = other.type;
	target = other.target;
	size = other.size;
	miplevels = other.miplevels;
	path = std::move(other.path);

	other.tex = 0;
	other.size = glm::ivec3(0);
	other.miplevels = 0;

	return *this;
}

bool Texture::bind(GLint tu)
{
	if (tex)
	{
		glActiveTexture(GL_TEXTURE0 + tu);
		glBindTexture(target, tex);
		return true;
	}
	return false;
}

void Texture::unbind()
{
	glBindTexture(target, 0);
}

void Texture::setTexParam(GLenum param, GLint value)
{
	glTexParameteri(target, param, value);
}

void Texture::setTexParam(GLenum param, GLfloat value)
{
	glTexParameterf(target, param, value);
}

void Texture::setTexParamArr(GLenum param, GLint * values, size_t n)
{
	glTexParameteriv(target, param, static_cast<GLint*>(values));
}

void Texture::setTexParamArr(GLenum param, GLfloat * values, size_t n)
{
	glTexParameterfv(target, param, static_cast<GLfloat*>(values));
}

std::unique_ptr<Texture> Texture::T2DFromData(GLenum internalformat, GLsizei width, GLsizei height, GLenum format, GLenum etype, void * data, bool genMipMaps)
{
	GLuint tex;
	glGenTextures(1, &tex);
	glBindTexture(GL_TEXTURE_2D, tex);
	glTexImage2D(GL_TEXTURE_2D, 0, internalformat, width, height, 0, format, etype, data);
	if (genMipMaps)
		glGenerateMipmap(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 0);
	return std::unique_ptr<Texture>(new Texture(tex, TexType::T2D, GL_TEXTURE_2D, { width, height, 1 }, 1 + static_cast<GLuint>(glm::floor(glm::log2(static_cast<float>(glm::max(width, height)))))));
}

std::unique_ptr<Texture> Texture::T2D(GLenum internalformat, GLsizei width, GLsizei height, GLsizei levels)
{
	GLuint tex;
	glGenTextures(1, &tex);
	glBindTexture(GL_TEXTURE_2D, tex);
	glTexStorage2D(GL_TEXTURE_2D, levels, internalformat, width, height);
	glBindTexture(GL_TEXTURE_2D, 0);
	return std::unique_ptr<Texture>(new Texture(tex, TexType::T2D, GL_TEXTURE_2D, { width, height, 1 }, levels));
}

std::unique_ptr<Texture> Texture::T2DArr(GLenum internalformat, GLsizei width, GLsizei height, GLsizei layers, GLsizei levels)
{
	GLuint tex;
	glGenTextures(1, &tex);
	glBindTexture(GL_TEXTURE_2D_ARRAY, tex);
	glTexStorage3D(GL_TEXTURE_2D_ARRAY, levels, internalformat, width, height, layers);
	glBindTexture(GL_TEXTURE_2D_ARRAY, 0);
	return std::unique_ptr<Texture>(new Texture(tex, TexType::T2DArr, GL_TEXTURE_2D_ARRAY, { width, height, layers }, levels));
}

std::unique_ptr<Texture> Texture::T3D(GLenum internalformat, GLsizei width, GLsizei height, GLsizei depth, GLsizei levels)
{
	GLuint tex;
	glGenTextures(1, &tex);
	glBindTexture(GL_TEXTURE_3D, tex);
	glTexStorage3D(GL_TEXTURE_3D, levels, internalformat, width, height, depth);
	glBindTexture(GL_TEXTURE_3D, 0);
	return std::unique_ptr<Texture>(new Texture(tex, TexType::T3D, GL_TEXTURE_3D, { width, height, depth }, levels));
}

std::unique_ptr<Texture> Texture::TCBFromData(GLenum internalformat, GLsizei width, GLsizei height, GLenum format, GLenum etype, void ** data, bool genMipMaps)
{
	GLuint tex;
	glGenTextures(1, &tex);
	glBindTexture(GL_TEXTURE_CUBE_MAP, tex);
	glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X, 0, internalformat, width, height, 0, format, etype, data[0]);
	glTexImage2D(GL_TEXTURE_CUBE_MAP_NEGATIVE_X, 0, internalformat, width, height, 0, format, etype, data[1]);
	glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_Y, 0, internalformat, width, height, 0, format, etype, data[2]);
	glTexImage2D(GL_TEXTURE_CUBE_MAP_NEGATIVE_Y, 0, internalformat, width, height, 0, format, etype, data[3]);
	glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_Z, 0, internalformat, width, height, 0, format, etype, data[4]);
	glTexImage2D(GL_TEXTURE_CUBE_MAP_NEGATIVE_Z, 0, internalformat, width, height, 0, format, etype, data[5]);
	if (genMipMaps)
		glGenerateMipmap(GL_TEXTURE_CUBE_MAP);
	glBindTexture(GL_TEXTURE_CUBE_MAP, 0);
	return std::unique_ptr<Texture>(new Texture(tex, TexType::TCB, GL_TEXTURE_CUBE_MAP, { width, height, 1 }, 1 + static_cast<GLuint>(glm::floor(glm::log2(static_cast<GLfloat>(glm::max(width, height)))))));
}

std::unique_ptr<Texture> Texture::TCB(GLenum internalformat, GLsizei width, GLsizei height, GLsizei levels)
{
	GLuint tex;
	glGenTextures(1, &tex);
	glBindTexture(GL_TEXTURE_CUBE_MAP, tex);
	glTexStorage2D(GL_TEXTURE_CUBE_MAP, levels, internalformat, width, height);
	glBindTexture(GL_TEXTURE_CUBE_MAP, 0);
	return std::unique_ptr<Texture>(new Texture(tex, TexType::TCB, GL_TEXTURE_CUBE_MAP, { width, height, 1 }, levels));
}

Texture::Texture(GLint _tex, TexType _type, GLenum _target, const glm::ivec3 & _size, GLuint _miplevels) :
	tex(_tex),
	type(_type),
	target(_target),
	size(_size),
	miplevels(_miplevels)
{
}
