/** \addtogroup resources
Texture class
*  @{
*/

/*!
\file Texture.h
*/

#ifndef _TEXTURE_H_
#define _TEXTURE_H_

#include <libheaders.h>
#include <memory>
#include <string>

/*!
\brief Enum class listing the supported texture types
*/
enum class TexType
{
	//! 2D Texture
	T2D,
	//! 2D Array Texture
	T2DArr,
	//! 3D Texture
	T3D,
	//! Cube Texture
	TCB	
};

enum class TexFileType
{
	PNG,
	HDR,
	JPG
};

/*!
\brief Implements a wrapper to create and handle OpenGL textures of various formats
*/
class Texture
{
public:
	//! Creates an empty texture object
	Texture();
	// dtor
	~Texture();

	//! No copying
	Texture(const Texture&) = delete;
	//! No copying
	Texture& operator=(const Texture&) = delete;

	//! Move ctor
	Texture(Texture&& other);
	//! Move assignment operator
	Texture& operator=(Texture&& other);

	//! Binds the texture to the given texture id
	bool bind(GLint tu);
	//! Unbinds the texture
	void unbind();

	GLenum getTarget() const { return target; }

	//! Sets an OpenGL texture parameter for this texture
	void setTexParam(GLenum param, GLint value);
	//! Sets an OpenGL texture parameter for this texture
	void setTexParam(GLenum param, GLfloat value);
	//! Sets an array of OpenGL texture parameter for this texture
	void setTexParamArr(GLenum param, GLint* values, size_t n);
	//! Sets an array of OpenGL texture parameter for this texture
	void setTexParamArr(GLenum param, GLfloat* values, size_t n);

	//factory functions

	/**
	\brief Creates a 2D Texture with the given width, height and format and fills it with some texture data

	\param internalFormat OpenGl texture format
	\param width Desired width of the texture
	\param height Desired height of the texture
	\param format Format of the input data (as one would pass it to TexImage*D...)
	\param etype Basic data type, GL_FLOAT, GL_UNSIGNED_BYTE  and so on
	\param data Pointer so some raw texture data
	\param genMipMaps If true, mipmaps are generated after the texture is created
	*/
	static std::unique_ptr<Texture> T2DFromData(GLenum internalformat, GLsizei width, GLsizei height, GLenum format, GLenum etype, void* data, bool genMipMaps);


	/**
	\brief Creates a 2D Texture from a file

	\param internalFormat OpenGl texture format
	\param format Format of the input data (as one would pass it to TexImage*D...)
	\param etype Basic data type, GL_FLOAT, GL_UNSIGNED_BYTE  and so on
	\param data Pointer so some raw texture data
	\param genMipMaps If true, mipmaps are generated after the texture is created
	*/
	static std::unique_ptr<Texture> T2DFromFile(const std::string& path, TexFileType filetype, GLenum internalformat, GLenum format, bool genMipMaps);

	/**
	\brief Creates an empty 2D texture

	\param internalFormat OpenGl texture format
	\param width Desired width of the texture
	\param height Desired height of the texture
	\param levels How many mipmap levels to reserve memory for
	*/
	static std::unique_ptr<Texture> T2D(GLenum internalformat, GLsizei width, GLsizei height, GLsizei levels);
	/**
	\brief Creates an empty 2D array texture

	\param internalFormat OpenGl texture format
	\param width Desired width of the texture
	\param height Desired height of the texture
	\param layers Desired number of texture layers
	\param levels How many mipmap levels to reserve memory for
	*/
	static std::unique_ptr<Texture> T2DArr(GLenum internalformat, GLsizei width, GLsizei height, GLsizei layers, GLsizei levels);
	/**
	\brief Creates an empty 3D texture

	\param internalFormat OpenGl texture format
	\param width Desired width of the texture
	\param height Desired height of the texture
	\param depth Desired depth of the texture
	\param levels How many mipmap levels to reserve memory for
	*/
	static std::unique_ptr<Texture> T3D(GLenum internalformat, GLsizei width, GLsizei height, GLsizei depth, GLsizei levels);		

	/**
	\brief Creates a cube map and fills it with some texture data

	\param internalFormat OpenGl texture format
	\param width Desired width of the texture
	\param height Desired height of the texture
	\param layers Desired number of texture layers
	\param format Format of the input data (as one would pass it to TexImage*D...)
	\param etype Basic data type, GL_FLOAT, GL_UNSIGNED_BYTE  and so on
	\param data An array of 6 pointers pointing for some raw data for every face, in the order +X, -X, +Y, -Y, +Z, -Z. Beware: cube maps have their own logic when sampling them
	\param genMipMaps If true, mipmaps are generated after the texture is created
	*/
	static std::unique_ptr<Texture> TCBFromData(GLenum internalformat, GLsizei width, GLsizei height, GLenum format, GLenum etype, void** data, bool genMipMaps);

	/**
	\brief Creates an empty cube texture

	\param internalFormat OpenGl texture format
	\param width Desired width of the texture
	\param height Desired height of the texture
	\param levels How many mipmap levels to reserve memory for
	*/
	static std::unique_ptr<Texture> TCB(GLenum internalformat, GLsizei width, GLsizei height, GLsizei levels);

	static GLsizei calculateMipMapLevels(GLsizei width, GLsizei height);
	static GLsizei calculateMipMapLevels(const glm::ivec3& sz);

	void generateMipMap();

	//! Can be set to identify the texture later
	std::string path;

	GLuint tex;
	TexType type;
	GLenum target;
	glm::ivec3 size;
	GLuint miplevels;
private:
	Texture(GLint _tex, TexType _type, GLenum _target, const glm::ivec3& _size, GLuint _miplevels);	
};

#endif

/** @}*/