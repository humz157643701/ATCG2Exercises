/** \addtogroup resources
Abstraction class for shader programs
*  @{
*/

/*!
\file Shader.h
*/

#ifndef _SHADER_H_
#define _SHADER_H_
#include <libheaders.h>
#include <Texture.h>
#include <string>

/*!
\brief Wraps OpenGL shader program creatiion and use into a more convenient C++ interface.
*/
class ShaderProgram
{
public:
	//! Create a ShaderProgram object that wraps the OpenGL Shader Program with the given id.
	ShaderProgram(GLuint program);
	//! Dtor
	~ShaderProgram();
	//! No copying
	ShaderProgram(const ShaderProgram& other) = delete;
	//! No copying
	ShaderProgram& operator=(const ShaderProgram& other) = delete;
	//! Move ctor
	ShaderProgram(ShaderProgram&& other);
	//! Move assignment operator
	ShaderProgram& operator=(ShaderProgram&& other);
	//! "Uses" the shader program, aka binds it to the pipeline
	void use();
	//! Program id
	GLuint prog;
	//! Used to store the current texture unit
	GLint currentTu;
	//! Used to store a previous texture unit
	GLint tuCheckpoint;

	//! Returns true if the program is currently bound to the pipeline
	bool isActive();
	//! Binds a texture, sets the sampler uniform called name and increases the current texture unit
	bool bindTex(const char* name, Texture* tex);
	//! Only sets the sampler uniform name to the current texture unit to prevent invalid operation errors
	bool occupyTex(const char* name);
	//! Resets the current texture unit to the given index
	void resetTU(int newCurrentTU = 0);
	//! Returns the current texture unit
	int getCurrentTU();
	//! Stores the current texture unit in the checkpoint member
	void saveTU();
	//! Restores the current texture unit to the index stored in checkpoint
	void restoreTU();

	//! Sets a uniform, returns true on success
	bool setUniform(const char* name, GLfloat value);
	//! Sets a uniform, returns true on success
	bool setUniform(const char* name, const glm::vec2& value);
	//! Sets a uniform, returns true on success
	bool setUniform(const char* name, const glm::vec3& value);
	//! Sets a uniform, returns true on success
	bool setUniform(const char* name, const glm::vec4& value);

	//! Sets a uniform, returns true on success
	bool setUniform(const char* name, GLint value);
	//! Sets a uniform, returns true on success
	bool setUniform(const char* name, const glm::ivec2& value);
	//! Sets a uniform, returns true on success
	bool setUniform(const char* name, const glm::ivec3& value);
	//! Sets a uniform, returns true on success
	bool setUniform(const char* name, const glm::ivec4& value);

	//! Sets a uniform, returns true on success
	bool setUniform(const char* name, GLuint value);
	//! Sets a uniform, returns true on success
	bool setUniform(const char* name, const glm::uvec2& value);
	//! Sets a uniform, returns true on success
	bool setUniform(const char* name, const glm::uvec3& value);
	//! Sets a uniform, returns true on success
	bool setUniform(const char* name, const glm::uvec4& value);

	//! Sets a uniform, returns true on success
	bool setUniform(const char* name, const glm::mat2& value, bool transpose);
	//! Sets a uniform, returns true on success
	bool setUniform(const char* name, const glm::mat3& value, bool transpose);
	//! Sets a uniform, returns true on success
	bool setUniform(const char* name, const glm::mat4& value, bool transpose);

	/**
	\brief Creates a shader program consisting of a vertex and a fragment shader.

	\param vspath Path to the vertex shader file
	\param fspath Path to the fragment shader file
	*/
	static ShaderProgram createShaderProgram(const std::string& vspath, const std::string& fspath);
	/**
	\brief Creates a shader program consisting of a vertex, fragment and a geometry shader.

	\param vspath Path to the vertex shader file
	\param fspath Path to the fragment shader file
	\param gspath Path to the geometry shader file
	*/
	static ShaderProgram createShaderProgram(const std::string& vspath, const std::string& fspath, const std::string& gspath);
	/**
	\brief Creates a compute shader program

	The strings %LOCAL_SIZE_X%, %LOCAL_SIZE_Y% and %LOCAL_SIZE_Z% in the shader file are replaced woth the respective size parameter given.
	Can be used to configure local work group sizes more conveniently.

	\param cspath Path to the compute shader file
	\param localSizeX Work group size in dimension x
	\param localSizeY Work group size in dimension y
	\param localSizeZ Work group size in dimension z
	*/
	static ShaderProgram createComputeShaderProgram(const std::string& cspath, size_t localSizeX = 1, size_t localSizeY = 1, size_t localSizeZ = 1);

private:
	GLint getUniformLocation(const char* name);
	static void setPlaceholder(std::string& templ, const std::string& placeholder, const std::string& value);
};

inline void ShaderProgram::saveTU()
{
	tuCheckpoint = currentTu;
}

inline void ShaderProgram::restoreTU()
{
	currentTu = tuCheckpoint;
}

inline void ShaderProgram::use()
{	
	currentTu = 0;
	tuCheckpoint = 0;
	GLint current;
	glGetIntegerv(GL_CURRENT_PROGRAM, &current);
		if (current != prog && prog != 0)
			glUseProgram(prog);
}

inline GLint ShaderProgram::getUniformLocation(const char* name)
{
	return glGetUniformLocation(this->prog, name);
}

inline void ShaderProgram::setPlaceholder(std::string & templ, const std::string & placeholder, const std::string& value)
{
	size_t f = templ.find(placeholder, 0);
	while (f != templ.npos)
	{
		templ.replace(f, placeholder.length(), value);
		f = templ.find(placeholder, f + placeholder.length());
	}
}

inline bool ShaderProgram::isActive()
{
	GLint progName = 0;
	glGetIntegerv(GL_CURRENT_PROGRAM, &progName);
	if (progName != this->prog)
		return false;
	return true;
}

inline bool ShaderProgram::occupyTex(const char * name)
{
	if (isActive() && getUniformLocation(name) != -1)
	{
		setUniform(name, currentTu);
		++currentTu;
		return true;
	}
	return false;
}

inline bool ShaderProgram::bindTex(const char * name, Texture * tex)
{
	if (isActive() && getUniformLocation(name) != -1)
	{
		if (tex)
			tex->bind(currentTu);
		setUniform(name, currentTu);
		++currentTu;
		return true;
	}
	return false;
}

inline void ShaderProgram::resetTU(int newCurrentTU)
{
	currentTu = newCurrentTU;
}

inline int ShaderProgram::getCurrentTU()
{
	return currentTu;
}



inline bool ShaderProgram::setUniform(const char* name, GLfloat value)
{
	GLint loc = getUniformLocation(name);
	if (loc == -1)
		return false;
	if (!isActive())
		return false;
	glUniform1f(loc, value); 
		return true;
}

inline bool ShaderProgram::setUniform(const char* name, const glm::vec2& value)
{
	GLint loc = getUniformLocation(name);
	if (loc == -1)
		return false;
	if (!isActive())
		return false;
	glUniform2fv(loc, 1, glm::value_ptr(value)); 
		return true;
}

inline bool ShaderProgram::setUniform(const char* name, const glm::vec3& value)
{
	GLint loc = getUniformLocation(name);
	if (loc == -1)
		return false;
	if (!isActive())
		return false;
	glUniform3fv(loc, 1, glm::value_ptr(value)); 
		return true;
}

inline bool ShaderProgram::setUniform(const char* name, const glm::vec4& value)
{
	GLint loc = getUniformLocation(name);
	if (loc == -1)
		return false;
	if (!isActive())
		return false;
	glUniform4fv(loc, 1, glm::value_ptr(value)); 
		return true;
}

inline bool ShaderProgram::setUniform(const char* name, GLint value)
{
	GLint loc = getUniformLocation(name);
	if (loc == -1)
		return false;
	if (!isActive())
		return false;
	glUniform1i(loc, value); 
		return true;
}

inline bool ShaderProgram::setUniform(const char* name, const glm::ivec2& value)
{
	GLint loc = getUniformLocation(name);
	if (loc == -1)
		return false;
	if (!isActive())
		return false;
	glUniform2iv(loc, 1, glm::value_ptr(value)); 
		return true;
}

inline bool ShaderProgram::setUniform(const char* name, const glm::ivec3& value)
{
	GLint loc = getUniformLocation(name);
	if (loc == -1)
		return false;
	if (!isActive())
		return false;
	glUniform3iv(loc, 1, glm::value_ptr(value)); 
		return true;
}

inline bool ShaderProgram::setUniform(const char* name, const glm::ivec4& value)
{
	GLint loc = getUniformLocation(name);
	if (loc == -1)
		return false;
	if (!isActive())
		return false;
	glUniform4iv(loc, 1, glm::value_ptr(value)); 
		return true;
}

inline bool ShaderProgram::setUniform(const char* name, GLuint value)
{
	GLint loc = getUniformLocation(name);
	if (loc == -1)
		return false;
	if (!isActive())
		return false;
	glUniform1ui(loc, value); 
		return true;
}

inline bool ShaderProgram::setUniform(const char* name, const glm::uvec2& value)
{
	GLint loc = getUniformLocation(name);
	if (loc == -1)
		return false;
	if (!isActive())
		return false;
	glUniform2uiv(loc, 1, glm::value_ptr(value)); 
		return true;
}

inline bool ShaderProgram::setUniform(const char* name, const glm::uvec3& value)
{
	GLint loc = getUniformLocation(name);
	if (loc == -1)
		return false;
	if (!isActive())
		return false;
	glUniform3uiv(loc, 1, glm::value_ptr(value)); 
		return true;
}

inline bool ShaderProgram::setUniform(const char* name, const glm::uvec4& value)
{
	GLint loc = getUniformLocation(name);
	if (loc == -1)
		return false;
	if (!isActive())
		return false;
	glUniform4uiv(loc, 1, glm::value_ptr(value)); 
		return true;
}

inline bool ShaderProgram::setUniform(const char*name, const glm::mat2& value, bool transpose)
{
	GLint loc = getUniformLocation(name);
	if (loc == -1)
		return false;
	if (!isActive())
		return false;
	glUniformMatrix2fv(loc, 1, transpose ? GL_TRUE : GL_FALSE, glm::value_ptr(value)); 
		return true;
}

inline bool ShaderProgram::setUniform(const char* name, const glm::mat3& value, bool transpose)
{
	GLint loc = getUniformLocation(name);
	if (loc == -1)
		return false;
	if (!isActive())
		return false;
	glUniformMatrix3fv(loc, 1, transpose ? GL_TRUE : GL_FALSE, glm::value_ptr(value)); 
		return true;
}

inline bool ShaderProgram::setUniform(const char* name, const glm::mat4& value, bool transpose)
{
	GLint loc = getUniformLocation(name);
	if (loc == -1)
		return false;
	if (!isActive())
		return false;
	glUniformMatrix4fv(loc, 1, transpose ? GL_TRUE : GL_FALSE, glm::value_ptr(value)); 
		return true;
}

#endif

/** @}*/