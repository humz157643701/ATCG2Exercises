#include <Shader.h>
#include <fstream>
#include <sstream>

ShaderProgram::ShaderProgram(GLuint program) :
	prog(program)
{

}

ShaderProgram::~ShaderProgram()
{
	if(prog)
		glDeleteProgram(prog);
}

ShaderProgram::ShaderProgram(ShaderProgram&& other) :
	prog(other.prog)
{
	other.prog = 0;
}

ShaderProgram& ShaderProgram::operator=(ShaderProgram&& other)
{
	if (this == &other)
		return *this;

	if (prog)
		glDeleteProgram(prog);

	prog = other.prog;
	other.prog = 0;

	return *this;
}

ShaderProgram ShaderProgram::createShaderProgram(const std::string& vspath, const std::string& fspath)
{
	GLuint vertexShader;
	GLuint fragmentShader;
	GLuint program;
	std::string vertexCode;
	std::string fragmentCode;
	std::ifstream vShaderFile;
	std::ifstream fShaderFile;
	vShaderFile.exceptions(std::ifstream::badbit);
	fShaderFile.exceptions(std::ifstream::badbit);
	try
	{
		vShaderFile.open(vspath);
		if (!vShaderFile.is_open())
			throw std::invalid_argument("Vertex shader file not found.");
		fShaderFile.open(fspath);
		if (!fShaderFile.is_open())
			throw std::invalid_argument("Fragment shader file not found.");
		std::stringstream vShaderStream, fShaderStream;
		vShaderStream << vShaderFile.rdbuf();
		fShaderStream << fShaderFile.rdbuf();
		vShaderFile.close();
		fShaderFile.close();
		vertexCode = vShaderStream.str();
		fragmentCode = fShaderStream.str();
	}
	catch (const std::exception& ex)
	{
		std::string errmsg;
		errmsg.append("Error: Shader files couldn't be read:\n");
		errmsg.append(ex.what());
		throw std::logic_error(errmsg.c_str());
	}
	const GLchar* vShaderCode = vertexCode.c_str();
	const GLchar* fShaderCode = fragmentCode.c_str();
	GLint success;
	GLchar infoLog[512];
	vertexShader = glCreateShader(GL_VERTEX_SHADER); 
		glShaderSource(vertexShader, 1, &vShaderCode, NULL); 
		glCompileShader(vertexShader); 
		glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success); 
		if (!success)
		{
			glGetShaderInfoLog(vertexShader, 512, NULL, infoLog); 
				std::string errmsg;
			errmsg.append("Compiler error in vertex shader '" + vspath + "':\n");
			errmsg.append(infoLog);
			glDeleteShader(vertexShader); 
				throw std::logic_error(errmsg.c_str());
		}
	fragmentShader = glCreateShader(GL_FRAGMENT_SHADER); 
		glShaderSource(fragmentShader, 1, &fShaderCode, NULL); 
		glCompileShader(fragmentShader); 
		glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success); 
		if (!success)
		{
			glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog); 
				std::string errmsg;
			errmsg.append("Compiler error in fragment shader '" + fspath + "':\n");
			errmsg.append(infoLog);
			glDeleteShader(vertexShader); 
				glDeleteShader(fragmentShader); 
				throw std::logic_error(errmsg.c_str());
		}
	program = glCreateProgram(); 
		glAttachShader(program, vertexShader); 
		glAttachShader(program, fragmentShader); 
		glLinkProgram(program); 
		glGetProgramiv(program, GL_LINK_STATUS, &success); 
		if (!success)
		{
			glGetProgramInfoLog(program, 512, NULL, infoLog); 
				std::string errmsg;
			errmsg.append("Linker error in program:\n");
			errmsg.append(infoLog);
			glDetachShader(program, vertexShader); 
				glDetachShader(program, fragmentShader); 
				glDeleteShader(vertexShader); 
				glDeleteShader(fragmentShader); 
				throw std::logic_error(errmsg.c_str());
		}
	glDetachShader(program, vertexShader); 
	glDetachShader(program, fragmentShader); 
	glDeleteShader(vertexShader); 
	glDeleteShader(fragmentShader); 
	return ShaderProgram(program);
}

ShaderProgram ShaderProgram::createShaderProgram(const std::string & vspath, const std::string & fspath, const std::string & gspath)
{
	GLuint vertexShader;
	GLuint fragmentShader;
	GLuint geometryShader;
	GLuint program;
	std::string vertexCode;
	std::string fragmentCode;
	std::string geometryCode;
	std::ifstream vShaderFile;
	std::ifstream fShaderFile;
	std::ifstream gShaderFile;
	vShaderFile.exceptions(std::ifstream::badbit);
	fShaderFile.exceptions(std::ifstream::badbit);
	gShaderFile.exceptions(std::ifstream::badbit);
	try
	{
		vShaderFile.open(vspath);
		if (!vShaderFile.is_open())
			throw std::invalid_argument("Vertex shader file not found.");
		fShaderFile.open(fspath);
		if (!fShaderFile.is_open())
			throw std::invalid_argument("Fragment shader file not found.");
		gShaderFile.open(gspath);
		if (!gShaderFile.is_open())
			throw std::invalid_argument("Geometry shader file not found.");
		std::stringstream vShaderStream, fShaderStream, gShaderStream;
		vShaderStream << vShaderFile.rdbuf();
		fShaderStream << fShaderFile.rdbuf();
		gShaderStream << gShaderFile.rdbuf();
		vShaderFile.close();
		fShaderFile.close();
		gShaderFile.close();
		vertexCode = vShaderStream.str();
		fragmentCode = fShaderStream.str();
		geometryCode = gShaderStream.str();
	}
	catch (const std::exception& ex)
	{
		std::string errmsg;
		errmsg.append("Error: Shader files couldn't be read:\n");
		errmsg.append(ex.what());
		throw std::logic_error(errmsg.c_str());
	}
	const GLchar* vShaderCode = vertexCode.c_str();
	const GLchar* fShaderCode = fragmentCode.c_str();
	const GLchar* gShaderCode = geometryCode.c_str();
	GLint success;
	GLchar infoLog[512];
	vertexShader = glCreateShader(GL_VERTEX_SHADER); 
		glShaderSource(vertexShader, 1, &vShaderCode, NULL); 
		glCompileShader(vertexShader); 
		glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success); 
		if (!success)
		{
			glGetShaderInfoLog(vertexShader, 512, NULL, infoLog); 
				std::string errmsg;
			errmsg.append("Compiler error in vertex shader:\n");
			errmsg.append(infoLog);
			glDeleteShader(vertexShader); 
				throw std::logic_error(errmsg.c_str());
		}
	fragmentShader = glCreateShader(GL_FRAGMENT_SHADER); 
		glShaderSource(fragmentShader, 1, &fShaderCode, NULL); 
		glCompileShader(fragmentShader); 
		glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success); 
		if (!success)
		{
			glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog); 
				std::string errmsg;
			errmsg.append("Compiler error in fragment shader:\n");
			errmsg.append(infoLog);
			glDeleteShader(vertexShader); 
				glDeleteShader(fragmentShader); 
				throw std::logic_error(errmsg.c_str());
		}
	geometryShader = glCreateShader(GL_GEOMETRY_SHADER); 
		glShaderSource(geometryShader, 1, &gShaderCode, NULL); 
		glCompileShader(geometryShader); 
		glGetShaderiv(geometryShader, GL_COMPILE_STATUS, &success); 
		if (!success)
		{
			glGetShaderInfoLog(geometryShader, 512, NULL, infoLog); 
				std::string errmsg;
			errmsg.append("Compiler error in geometry shader:\n");
			errmsg.append(infoLog);
			glDeleteShader(geometryShader); 
				glDeleteShader(vertexShader); 
				glDeleteShader(fragmentShader); 
				throw std::logic_error(errmsg.c_str());
		}
	program = glCreateProgram(); 
		glAttachShader(program, vertexShader); 
		glAttachShader(program, fragmentShader); 
		glAttachShader(program, geometryShader); 
		glLinkProgram(program); 
		glGetProgramiv(program, GL_LINK_STATUS, &success); 
		if (!success)
		{
			glGetProgramInfoLog(program, 512, NULL, infoLog); 
				std::string errmsg;
			errmsg.append("Linker error in program:\n");
			errmsg.append(infoLog);
			glDetachShader(program, vertexShader); 
				glDetachShader(program, fragmentShader); 
				glDetachShader(program, geometryShader); 
				glDeleteShader(vertexShader); 
				glDeleteShader(fragmentShader); 
				glDeleteShader(geometryShader); 
				throw std::logic_error(errmsg.c_str());
		}
	glDetachShader(program, vertexShader); 
		glDetachShader(program, fragmentShader); 
		glDetachShader(program, geometryShader); 
		glDeleteShader(vertexShader); 
		glDeleteShader(fragmentShader); 
		glDeleteShader(geometryShader); 
		return ShaderProgram(program);
}

ShaderProgram ShaderProgram::createComputeShaderProgram(const std::string & cspath, size_t localSizeX, size_t localSizeY, size_t localSizeZ)
{
	GLuint computeShader;
	GLuint program;
	std::string computeCode;
	std::ifstream cShaderFile;
	cShaderFile.exceptions(std::ifstream::badbit);
	try
	{
		cShaderFile.open(cspath);
		if (!cShaderFile.is_open())
			throw std::invalid_argument("Compute shader file not found.");
		std::stringstream cShaderStream;
		cShaderStream << cShaderFile.rdbuf();
		cShaderFile.close();
		computeCode = cShaderStream.str();
	}
	catch (const std::exception& ex)
	{
		std::string errmsg;
		errmsg.append("Error: Shader files couldn't be read:\n");
		errmsg.append(ex.what());
		throw std::logic_error(errmsg.c_str());
	}
	setPlaceholder(computeCode, "%LOCAL_SIZE_X%", std::to_string(localSizeX));
	setPlaceholder(computeCode, "%LOCAL_SIZE_Y%", std::to_string(localSizeY));
	setPlaceholder(computeCode, "%LOCAL_SIZE_Z%", std::to_string(localSizeZ));
	const GLchar* cShaderCode = computeCode.c_str();
	GLint success;
	GLchar infoLog[512];
	computeShader = glCreateShader(GL_COMPUTE_SHADER);
	glShaderSource(computeShader, 1, &cShaderCode, NULL);
	glCompileShader(computeShader);
	glGetShaderiv(computeShader, GL_COMPILE_STATUS, &success);
	if (!success)
	{
		glGetShaderInfoLog(computeShader, 512, NULL, infoLog);
		std::string errmsg;
		errmsg.append("Compiler error in compute shader '" + cspath + "':\n");
		errmsg.append(infoLog);
		glDeleteShader(computeShader);
		throw std::logic_error(errmsg.c_str());
	}
	program = glCreateProgram();
	glAttachShader(program, computeShader);
	glLinkProgram(program);
	glGetProgramiv(program, GL_LINK_STATUS, &success);
	if (!success)
	{
		glGetProgramInfoLog(program, 512, NULL, infoLog);
		std::string errmsg;
		errmsg.append("Linker error in program:\n");
		errmsg.append(infoLog);
		glDetachShader(program, computeShader);
		glDeleteShader(computeShader);
		throw std::logic_error(errmsg.c_str());
	}
	glDetachShader(program, computeShader);
	glDeleteShader(computeShader);
	return ShaderProgram(program);
}
