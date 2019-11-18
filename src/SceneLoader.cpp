#include "SceneLoader.h"
#include <sstream>
#include <cctype>
#include <fstream>
#include <iostream>
#include <stack>
#include <queue>
#include <stb_image.h>
#include <OBJLoader.h>
#include <AxfReader.hpp>
#include <cstddef>

//internal types --------------------------------------------------------------------------------
enum class TokType
{
	String,
	Name,
	Float,
	Int,
	Invalid,
	LineBreak,
	EndOfFile,
	Boolean
};

struct Token
{
	Token() :
		type(TokType::Invalid),
		value(""),
		intvalue(0)
	{
	}

	Token(TokType _type, const std::string& sval) :
		type(_type),
		value(sval),
		intvalue(0)
	{}

	Token(TokType _type, const std::string& sval, int64_t intval) :
		type(_type),
		value(sval),
		intvalue(intval)
	{
	}

	Token(TokType _type, const std::string& sval, double floatval) :
		type(_type),
		value(sval),
		floatvalue(floatval)
	{
	}

	Token(TokType _type, const std::string& sval, bool boolval) :
		type(_type),
		value(sval),
		boolvalue(boolval)
	{}

	TokType type;
	std::string value;
	union
	{
		bool boolvalue;
		int64_t intvalue;
		double floatvalue;
	};
};

class Scanner
{
public:
	std::string file;
	size_t pos;
	size_t filesize;
	size_t row;
	size_t col;
	Token lh;

	Scanner(std::istream& _stream)
	{
		std::stringstream ss;
		ss << _stream.rdbuf();
		file = ss.str();
		filesize = file.length();
		pos = 0;
		row = 0;
		col = 0;
		lh = Token{ TokType::Invalid, "" };
		next();
	}

	std::string getErrorString(const std::string& msg)
	{
		std::string estr = "Config file parsing error at line: ";
		estr += std::to_string(row) + ", column: ";
		estr += std::to_string(col + 1) + "\n";
		estr += msg;
		return estr;
	}

	bool eof()
	{
		return pos >= filesize;
	}

	const Token& next()
	{
		if (lh.type == TokType::EndOfFile)
			return lh;

		if (eof())
		{
			lh = Token(TokType::EndOfFile, "<eof>");
			return lh;
		}

		
		while (!eof() && (file[pos] == ' ' || file[pos] == '\t'))
		{
			++pos;
		}
		

		//comments
		if (file[pos] == '#')
		{
			while (!eof() && file[pos] != '\n') ++pos;
		}

		//linebreaks
		if (file[pos] == '\n' ||
			((pos + 1 < filesize) &&
				file[pos] == '\r' &&
				file[pos + 1] == '\n'))
		{
			col = 0;
			//position
			if (file[pos] == '\r' && file[pos + 1] == '\n')
				pos += 2;
			else
				pos++;
			++row;
			lh = Token(TokType::LineBreak, "<linebreak>");
			return lh;
		}

		if (file[pos] == '.' || file[pos] == '-' || std::isdigit(static_cast<unsigned char>(file[pos])))
		{
			//test for float or integer, maybe include e notation
			std::string num("");
			std::string epart("");
			int pot = 0;
			size_t curpos = pos;
			bool dot = false;
			bool sign = false;
			if (file[pos] == '-')
			{
				sign = true;
				++curpos;
				++col;
			}
			while (!eof() && (file[curpos] == '.' || std::isdigit(static_cast<unsigned char>(file[curpos]))))
			{
				if (std::isdigit(static_cast<unsigned char>(file[curpos])))
				{
					num += file[curpos];
					if (dot)
						pot--;
					curpos++;
					col++;
					continue;
				}

				if (file[curpos] == '.' && !dot)
				{
					dot = true;
					num += '.';
					curpos++;
					col++;
					continue;
				}

				if (file[curpos] == '.' && dot)
				{
					throw std::logic_error(getErrorString("Wrong floating point format"));
				}
			}

			pos = curpos;

			if (!eof() && file[pos] == 'e')
			{
				++pos;
				++col;
				epart += 'e';
				bool esign = false;
				if (!eof() && file[pos] == '-')
				{
					esign = true;
					epart += '-';
					++pos;
					++col;
				}
				int epot = 0;
				while (!eof() && std::isdigit(static_cast<unsigned char>(file[pos])))
				{
					epart += file[pos];
					switch (file[pos])
					{
					case '0':
						epot = epot * 10 + 0;
						break;
					case '1':
						epot = epot * 10 + 1;
						break;
					case '2':
						epot = epot * 10 + 2;
						break;
					case '3':
						epot = epot * 10 + 3;
						break;
					case '4':
						epot = epot * 10 + 4;
						break;
					case '5':
						epot = epot * 10 + 5;
						break;
					case '6':
						epot = epot * 10 + 6;
						break;
					case '7':
						epot = epot * 10 + 7;
						break;
					case '8':
						epot = epot * 10 + 8;
						break;
					case '9':
						epot = epot * 10 + 9;
					}
					++pos;
					++col;
				}
				epot = (esign ? -epot : epot);
				pot += epot;
			}

			int64_t rawnum = 0;
			for (size_t i = 0; i < num.length(); i++)
			{
				switch (num[i])
				{
				case '0':
					rawnum = rawnum * 10 + 0;
					break;
				case '1':
					rawnum = rawnum * 10 + 1;
					break;
				case '2':
					rawnum = rawnum * 10 + 2;
					break;
				case '3':
					rawnum = rawnum * 10 + 3;
					break;
				case '4':
					rawnum = rawnum * 10 + 4;
					break;
				case '5':
					rawnum = rawnum * 10 + 5;
					break;
				case '6':
					rawnum = rawnum * 10 + 6;
					break;
				case '7':
					rawnum = rawnum * 10 + 7;
					break;
				case '8':
					rawnum = rawnum * 10 + 8;
					break;
				case '9':
					rawnum = rawnum * 10 + 9;
				case '.':
					break;
				}
			}

			if (dot && num.length() >= 2) //float
			{
				double res = (sign ? -1.0 : 1.0) * static_cast<double>(rawnum) * std::pow(10.0, static_cast<double>(pot));
				lh = Token(TokType::Float, num + epart, res);
				return lh;
			}
			else if (!dot && num.length() >= 1)
			{
				int64_t res = (sign ? -1 : 1) * rawnum * static_cast<int64_t>(std::pow(int64_t(10), static_cast<int64_t>(pot)));
				lh = Token(TokType::Int, num + epart, res);
				return lh;
			}
			else
			{
				throw std::logic_error(getErrorString("Wrong number format"));
			}
		}

		if (file[pos] == '"')
		{
			std::string str = "";
			std::string::value_type c;
			++pos;
			++col;
			bool stringclosed = false;
			while (!eof())
			{
				c = file[pos++];
				++col;
				if (c != '"')
					str += c;
				else
				{
					stringclosed = true;
					break;
				}
			}
			if (!stringclosed)
			{
				throw std::logic_error(getErrorString("Unclosed string literal"));
			}
			lh = Token(TokType::String, str);
			return lh;
		}

		//name or true/false
		if (std::isalpha(static_cast<unsigned char>(file[pos])) || file[pos] == '_')
		{
			std::string idname = "";
			idname += file[pos++];
			++col;
			while (!eof() && (std::isalnum(static_cast<unsigned char>(file[pos])) || file[pos] == '_'))
			{
				idname += file[pos++];
				++col;
			}
			if (idname == "true")
			{
				lh = Token(TokType::Boolean, idname, true);
			}
			else if (idname == "false")
			{
				lh = Token(TokType::Boolean, idname, false);
			}
			else
			{
				lh = Token(TokType::Name, idname);
			}
			return lh;
		}

		std::string estr = "Unexpected '" + lh.value + "'";
		throw std::logic_error(getErrorString(estr));
	}

	Scanner& match(TokType mt)
	{
		if (lh.type == mt)
		{
			next();
			return *this;
		}
		std::string estr = "Unexpected '" + lh.value + "'";
		throw std::logic_error(getErrorString(estr));
	}

	const Token& lookahead()
	{
		return lh;
	}
};

std::unique_ptr<Scene> SceneLoader::loadScene(const std::string & scenepath)
{
	Scene* scn = new Scene();
	try
	{
		std::ifstream fs(scenepath, std::ios::in | std::ios::binary);
		if (!fs.is_open())
			throw std::logic_error("Scene file '" + scenepath + "' couldn't be read.");
		Scanner scan(fs);

		//create default material
		scn->clearColor = { 0.0f, 0.0f, 0.0f, 1.0f };
		
		parseFile(&scan, scn);		

		return std::unique_ptr<Scene>(scn);
	}
	catch (const std::exception& ex)
	{
		delete scn;
		throw ex;
	}
}

void SceneLoader::parseFile(Scanner * scan, Scene * scn)
{
	while (scan->lookahead().type != TokType::EndOfFile)
	{
		if (scan->lookahead().value == "CAMERA")
		{
			scan->match(TokType::Name);
			parseCamera(scan, scn);
		}		
		else if (scan->lookahead().value == "MODEL_OPAQUE")
		{
			scan->match(TokType::Name);
			parseModelOpaque(scan, scn);
		}
		else if (scan->lookahead().value == "MODEL_TRANSPARENT")
		{
			scan->match(TokType::Name);
			parseModelTransparent(scan, scn);
		}
		else if (scan->lookahead().value == "DIRECTIONAL_LIGHT")
		{
			scan->match(TokType::Name);
			parseDirLight(scan, scn);
		}
		else if (scan->lookahead().value == "POINT_LIGHT")
		{
			scan->match(TokType::Name);
			parsePointLight(scan, scn);
		}
		else if (scan->lookahead().value == "AMBIENT_LIGHT")
		{
			scan->match(TokType::Name);
			parseAmbientLight(scan, scn);
		}
		else if (scan->lookahead().value == "ENVIRONMENT_MAP")
		{
			scan->match(TokType::Name);
			parseEnvMap(scan, scn);
		}
		else if (scan->lookahead().value == "MATERIAL")
		{
			scan->match(TokType::Name);
			parseMaterial(scan, scn);
		}
		else if (scan->lookahead().value == "TM_EXPOSURE")
		{
			scan->match(TokType::Name);
			parseExposure(scan, scn);
		}
		else if (scan->lookahead().value == "SKYBOX_RES")
		{
			scan->match(TokType::Name);
			parseSkyboxRes(scan, scn);
		}
		else if (scan->lookahead().value == "CLEAR_COLOR")
		{
			scan->match(TokType::Name);
			parseClearColor(scan, scn);
		}
		else if (scan->lookahead().type == TokType::LineBreak)
		{
			scan->match(TokType::LineBreak);
		}
		else
		{
			throw std::logic_error("Unexpected sympbol: " + scan->lookahead().value);
		}
	}
}

void SceneLoader::parseCamera(Scanner * scan, Scene * scn)
{
	glm::vec3 p, l, u;
	float fov, near, far;
	std::size_t id;
	std::size_t pid;
	// object id
	id = static_cast<std::size_t>(scan->lookahead().intvalue); scan->match(TokType::Int);
	// parentid. no parent => 0
	pid = static_cast<std::size_t>(scan->lookahead().intvalue); scan->match(TokType::Int);

	p.x = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	p.y = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	p.z = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	l.x = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	l.y = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	l.z = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	u.x = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	u.y = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	u.z = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	fov = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	near = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	far = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);	

	glm::mat4 viewmat = glm::lookAt(p, l, u);

	scn->m_camera.transform = Transform(glm::inverse(viewmat), id);
	scn->m_camera.setParameters(800, 600, glm::radians(fov), near, far);	
	if(pid != 0)
		scn->m_camera.transform.setParent(scn->getTransformByOID(pid));
}

void SceneLoader::parseModelOpaque(Scanner * scan, Scene * scn)
{	
	std::string path;
	glm::vec3 pos;
	glm::vec3 scale;
	glm::quat rot;
	std::size_t id;
	std::size_t pid;
	std::size_t mid;
	// object id
	id = static_cast<std::size_t>(scan->lookahead().intvalue); scan->match(TokType::Int);
	// parentid. no parent => 0
	pid = static_cast<std::size_t>(scan->lookahead().intvalue); scan->match(TokType::Int);
	// material id
	mid = static_cast<std::size_t>(scan->lookahead().intvalue); scan->match(TokType::Int);
	// material scale x
	float sx = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	// material scale y
	float sy = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	// obj path
	path = scan->lookahead().value; scan->match(TokType::String);
	// transform
	pos.x =		static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	pos.y =		static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	pos.z =		static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	rot.w =		static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	rot.x =		static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	rot.y =		static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	rot.z =		static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	scale.x =	static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	scale.y =	static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	scale.z =	static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);

	OBJResult res = OBJLoader::loadOBJ(path, false, true);

	// for this simple application we simply merge all objects/meshes from one file into one model entity.
	Model model;
	Material* mat = scn->getMaterialByID(mid);
	if (mat == nullptr)
		throw std::logic_error("SCENE LOADER: Material for model couldn't be found.");
	for(OBJObject& o : res.objects)
	{
		for(OBJMesh& m : o.meshes)
		{
			if(!m.hasNormals)
				OBJLoader::recalculateNormals(m);
			scn->m_meshes.push_back(std::move(Mesh::createMesh(m.vertices, m.indices, m.atts, mat, glm::vec2(sx, sy))));
			model.meshes.push_back(scn->m_meshes.back().get());
		}
	}

	// set transform
	model.transform = Transform(pos, rot, scale, id);
	if(pid != 0)
		model.transform.setParent(scn->getTransformByOID(pid), false);

	// add model to scene
	scn->m_opaque_models.push_back(std::move(model));
}

void SceneLoader::parseModelTransparent(Scanner * scan, Scene * scn)
{
	std::string path;
	glm::vec3 pos;
	glm::vec3 scale;
	glm::quat rot;
	std::size_t id;
	std::size_t pid;
	std::size_t mid;
	// object id
	id = static_cast<std::size_t>(scan->lookahead().intvalue); scan->match(TokType::Int);
	// parentid. no parent => 0
	pid = static_cast<std::size_t>(scan->lookahead().intvalue); scan->match(TokType::Int);
	// material id
	mid = static_cast<std::size_t>(scan->lookahead().intvalue); scan->match(TokType::Int);
	// material scale x
	float sx = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	// material scale y
	float sy = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	// obj path
	path = scan->lookahead().value; scan->match(TokType::String);
	// transform
	pos.x =		static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	pos.y =		static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	pos.z =		static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	rot.w =		static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	rot.x =		static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	rot.y =		static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	rot.z =		static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	scale.x =	static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	scale.y =	static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	scale.z =	static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);

	OBJResult res = OBJLoader::loadOBJ(path, false, true);

	// for this simple application we simply merge all objects/meshes from one file into one model entity.
	Model model;
	Material* mat = scn->getMaterialByID(mid);
	if (mat == nullptr)
		throw std::logic_error("SCENE LOADER: Material for model couldn't be found.");
	for(OBJObject& o : res.objects)
	{
		for(OBJMesh& m : o.meshes)
		{
			if(!m.hasNormals)
				OBJLoader::recalculateNormals(m);
			scn->m_meshes.push_back(std::move(Mesh::createMesh(m.vertices, m.indices, m.atts, mat, glm::vec2(sx, sy))));
			model.meshes.push_back(scn->m_meshes.back().get());
		}
	}

	// set transform
	model.transform = Transform(pos, rot, scale, id);
	if(pid != 0)
		model.transform.setParent(scn->getTransformByOID(pid), false);

	// add model to scene
	scn->m_transparent_models.push_back(std::move(model));
}

void SceneLoader::parseDirLight(Scanner * scan, Scene * scn)
{
	glm::vec3 c;
	glm::vec3 d;

	c.x =	static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	c.y =	static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	c.z =	static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	d.x =	static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	d.y =	static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	d.z =	static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	

	glm::mat3 rot;
	rot[2] = glm::normalize(-d);
	rot[0] = glm::normalize(glm::cross(glm::vec3(0.0f, 1.0f, 0.0f), rot[2]));
	rot[1] = glm::normalize(glm::cross(rot[2], rot[0]));

	scn->m_dirlights.push_back(DirectionalLight(Transform(glm::vec3(0.0f), glm::quat_cast(rot), glm::vec3(1.0f)), c));
}

void SceneLoader::parsePointLight(Scanner * scan, Scene * scn)
{
	glm::vec3 c;
	glm::vec3 p;
	std::size_t id;
	std::size_t pid;
	// object id
	id = static_cast<std::size_t>(scan->lookahead().intvalue); scan->match(TokType::Int);
	// parentid. no parent => 0
	pid = static_cast<std::size_t>(scan->lookahead().intvalue); scan->match(TokType::Int);
	c.x = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	c.y = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	c.z = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	p.x = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	p.y = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	p.z = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	
	Transform t = Transform(p, glm::quat(), glm::vec3(1.0f), id);
	if(pid != 0)
		t.setParent(scn->getTransformByOID(pid));
	scn->m_pointlights.push_back(PointLight(t, c));
}

void SceneLoader::parseAmbientLight(Scanner * scan, Scene * scn)
{
	glm::vec3 color;
	color.x = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	color.y = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	color.z = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	scn->m_ambientlights.push_back(AmbientLight(color));
}

void SceneLoader::parseExposure(Scanner * scan, Scene * scn)
{
	float exposure = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	scn->tm_exposure = exposure;
}

void SceneLoader::parseSkyboxRes(Scanner * scan, Scene * scn)
{
	GLsizei res = static_cast<GLsizei>(scan->lookahead().intvalue); scan->match(TokType::Int);
	scn->m_skyboxres = res;
}

void SceneLoader::parseClearColor(Scanner * scan, Scene * scn)
{
	glm::vec4 color;
	color.x = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	color.y = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	color.z = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	color.w = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	scn->clearColor = std::move(color);
}

void SceneLoader::parseEnvMap(Scanner * scan, Scene * scn)
{
	std::string path;
	bool isEquirectangularMap;

	path = scan->lookahead().value; scan->match(TokType::String);
	isEquirectangularMap = scan->lookahead().boolvalue; scan->match(TokType::Boolean);
	
	if (isEquirectangularMap)
	{
		std::unique_ptr<Texture> tex = Texture::T2DFromFile(path, TexFileType::HDR, GL_RGB32F, GL_RGB, true);

		float max_aniso;
		glGetFloatv(GL_MAX_TEXTURE_MAX_ANISOTROPY, &max_aniso);

		tex->bind(0);
		tex->setTexParam(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		tex->setTexParam(GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		tex->setTexParam(GL_TEXTURE_WRAP_S, GL_REPEAT);
		tex->setTexParam(GL_TEXTURE_WRAP_T, GL_REPEAT);
		tex->setTexParamArr(GL_TEXTURE_MAX_ANISOTROPY, &max_aniso, 1);
		tex->unbind();

		scn->m_er_env_map = std::move(tex);
	}
	else
	{
		throw std::logic_error("Cube maps are not yet supported by the scene loader.");
	}	
}

void SceneLoader::parseMaterial(Scanner * scan, Scene * scn)
{
	std::size_t mid;
	std::string path;

	float max_aniso;
	glGetFloatv(GL_MAX_TEXTURE_MAX_ANISOTROPY, &max_aniso);

	mid = static_cast<std::size_t>(scan->lookahead().intvalue); scan->match(TokType::Int);
	path = scan->lookahead().value; scan->match(TokType::String);

	AxfReader::TextureType textures;
	AxfReader::TextureDimType texture_dims;

	if (!AxfReader::readAxF(path, textures, texture_dims))
	{
		throw std::logic_error(("SCENE_LOADER: Error reading AxF " + path).c_str());
	}
	int dummy = 0;

	Material mat;
	GLenum internalformat;
	GLenum format;

	mat.m_mid = mid;

	if (textures.find("diffuse") == textures.end())
		throw std::logic_error("SCENE_LOADER: No diffuse texture found in .axf file.");
	switch (texture_dims["diffuse"].num_channels)
	{
	case 1:
		internalformat = GL_R32F;
		format = GL_RED;
		break;
	case 2:
		internalformat = GL_RG32F;
		format = GL_RG;
		break;
	case 3:
		internalformat = GL_RGB32F;
		format = GL_RGB;
		break;
	case 4:
		internalformat = GL_RGBA32F;
		format = GL_RGBA;
	default:
		throw std::logic_error("SCENE_LOADER: Invalid channel count for diffuse albedo.");
		break;
	}
	std::unique_ptr<Texture> diff_albedo = Texture::T2DFromData(	
		internalformat,
		texture_dims["diffuse"].width,
		texture_dims["diffuse"].height,
		format,
		GL_FLOAT,
		textures["diffuse"].data(),
		true
	);

	diff_albedo->bind(0);
	diff_albedo->setTexParam(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	diff_albedo->setTexParam(GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	diff_albedo->setTexParam(GL_TEXTURE_WRAP_S, GL_REPEAT);
	diff_albedo->setTexParam(GL_TEXTURE_WRAP_T, GL_REPEAT);
	diff_albedo->setTexParamArr(GL_TEXTURE_MAX_ANISOTROPY, &max_aniso, 1);
	diff_albedo->unbind();

	scn->m_textures.push_back(std::move(diff_albedo));
	mat.m_diffuse_albedo = scn->m_textures.back().get();

	if (textures.find("specular") == textures.end())
		throw std::logic_error("SCENE_LOADER: No specular texture found in .axf file.");
	switch (texture_dims["specular"].num_channels)
	{
	case 1:
		internalformat = GL_R32F;
		format = GL_RED;
		break;
	case 2:
		internalformat = GL_RG32F;
		format = GL_RG;
		break;
	case 3:
		internalformat = GL_RGB32F;
		format = GL_RGB;
		break;
	case 4:
		internalformat = GL_RGBA32F;
		format = GL_RGBA;
	default:
		throw std::logic_error("SCENE_LOADER: Invalid channel count for specular albedo.");
		break;
	}
	std::unique_ptr<Texture> specular_albedo = Texture::T2DFromData(
		internalformat,
		texture_dims["specular"].width,
		texture_dims["specular"].height,
		format,
		GL_FLOAT,
		textures["specular"].data(),
		true
	);

	specular_albedo->bind(0);
	specular_albedo->setTexParam(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	specular_albedo->setTexParam(GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	specular_albedo->setTexParam(GL_TEXTURE_WRAP_S, GL_REPEAT);
	specular_albedo->setTexParam(GL_TEXTURE_WRAP_T, GL_REPEAT);
	specular_albedo->setTexParamArr(GL_TEXTURE_MAX_ANISOTROPY, &max_aniso, 1);
	specular_albedo->unbind();

	scn->m_textures.push_back(std::move(specular_albedo));
	mat.m_specular_albedo = scn->m_textures.back().get();

	if (textures.find("aniso") == textures.end())
		throw std::logic_error("SCENE_LOADER: No aniso texture found in .axf file.");
	switch (texture_dims["aniso"].num_channels)
	{
	case 1:
		internalformat = GL_R32F;
		format = GL_RED;
		break;
	case 2:
		internalformat = GL_RG32F;
		format = GL_RG;
		break;
	case 3:
		internalformat = GL_RGB32F;
		format = GL_RGB;
		break;
	case 4:
		internalformat = GL_RGBA32F;
		format = GL_RGBA;
	default:
		throw std::logic_error("SCENE_LOADER: Invalid channel count for aniso.");
		break;
	}
	std::unique_ptr<Texture> aniso = Texture::T2DFromData(
		internalformat,
		texture_dims["aniso"].width,
		texture_dims["aniso"].height,
		format,
		GL_FLOAT,
		textures["aniso"].data(),
		true
	);

	aniso->bind(0);
	aniso->setTexParam(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	aniso->setTexParam(GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	aniso->setTexParam(GL_TEXTURE_WRAP_S, GL_REPEAT);
	aniso->setTexParam(GL_TEXTURE_WRAP_T, GL_REPEAT);
	aniso->setTexParamArr(GL_TEXTURE_MAX_ANISOTROPY, &max_aniso, 1);
	aniso->unbind();

	scn->m_textures.push_back(std::move(aniso));
	mat.m_aniso_rotation = scn->m_textures.back().get();

	if (textures.find("displacement") != textures.end())
	{
		switch (texture_dims["displacement"].num_channels)
		{
		case 1:
			internalformat = GL_R32F;
			format = GL_RED;
			break;
		case 2:
			internalformat = GL_RG32F;
			format = GL_RG;
			break;
		case 3:
			internalformat = GL_RGB32F;
			format = GL_RGB;
			break;
		case 4:
			internalformat = GL_RGBA32F;
			format = GL_RGBA;
		default:
			throw std::logic_error("SCENE_LOADER: Invalid channel count for displacement.");
			break;
		}
		std::unique_ptr<Texture> displacement = Texture::T2DFromData(
			internalformat,
			texture_dims["displacement"].width,
			texture_dims["displacement"].height,
			format,
			GL_FLOAT,
			textures["displacement"].data(),
			true
		);

		displacement->bind(0);
		displacement->setTexParam(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		displacement->setTexParam(GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		displacement->setTexParam(GL_TEXTURE_WRAP_S, GL_REPEAT);
		displacement->setTexParam(GL_TEXTURE_WRAP_T, GL_REPEAT);
		displacement->setTexParamArr(GL_TEXTURE_MAX_ANISOTROPY, &max_aniso, 1);
		displacement->unbind();

		scn->m_textures.push_back(std::move(displacement));
		mat.m_displacement = scn->m_textures.back().get();
	}
	else
	{
		mat.m_displacement = nullptr;
	}

	if (textures.find("fresnel") == textures.end())
		throw std::logic_error("SCENE_LOADER: No fresnel texture found in .axf file.");
	switch (texture_dims["fresnel"].num_channels)
	{
	case 1:
		internalformat = GL_R32F;
		format = GL_RED;
		break;
	case 2:
		internalformat = GL_RG32F;
		format = GL_RG;
		break;
	case 3:
		internalformat = GL_RGB32F;
		format = GL_RGB;
		break;
	case 4:
		internalformat = GL_RGBA32F;
		format = GL_RGBA;
	default:
		throw std::logic_error("SCENE_LOADER: Invalid channel count for fresnel.");
		break;
	}
	std::unique_ptr<Texture> fresnel = Texture::T2DFromData(
		internalformat,
		texture_dims["fresnel"].width,
		texture_dims["fresnel"].height,
		format,
		GL_FLOAT,
		textures["fresnel"].data(),
		true
	);

	fresnel->bind(0);
	fresnel->setTexParam(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	fresnel->setTexParam(GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	fresnel->setTexParam(GL_TEXTURE_WRAP_S, GL_REPEAT);
	fresnel->setTexParam(GL_TEXTURE_WRAP_T, GL_REPEAT);
	fresnel->setTexParamArr(GL_TEXTURE_MAX_ANISOTROPY, &max_aniso, 1);
	fresnel->unbind();

	scn->m_textures.push_back(std::move(fresnel));
	mat.m_fresnel_f0 = scn->m_textures.back().get();

	if (textures.find("lobes") == textures.end())
		throw std::logic_error("SCENE_LOADER: No roughness texture found in .axf file.");
	switch (texture_dims["lobes"].num_channels)
	{
	case 1:
		internalformat = GL_R32F;
		format = GL_RED;
		break;
	case 2:
		internalformat = GL_RG32F;
		format = GL_RG;
		break;
	case 3:
		internalformat = GL_RGB32F;
		format = GL_RGB;
		break;
	case 4:
		internalformat = GL_RGBA32F;
		format = GL_RGBA;
	default:
		throw std::logic_error("SCENE_LOADER: Invalid channel count for roughness.");
		break;
	}
	std::unique_ptr<Texture> roughness = Texture::T2DFromData(
		internalformat,
		texture_dims["lobes"].width,
		texture_dims["lobes"].height,
		format,
		GL_FLOAT,
		textures["lobes"].data(),
		true
	);

	roughness->bind(0);
	roughness->setTexParam(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	roughness->setTexParam(GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	roughness->setTexParam(GL_TEXTURE_WRAP_S, GL_REPEAT);
	roughness->setTexParam(GL_TEXTURE_WRAP_T, GL_REPEAT);
	roughness->setTexParamArr(GL_TEXTURE_MAX_ANISOTROPY, &max_aniso, 1);
	roughness->unbind();

	scn->m_textures.push_back(std::move(roughness));
	mat.m_roughness = scn->m_textures.back().get();

	if (textures.find("normal") == textures.end())
		throw std::logic_error("SCENE_LOADER: No normal texture found in .axf file.");
	switch (texture_dims["normal"].num_channels)
	{
	case 1:
		internalformat = GL_R32F;
		format = GL_RED;
		break;
	case 2:
		internalformat = GL_RG32F;
		format = GL_RG;
		break;
	case 3:
		internalformat = GL_RGB32F;
		format = GL_RGB;
		break;
	case 4:
		internalformat = GL_RGBA32F;
		format = GL_RGBA;
	default:
		throw std::logic_error("SCENE_LOADER: Invalid channel count for normal.");
		break;
	}
	std::unique_ptr<Texture> normal = Texture::T2DFromData(
		internalformat,
		texture_dims["normal"].width,
		texture_dims["normal"].height,
		format,
		GL_FLOAT,
		textures["normal"].data(),
		true
	);

	normal->bind(0);
	normal->setTexParam(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	normal->setTexParam(GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	normal->setTexParam(GL_TEXTURE_WRAP_S, GL_REPEAT);
	normal->setTexParam(GL_TEXTURE_WRAP_T, GL_REPEAT);
	normal->setTexParamArr(GL_TEXTURE_MAX_ANISOTROPY, &max_aniso, 1);
	normal->unbind();

	scn->m_textures.push_back(std::move(normal));
	mat.m_normals = scn->m_textures.back().get();

	if (textures.find("transparency") != textures.end())
	{
		switch (texture_dims["transparency"].num_channels)
		{
		case 1:
			internalformat = GL_R32F;
			format = GL_RED;
			break;
		case 2:
			internalformat = GL_RG32F;
			format = GL_RG;
			break;
		case 3:
			internalformat = GL_RGB32F;
			format = GL_RGB;
			break;
		case 4:
			internalformat = GL_RGBA32F;
			format = GL_RGBA;
		default:
			throw std::logic_error("SCENE_LOADER: Invalid channel count for transparency.");
			break;
		}
		std::unique_ptr<Texture> transparency = Texture::T2DFromData(
			internalformat,
			texture_dims["transparency"].width,
			texture_dims["transparency"].height,
			format,
			GL_FLOAT,
			textures["transparency"].data(),
			true
		);

		transparency->bind(0);
		transparency->setTexParam(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		transparency->setTexParam(GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		transparency->setTexParam(GL_TEXTURE_WRAP_S, GL_REPEAT);
		transparency->setTexParam(GL_TEXTURE_WRAP_T, GL_REPEAT);
		transparency->setTexParamArr(GL_TEXTURE_MAX_ANISOTROPY, &max_aniso, 1);
		transparency->unbind();

		scn->m_textures.push_back(std::move(transparency));
		mat.m_transparency = scn->m_textures.back().get();
	}
	else
	{
		mat.m_transparency = nullptr;
	}

	scn->m_materials.push_back(std::unique_ptr<Material>(new Material(mat)));
}