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

struct TexEntry
{
	Texture tex;
	std::string path;
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
		else if (scan->lookahead().value == "AXF_MATERIAL")
		{
			scan->match(TokType::Name);
			parseMaterial(scan, scn);
		}
		else if (scan->lookahead().value == "MAX_LUM")
		{
			scan->match(TokType::Name);
			parseMaxLum(scan, scn);
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

	scn->m_camera.setParameters(800, 600, glm::radians(fov), near, far);
	scn->m_camera.transform.setTransformMatrix(glm::inverse(viewmat));
}

void SceneLoader::parseModelOpaque(Scanner * scan, Scene * scn)
{	
	std::string path;
	glm::vec3 pos;
	glm::vec3 scale;
	glm::quat rot;

	path = scan->lookahead().value; scan->match(TokType::String);
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

	bool recalcnormals = scan->lookahead().boolvalue; scan->match(TokType::Boolean);

	OBJResult res = OBJLoader::loadOBJ(path, recalcnormals, true);

	


}

void SceneLoader::parseModelTransparent(Scanner * scan, Scene * scn)
{
	std::string path;
	glm::vec3 pos;
	glm::vec3 scale;
	glm::quat rot;

	path =		scan->lookahead().value; scan->match(TokType::String);
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

	c.x = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	c.y = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	c.z = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	p.x = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	p.y = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	p.z = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	

	scn->m_pointlights.push_back(PointLight(Transform(p, glm::quat(), glm::vec3(1.0f)), c));
}

void SceneLoader::parseAmbientLight(Scanner * scan, Scene * scn)
{
	glm::vec3 color;
	color.x = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	color.y = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	color.z = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	scn->m_ambientlights.push_back(AmbientLight(color));
}

void SceneLoader::parseMaxLum(Scanner * scan, Scene * scn)
{
	float maxlum = static_cast<float>(scan->lookahead().floatvalue); scan->match(TokType::Float);
	scn->maxlum = maxlum;
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

}

void SceneLoader::parseMaterial(Scanner * scan, Scene * scn)
{

}