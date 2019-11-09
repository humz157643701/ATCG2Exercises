#include "PerfDisplay.h"
// thats a bit weird, isn't it?
#include <ft2build.h>
#include FT_FREETYPE_H
#include <cmath>
#include <algorithm>
#include <vector>


TextRenderer::TextRenderer(const std::string & fontpath, size_t baseFontHeight, size_t maxlength) :
	m_charset_tex(nullptr),
	m_text_shader(ShaderProgram::createShaderProgram("assets/shaders/textrendering/textshader.vert", "assets/shaders/textrendering/textshader.frag")),
	m_max_length(maxlength)
{
	//--------------------------------------------- load font -------------------------------------------------
	FT_Library ft;
	if (FT_Init_FreeType(&ft))
		throw std::logic_error("FreeType library couldn't be initialized.");

	FT_Face face;
	if(FT_New_Face(ft, fontpath.c_str(), 0, &face))
		throw std::logic_error("Error while loading font.");

	if(FT_Set_Pixel_Sizes(face, 0, baseFontHeight))
		throw std::logic_error("Error while setting font height.");

	// load glyph format
	for (GLubyte c = 0; c < 128; ++c)
	{
		if(FT_Load_Char(face, c, FT_LOAD_RENDER))
			throw std::logic_error("Glyph couldn't be loaded.");
		m_glyphs[static_cast<size_t>(c)] = Char{
			{face->glyph->bitmap.width, face->glyph->bitmap.rows},
			{face->glyph->bitmap_left, face->glyph->bitmap_top},
			face->glyph->advance.x >> 6
		};
	}

	// generate texture array
	size_t maxwidth = 0;
	size_t maxheight = 0;
	for (size_t i = 0; i < 128; ++i)
	{
		if (m_glyphs[i].size.x > maxwidth)
			maxwidth = m_glyphs[i].size.x;
		if (m_glyphs[i].size.y > maxheight)
			maxheight = m_glyphs[i].size.y;
	}

	m_max_glyph_size = glm::ivec2{ maxwidth, maxheight };
	size_t texwidth = 1 << static_cast<size_t>(std::ceil(std::log2(maxwidth)));
	size_t texheight = 1 << static_cast<size_t>(std::ceil(std::log2(maxheight)));
	m_texsize = glm::ivec2{ texwidth, texheight };

	GLsizei miplevels = static_cast<GLsizei>(std::floor(std::log2(std::max(texwidth, texheight)))) + 1;
	m_charset_tex = std::move(Texture::T2DArr(GL_R8, texwidth, texheight, 128, miplevels));

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	m_charset_tex->bind(0);
	for (GLubyte c = 0; c < 128; ++c)
	{
		if (FT_Load_Char(face, c, FT_LOAD_RENDER))
			throw std::logic_error("Glyph couldn't be loaded.");
		glTexSubImage3D(m_charset_tex->getTarget(), 0, 0, 0, static_cast<GLint>(c),
			face->glyph->bitmap.width, face->glyph->bitmap.rows, 1,
			GL_RED, GL_UNSIGNED_BYTE, face->glyph->bitmap.buffer);
	}
	glGenerateMipmap(m_charset_tex->getTarget());
	m_charset_tex->setTexParam(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	m_charset_tex->setTexParam(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	m_charset_tex->setTexParam(GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	m_charset_tex->setTexParam(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	m_charset_tex->unbind();
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);

	FT_Done_Face(face);
	FT_Done_FreeType(ft);

	//----------------------------------------- buffers ----------------------------------------------
	
	glGenVertexArrays(1, &m_vao);
	if (m_vao == 0)
	{
		throw std::logic_error("VAO could not be created.");
	}
	glGenBuffers(1, &m_quad_buffer);
	if (m_quad_buffer == 0)
	{
		glDeleteVertexArrays(1, &m_vao);
		throw std::logic_error("Quad buffer could not be created.");
	}
	glGenBuffers(1, &m_pos_buffer);
	if (m_pos_buffer == 0)
	{
		glDeleteVertexArrays(1, &m_vao);
		glDeleteBuffers(1, &m_quad_buffer);
		throw std::logic_error("Pos buffer couldn't be created.");
	}
	glGenBuffers(1, &m_char_index_buffer);
	if (m_char_index_buffer == 0)
	{
		glDeleteVertexArrays(1, &m_vao);
		glDeleteBuffers(1, &m_quad_buffer);
		glDeleteBuffers(1, &m_pos_buffer);
		throw std::logic_error("Char index buffer couldn't be created.");
	}
	glGenBuffers(1, &m_align_buffer);
	if (m_align_buffer == 0)
	{
		glDeleteVertexArrays(1, &m_vao);
		glDeleteBuffers(1, &m_quad_buffer);
		glDeleteBuffers(1, &m_pos_buffer);
		glDeleteBuffers(1, &m_char_index_buffer);
		throw std::logic_error("Char align buffer couldn't be created.");
	}

	float quaddata[]{
		0.0f, 1.0f,
		0.0f, 0.0f,
		1.0f, 0.0f,

		1.0f, 0.0f,
		1.0f, 1.0f,
		0.0f, 1.0f
	};

	//allocate storage
	glBindBuffer(GL_ARRAY_BUFFER, m_quad_buffer);
	glBufferStorage(GL_ARRAY_BUFFER, sizeof(quaddata), &quaddata[0], 0);

	glBindBuffer(GL_ARRAY_BUFFER, m_pos_buffer);
	glBufferStorage(GL_ARRAY_BUFFER, sizeof(GLfloat) * 4 * maxlength, reinterpret_cast<const void*>(0), GL_DYNAMIC_STORAGE_BIT);

	glBindBuffer(GL_ARRAY_BUFFER, m_char_index_buffer);
	glBufferStorage(GL_ARRAY_BUFFER, sizeof(GLuint) * maxlength, reinterpret_cast<const void*>(0), GL_DYNAMIC_STORAGE_BIT);

	GLfloat uvaligndata[4 * 128];

	// calc uv scales
	for (size_t c = 0; c < 128; ++c)
	{
		uvaligndata[4 * c]		= static_cast<float>(m_glyphs[c].size.x) / static_cast<float>(texwidth);
		uvaligndata[4 * c + 1]	= static_cast<float>(m_glyphs[c].size.y) / static_cast<float>(texheight);
		// padding. Use some magic numbers for easy debugging.
		uvaligndata[4 * c + 2] = 42.6f;
		uvaligndata[4 * c + 3] = 42.9f;
	}

	glBindBuffer(GL_UNIFORM_BUFFER, m_align_buffer);
	glBufferStorage(GL_UNIFORM_BUFFER, sizeof(GLfloat) * 4 * 128, reinterpret_cast<const void*>(&uvaligndata[0]), 0);
	glBindBuffer(GL_UNIFORM_BUFFER, 0);

	// setup pipeline input state
	glBindVertexArray(m_vao);

	glBindBuffer(GL_ARRAY_BUFFER, m_quad_buffer);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 2, reinterpret_cast<const void*>(0));
	glVertexAttribDivisor(0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, m_pos_buffer);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 4, reinterpret_cast<const void*>(0));
	glVertexAttribDivisor(1, 1);

	glBindBuffer(GL_ARRAY_BUFFER, m_char_index_buffer);
	glEnableVertexAttribArray(2);
	glVertexAttribIPointer(2, 1, GL_UNSIGNED_INT, sizeof(GLuint), reinterpret_cast<const void*>(0));
	glVertexAttribDivisor(2, 1);

	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindVertexArray(0);
}

TextRenderer::~TextRenderer()
{
	if (m_vao)
		glDeleteVertexArrays(1, &m_vao);
	if (m_quad_buffer)
		glDeleteBuffers(1, &m_quad_buffer);
	if (m_pos_buffer)
		glDeleteBuffers(1, &m_pos_buffer);
	if (m_char_index_buffer)
		glDeleteBuffers(1, &m_char_index_buffer);
	if (m_align_buffer)
		glDeleteBuffers(1, &m_align_buffer);
}

TextRenderer::TextRenderer(TextRenderer && other) :
	m_text_shader(std::move(other.m_text_shader)),
	m_charset_tex(std::move(other.m_charset_tex)),
	m_current_string(other.m_current_string),
	m_font_color(other.m_font_color),
	m_max_length(other.m_max_length),
	m_texsize(other.m_texsize),
	m_max_glyph_size(other.m_max_glyph_size),
	m_vao(other.m_vao),
	m_quad_buffer(other.m_quad_buffer),
	m_align_buffer(other.m_align_buffer),
	m_origin(other.m_origin),
	m_vscale(other.m_vscale)
{
	other.m_vao = 0;
	other.m_quad_buffer = 0;
	other.m_pos_buffer = 0;
	other.m_char_index_buffer = 0;
	other.m_align_buffer = 0;
}

TextRenderer & TextRenderer::operator=(TextRenderer && other)
{
	if (&other == this)
		return *this;

	m_text_shader = std::move(other.m_text_shader);
	m_charset_tex = std::move(other.m_charset_tex);
	m_current_string = other.m_current_string;
	m_font_color = other.m_font_color;
	m_max_length = other.m_max_length;
	m_texsize = other.m_texsize;
	m_max_glyph_size = other.m_max_glyph_size;
	m_vao = other.m_vao;
	m_quad_buffer = other.m_quad_buffer;
	m_align_buffer = other.m_align_buffer;
	m_origin = other.m_origin;
	m_vscale = other.m_vscale;

	other.m_vao = 0;
	other.m_quad_buffer = 0;
	other.m_pos_buffer = 0;
	other.m_char_index_buffer = 0;
	other.m_align_buffer = 0;

	return *this;
}

void TextRenderer::setString(const std::string & str, const glm::vec2& origin, float scale, float screenaspect, const glm::vec4& color, float linespacing)
{
	m_current_string = str.substr(0, m_max_length);
	m_origin = origin;
	m_vscale = scale;
	m_font_color = color;
	m_linespacing = linespacing;

	updateBuffers(screenaspect);
}

void TextRenderer::draw()
{
	m_text_shader.use();

	glBindBufferRange(GL_UNIFORM_BUFFER, 1, m_align_buffer, 0, sizeof(GLfloat) * 4 * 128);
	m_text_shader.setUniform("fontcolor", m_font_color);
	m_text_shader.bindTex("glyphtex", m_charset_tex.get());
	glEnable(GL_BLEND);
	glBlendEquation(GL_FUNC_ADD);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDisable(GL_DEPTH_TEST);

	glBindVertexArray(m_vao);
	glDrawArraysInstanced(GL_TRIANGLES, 0, 6, m_numinstances);
	glBindVertexArray(0);
}

void TextRenderer::updateBuffers(float screenaspect)
{
	// generate position buffer (origin + bearing + advance stuff, relative to [0,1]²)
	// generate char index buffer (iterate string and store respective texture indices in the buffer)
	std::vector<GLuint> charindices;
	charindices.reserve(m_max_length);
	std::vector<glm::vec4> positions;
	positions.reserve(m_max_length);
	glm::vec2 cpos = { m_origin.x / screenaspect, m_origin.y };
	float hscale = m_vscale / m_max_glyph_size.y;
	float wscale = hscale / screenaspect;
	m_numinstances = 0;
	for (GLubyte c : m_current_string)
	{
		if (c == static_cast<GLubyte>('\n'))
		{
			cpos.x = m_origin.x / screenaspect;
			cpos.y -= m_vscale * (1.0f + m_linespacing);
		}
		else
		{
			size_t ci = static_cast<size_t>(c);
			glm::vec4 pos;
			pos.x = cpos.x + static_cast<float>(m_glyphs[ci].bearing.x) * wscale;
			pos.y = cpos.y - static_cast<float>(m_glyphs[ci].size.y - m_glyphs[ci].bearing.y) * hscale;
			pos.z = m_glyphs[ci].size.x * wscale;
			pos.w = m_glyphs[ci].size.y * hscale;
			positions.push_back(pos);
			charindices.push_back(static_cast<GLuint>(ci));
			cpos.x += static_cast<float>(m_glyphs[ci].advance) * wscale;
			m_numinstances++;
		}
	}

	glBindBuffer(GL_ARRAY_BUFFER, m_pos_buffer);
	glBufferSubData(GL_ARRAY_BUFFER, 0, positions.size() * sizeof(GLfloat) * 4, positions.data());

	glBindBuffer(GL_ARRAY_BUFFER, m_char_index_buffer);
	glBufferSubData(GL_ARRAY_BUFFER, 0, charindices.size() * sizeof(GLuint), charindices.data());

	glBindBuffer(GL_ARRAY_BUFFER, 0);
}
