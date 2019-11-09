#ifndef _PERF_DISPLAY_H_
#define _PERF_DISPLAY_H_
#include <libheaders.h>
#include <Shader.h>
#include <memory>

//assume utf8 encoding, using only the ascii range [0, 127]

class TextRenderer
{
private:
	struct Char
	{
		glm::ivec2 size;
		glm::ivec2 bearing;
		GLint advance;
	};

public:
	TextRenderer(const std::string& fontpath, size_t baseFontHeight, size_t maxlength);
	~TextRenderer();

	TextRenderer(const TextRenderer&) = delete;
	TextRenderer& operator=(const TextRenderer&) = delete;

	TextRenderer(TextRenderer&& other);
	TextRenderer& operator=(TextRenderer&& other);

	void setString(const std::string & str, const glm::vec2& origin, float scale, float screenaspect, const glm::vec4& color = {0.0f, 0.0f, 0.0f, 1.0f}, float linespacing = 0.5f);

	void draw();

private:
	void updateBuffers(float screenaspect);
private:
	// array texture with 2d texture per glyph
	std::unique_ptr<Texture> m_charset_tex;
	std::string m_current_string;
	glm::vec4 m_font_color;
	Char m_glyphs[128];
	size_t m_max_length;
	glm::ivec2 m_texsize;
	glm::ivec2 m_max_glyph_size;
	GLuint m_numinstances;

	glm::vec2 m_origin;
	float m_vscale;
	float m_linespacing;

	GLuint m_vao;
	// simple quad vertices with canonical texture coordinates
	GLuint m_quad_buffer;
	// position + glyph offset including advance and bearing
	GLuint m_pos_buffer;
	// texture index for every character
	GLuint m_char_index_buffer;
	// u scale, v scale
	GLuint m_align_buffer;

	ShaderProgram m_text_shader;
};
#endif