#include "Mesh.h"

Mesh::~Mesh()
{
	if (vao)
		glDeleteVertexArrays(1, &vao);
	if (vbo)
		glDeleteBuffers(1, &vbo);
	if (ibo)
		glDeleteBuffers(1, &ibo);
}

Mesh::Mesh(Mesh && other) :
	vbo(other.vbo),
	ibo(other.ibo),
	vao(other.vao),
	material(other.material),
	numIndices(other.numIndices),
	mscale(other.mscale)
{
	other.vbo = 0;
	other.ibo = 0;
	other.vao = 0;
	other.numIndices = 0;
}

Mesh & Mesh::operator=(Mesh && other)
{
	if (this == &other)
		return *this;

	if (vao)
		glDeleteVertexArrays(1, &vao);
	if (vbo)
		glDeleteBuffers(1, &vbo);
	if (ibo)
		glDeleteBuffers(1, &ibo);

	vao = other.vao;
	vbo = other.vbo;
	ibo = other.ibo;
	numIndices = other.numIndices;
	mscale = other.mscale;

	other.vbo = 0;
	other.ibo = 0;
	other.vao = 0;
	other.numIndices = 0;

	return *this;
}

void Mesh::drawWithMaterial(ShaderProgram * shader, size_t rid)
{
	assert(material);
	if (rid == 0 || material->renderer_id == rid)
	{
		shader->saveTU();
		material->bind(shader, mscale);
		draw();
		shader->restoreTU();
	}
}

void Mesh::draw(size_t rid)
{
	assert(material);
	if (material == nullptr || rid == 0 || material->renderer_id == rid)
	{
		if (vao)
		{
			glBindVertexArray(vao);
			glDrawElements(GL_TRIANGLES, numIndices, GL_UNSIGNED_INT, 0);
			glBindVertexArray(0);
		}
	}
}

void Mesh::setMaterial(Material* m)
{
	material = m;
}

std::unique_ptr<Mesh> Mesh::createMesh(const std::vector<Vertex>& vd, const std::vector<Index>& id, const std::vector<VertexAttribute>& att)
{
	GLuint vbo, ibo, vao;
	glGenVertexArrays(1, &vao);
	if (vao == 0)
	{
		throw std::logic_error("VAO could not be created.");
	}
	glGenBuffers(1, &vbo);
	if (vbo == 0)
	{
		glDeleteVertexArrays(1, &vao);
		throw std::logic_error("VBO could not be created.");
	}
	glGenBuffers(1, &ibo);
	if (ibo == 0)
	{
		glDeleteVertexArrays(1, &vao);
		glDeleteBuffers(1, &vbo);
		throw std::logic_error("IBO could not be created.");
	}

	glBindVertexArray(vao);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(Vertex) * vd.size(), reinterpret_cast<const GLvoid*>(vd.data()), GL_STATIC_DRAW);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(Index) * id.size(), reinterpret_cast<const GLvoid*>(id.data()), GL_STATIC_DRAW);
	for (size_t i = 0; i < att.size(); ++i)
	{
		glEnableVertexAttribArray(static_cast<GLuint>(i));
		glVertexAttribPointer(
			static_cast<GLuint>(i),
			att[i].n,
			att[i].type,
			false,
			att[i].stride,
			reinterpret_cast<const void*>(att[i].offset)
		);
	}
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
	return std::unique_ptr<Mesh>(new Mesh(vbo, ibo, vao, static_cast<GLuint>(id.size()), nullptr));
}

std::unique_ptr<Mesh> Mesh::createMesh(const std::vector<Vertex>& vd, const std::vector<Index>& id, const std::vector<VertexAttribute>& att, Material* material, const glm::vec2& material_scale)
{
	GLuint vbo, ibo, vao;
	glGenVertexArrays(1, &vao); 
	if (vao == 0)
	{
		throw std::logic_error("VAO could not be created.");
	}
	glGenBuffers(1, &vbo);
	if (vbo == 0)
	{
		glDeleteVertexArrays(1, &vao);
		throw std::logic_error("VBO could not be created.");
	}
	glGenBuffers(1, &ibo);
	if (ibo == 0)
	{
		glDeleteVertexArrays(1, &vao);
		glDeleteBuffers(1, &vbo);
		throw std::logic_error("IBO could not be created.");
	}

	glBindVertexArray(vao);	
	glBindBuffer(GL_ARRAY_BUFFER, vbo); 
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo); 
	glBufferData(GL_ARRAY_BUFFER, sizeof(Vertex) * vd.size(), reinterpret_cast<const GLvoid*>(vd.data()), GL_STATIC_DRAW); 
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(Index) * id.size(), reinterpret_cast<const GLvoid*>(id.data()), GL_STATIC_DRAW);	
	for (size_t i = 0; i < att.size(); ++i)
	{
		glEnableVertexAttribArray(static_cast<GLuint>(i));
		glVertexAttribPointer(
			static_cast<GLuint>(i),
			att[i].n,
			att[i].type,
			GL_FALSE,
			att[i].stride,
			reinterpret_cast<const void*>(att[i].offset)
		);
	}
	glBindBuffer(GL_ARRAY_BUFFER, 0); 
	glBindVertexArray(0); 
	return std::unique_ptr<Mesh>(new Mesh(vbo, ibo, vao, static_cast<GLuint>(id.size()), material, material_scale));
}

Mesh::Mesh(GLuint _vbo, GLuint _ibo, GLuint _vao, GLsizei _numIndices, Material * _material, const glm::vec2& material_scale) :
	vbo(_vbo),
	ibo(_ibo),
	vao(_vao),
	numIndices(_numIndices),
	material(_material),
	mscale(material_scale)
{
}
