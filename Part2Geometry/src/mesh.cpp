#include "..\include\mesh.h"
#include <igl/adjacency_list.h>

Mesh::Mesh() :
	m_kdtree(nullptr),
	m_vertices(),
	m_normals(),
	m_faces(),
	m_colors()
{
}

Mesh::Mesh(const Eigen::MatrixXd & _vertices, const Eigen::MatrixXd & _normals, const Eigen::MatrixXi & _faces)	:
	m_kdtree(nullptr),
	m_vertices(_vertices),
	m_normals(_normals),
	m_faces(_faces),
	m_colors()
{
	// build kdtree
	m_kdtree.reset(new kdtree_t(3, m_vertices));
	m_kdtree->index->buildIndex();

	// build adjencency list (hope that thing works. documentation is GREAT!)
	igl::adjacency_list(_faces, m_adjacency_list);
}

Mesh::Mesh(const Eigen::MatrixXd& _vertices, const Eigen::MatrixXd& _normals, const Eigen::MatrixXi& _faces, const Eigen::MatrixXd& _colors) :
	m_kdtree(nullptr),
	m_vertices(_vertices),
	m_normals(_normals),
	m_faces(_faces),
	m_colors(_colors)
{
	// build kdtree
	m_kdtree.reset(new kdtree_t(3, m_vertices));
	m_kdtree->index->buildIndex();

	// build adjencency list (hope that thing works. documentation is GREAT!)
	igl::adjacency_list(_faces, m_adjacency_list);
}

Mesh::Mesh(Eigen::MatrixXd && _vertices, Eigen::MatrixXd && _normals, Eigen::MatrixXi && _faces) :
	m_kdtree(nullptr),
	m_vertices(std::move(_vertices)),
	m_normals(std::move(_normals)),
	m_faces(std::move(_faces))
{
	// build kdtree
	m_kdtree.reset(new kdtree_t(3, m_vertices));
	m_kdtree->index->buildIndex();

	// build adjencency list (hope that thing works. documentation is GREAT!)
	igl::adjacency_list(_faces, m_adjacency_list);
}

Mesh::Mesh(Eigen::MatrixXd&& _vertices, Eigen::MatrixXd&& _normals, Eigen::MatrixXi&& _faces, Eigen::MatrixXd& _colors) :
	m_kdtree(nullptr),
	m_vertices(std::move(_vertices)),
	m_normals(std::move(_normals)),
	m_faces(std::move(_faces)),
	m_colors(std::move(_colors))
{
	// build kdtree
	m_kdtree.reset(new kdtree_t(3, m_vertices));
	m_kdtree->index->buildIndex();

	// build adjencency list (hope that thing works. documentation is GREAT!)
	igl::adjacency_list(_faces, m_adjacency_list);
}

Mesh::Mesh(const Mesh& _other) :
	m_kdtree(nullptr),
	m_vertices(_other.m_vertices),
	m_normals(_other.m_normals),
	m_faces(_other.m_faces),
	m_colors(_other.m_colors)
{
	m_kdtree.reset(new kdtree_t(3, m_vertices));
	m_kdtree->index->buildIndex();

	// build adjencency list (hope that thing works. documentation is GREAT!)
	igl::adjacency_list(m_faces, m_adjacency_list);
}

Mesh::Mesh(Mesh&& _other) :
	m_vertices(std::move(_other.m_vertices)),
	m_normals(std::move(_other.m_normals)),
	m_faces(std::move(_other.m_faces)),
	m_colors(std::move(_other.m_colors)),
	m_kdtree(std::move(_other.m_kdtree)),
	m_adjacency_list(std::move(_other.m_adjacency_list))
{
	_other.m_kdtree.reset();
}

Mesh& Mesh::operator=(const Mesh& _other)
{
	if (this == &_other)
		return *this;

	m_vertices = _other.m_vertices;
	m_normals = _other.m_normals;
	m_faces = _other.m_faces;
	m_colors = _other.m_colors;

	m_kdtree.reset(new kdtree_t(3, m_vertices));
	m_kdtree->index->buildIndex();

	// build adjencency list (hope that thing works. documentation is GREAT!)
	igl::adjacency_list(m_faces, m_adjacency_list);

	return *this;
}

Mesh& Mesh::operator=(Mesh&& _other)
{
	if (this == &_other)
		return *this;

	m_vertices = std::move(_other.m_vertices);
	m_normals = std::move(_other.m_normals);
	m_faces = std::move(_other.m_faces);
	m_colors = std::move(_other.m_colors);
	m_kdtree = std::move(_other.m_kdtree);
	m_adjacency_list = std::move(_other.m_adjacency_list);

	_other.m_kdtree.reset();
	 
	return *this;
}

void Mesh::recalculateKdTree()
{
	m_kdtree.reset(new kdtree_t(3, m_vertices));
	m_kdtree->index->buildIndex();
}
