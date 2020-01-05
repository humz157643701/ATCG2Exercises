#include "..\include\mesh.h"
#include <igl/adjacency_list.h>

Mesh::Mesh(const Eigen::MatrixXd & _vertices, const Eigen::MatrixXd & _normals, const Eigen::MatrixXi & _faces)	:
	m_kdtree(nullptr),
	m_vertices(_vertices),
	m_normals(_normals),
	m_faces(_faces)
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
