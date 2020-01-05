#include "..\include\mesh.h"
#include <igl/adjacency_list.h>

Mesh::Mesh(const Eigen::MatrixXd & _vertices, const Eigen::MatrixXd & _normals, const Eigen::MatrixXi & _faces)	:
	m_octree(nullptr),
	m_vertices(_vertices),
	m_normals(_normals),
	m_faces(_faces)
{
	// build octree
	m_octree_data.reserve(static_cast<std::size_t>(_vertices.rows()));
	for (Eigen::DenseIndex i = 0; i < _vertices.rows(); ++i)
		m_octree_data.push_back(Vec3{ static_cast<float>(_vertices(i, 0)), static_cast<float>(_vertices(i, 1)), static_cast<float>(_vertices(i, 2)) });
	m_octree.reset(new Octree(m_octree_data));
	m_octree->build();

	// build adjencency list (hope that thing works. documentation is GREAT!)
	igl::adjacency_list(_faces, m_adjacency_list);
}

Mesh::Mesh(Eigen::MatrixXd && _vertices, Eigen::MatrixXd && _normals, Eigen::MatrixXi && _faces) :
	m_octree(nullptr),
	m_vertices(std::move(_vertices)),
	m_normals(std::move(_normals)),
	m_faces(std::move(_faces))
{
	// build octree
	m_octree_data.reserve(static_cast<std::size_t>(_vertices.rows()));
	for (Eigen::DenseIndex i = 0; i < _vertices.rows(); ++i)
		m_octree_data.push_back(Vec3{ static_cast<float>(_vertices(i, 0)), static_cast<float>(_vertices(i, 1)), static_cast<float>(_vertices(i, 2)) });
	m_octree.reset(new Octree(m_octree_data));
	m_octree->build(MESH_OCTREE_LEAF_SIZE);

	// build adjencency list (hope that thing works. documentation is GREAT!)
	igl::adjacency_list(_faces, m_adjacency_list);
}
