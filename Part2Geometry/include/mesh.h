#ifndef _MESH_H_
#define _MESH_H_
#include <Eigen/Dense>
#include <Octree.h>
#include <vector>
#include <igl/adjacency_list.h>
#define MESH_OCTREE_LEAF_SIZE 5
class Mesh
{
public:
	Mesh(const Eigen::MatrixXd& _vertices,
		const Eigen::MatrixXd& _normals,
		const Eigen::MatrixXi& _faces);

	Mesh(Eigen::MatrixXd&& _vertices,
		Eigen::MatrixXd&& _normals,
		Eigen::MatrixXi&& _faces);

	const Eigen::MatrixXd& vertices() const { return m_vertices; }
	const Eigen::MatrixXd& normals() const { return m_normals; }
	const Eigen::MatrixXi& faces() const { return m_faces; }
	const std::vector<std::vector<Eigen::DenseIndex>>& adjacency_list() const { return m_adjacency_list; }
	const Octree& octree() const { return *m_octree.get(); };
private:
	Eigen::MatrixXd m_vertices;
	Eigen::MatrixXi m_faces;
	Eigen::MatrixXd m_normals;
	std::vector<std::vector<Eigen::DenseIndex>> m_adjacency_list;
	std::unique_ptr<Octree> m_octree;
	std::vector<Vec3> m_octree_data;
};


#endif