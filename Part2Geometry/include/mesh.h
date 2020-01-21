#ifndef _MESH_H_
#define _MESH_H_
#include <Eigen/Dense>
//#include <Octree.h>
#include <vector>
#include <igl/adjacency_list.h>
#include <nanoflann.hpp>

//#define MESH_OCTREE_LEAF_SIZE 5
using kdtree_t = nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixXd, 3, nanoflann::metric_L2>;

class Mesh
{
public:
	Mesh();

	Mesh(const Eigen::MatrixXd& _vertices,
		const Eigen::MatrixXd& _normals,
		const Eigen::MatrixXi& _faces);

	Mesh(const Eigen::MatrixXd& _vertices,
		const Eigen::MatrixXd& _normals,
		const Eigen::MatrixXi& _faces,
		const Eigen::MatrixXd& _colors);

	Mesh(Eigen::MatrixXd&& _vertices,
		Eigen::MatrixXd&& _normals,
		Eigen::MatrixXi&& _faces);

	Mesh(Eigen::MatrixXd&& _vertices,
		Eigen::MatrixXd&& _normals,
		Eigen::MatrixXi&& _faces,
		Eigen::MatrixXd& _colors);

	Mesh(const Mesh& _other);
	Mesh(Mesh&& _other);
	Mesh& operator=(const Mesh& _other);
	Mesh& operator=(Mesh&& _other);

	const Eigen::MatrixXd& vertices() const { return m_vertices; }
	const Eigen::MatrixXd& normals() const { return m_normals; }
	const Eigen::MatrixXi& faces() const { return m_faces; }
	const Eigen::MatrixXd& colors() const { return m_colors; }
	Eigen::MatrixXd& vertices() { return m_vertices; }
	Eigen::MatrixXd& normals() { return m_normals; }
	Eigen::MatrixXi& faces() { return m_faces; }
	Eigen::MatrixXd& colors() { return m_colors; }
	const std::vector<std::vector<Eigen::DenseIndex>>& adjacency_list() const { return m_adjacency_list; }
	const kdtree_t& kdtree() const { return *m_kdtree.get(); };
	const std::vector<std::vector<Eigen::DenseIndex>>& triangle_list() const { return m_triangle_list; }
	void recalculateKdTree();
	void recalculateTriangleList();
private:
	Eigen::MatrixXd m_vertices;
	Eigen::MatrixXi m_faces;
	Eigen::MatrixXd m_normals;
	Eigen::MatrixXd m_colors;
	std::vector<std::vector<Eigen::DenseIndex>> m_adjacency_list;
	std::vector<std::vector<Eigen::DenseIndex>> m_triangle_list;
	std::unique_ptr<kdtree_t> m_kdtree;
};


#endif