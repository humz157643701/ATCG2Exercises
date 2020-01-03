#ifndef _SYMMETRY_H_
#define _SYMMETRY_H_
#include <memory>
#include <Eigen/Dense>
#include <limits>
#include <icp.h>
#include <chrono>

struct SymmetryResult
{
	Eigen::Vector3d center;
	Eigen::Vector3d normal;
	Eigen::Vector3d trans;
	Eigen::Matrix3d rot;
};

template <typename MeshSampler>
class SymmetryDetector : private MeshSampler
{
public:
	SymmetryResult findMainSymmetryPlane(const Eigen::MatrixXd& vertices, const Eigen::MatrixXd& normals, const Eigen::MatrixXi& faces, std::size_t max_points_for_symmetry_calculation = 0)
	{
		// --- filter point cloud ---
		Eigen::MatrixXd Vt;
		Eigen::MatrixXd Nt;
		sampleMeshPoints(vertices, normals, faces, max_points_for_symmetry_calculation, Vt, Nt);

		// --- find symmetry plane ---
		Eigen::Vector3d center_of_mass = Vt.colwise().mean().transpose();
		Eigen::Vector3d plane_normal{ 1.0, 0.0, 0.0 };

		// icp instance
		ICPAligner icp(Vt);
		// query set copy
		Eigen::MatrixXd Vq(Vt);

		// construct initial plane and reflection matrix
		double origin_plane_distance = 0;
		Eigen::Matrix3d reflection_matrix;
		reflectionmatreflection_matrixrix << 1, 0, 0, 0, 1, 0, 0, 0, 1;
		reflectionmatrix = reflectionmatrix - 2 * (plane_normal * plane_normal.transpose());
		double origin_plane_distance = (-centerofmass).dot(plane_normal);
		
		// reflect query set across initial plane		
		Vq = Vq * reflectionmatrix.transpose();
		Vq += 2.0 * origin_plane_distance * plane_normal;

		// return result
	}
};

#endif