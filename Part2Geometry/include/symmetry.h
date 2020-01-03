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
		Eigen::MatrixXd Nq(Nt);

		// construct initial plane and reflection matrix
		double origin_plane_distance = 0;
		Eigen::Matrix3d reflection_matrix;
		reflection_matrix << 1, 0, 0, 0, 1, 0, 0, 0, 1;
		reflection_matrix = reflection_matrix - 2 * (plane_normal * plane_normal.transpose());
		double origin_plane_distance = (-centerofmass).dot(plane_normal);
		
		// reflect query set across initial plane		
		Vq *= reflection_matrix.transpose();
		Vq += 2.0 * origin_plane_distance * plane_normal;

		Nq *= reflection_matrix.transpose();
		Nq.rowwise().normalize();

		// align target and reflected query set
		Eigen::Matrix3d optimal_rotation;
		Eigen::Vector3d optimal_translation;
		icp.align(optimal_rotation, optimal_translation, Vq, Vt, Nt, Nq, 50.0, 0.2, 1e-2, 1e-4, 100);

		auto reflected_rot = reflection_matrix * optimal_rotation.transpose();
		Eigen::EigenSolver<Eigen::MatrixXd> es(reflected_rot);
		int smallesteigenidx = 0;
		double smallesteigendiff = 1;
		for (size_t x = 0; x < es.eigenvalues().rows(); ++x)
		{
			auto diff = es.eigenvalues()(x, 0).real() + 1;
			if (diff < smallesteigendiff)
			{
				smallesteigenidx = x;
				smallesteigendiff = diff;
			}
		}

		Eigen::Vector3d newplanepoint = 0.5 * (optimal_rotation * (2 * origin_plane_distance * plane_normal) + optimal_translation);
		Eigen::Vector3d newnormal = es.eigenvectors().col(smallesteigenidx)(Eigen::all, 0).real();

		// return result
		return { newplanepoint, newnormal, optimal_translation, optimal_rotation };
	}
};

#endif