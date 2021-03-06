#ifndef _SYMMETRY_H_
#define _SYMMETRY_H_
#include <memory>
#include <Eigen/Dense>
#include <limits>
#include <icp.h>
#include <chrono>
#include <mesh.h>

struct SymmetryResult
{
	Eigen::Vector3d center;
	Eigen::Vector3d normal;
	Eigen::Vector3d trans;
	Eigen::Matrix3d rot;
	Eigen::MatrixXd refVt;
	Eigen::MatrixXd refVq;
	double time_for_prefiltering;
	double time_for_alignment;
	double time_for_sym_calculation;
	double time_total;
	double whole_mesh_mae;
};

void reflectAlongPlane(const SymmetryResult &sres, Eigen::MatrixXd& vertices) {
	Eigen::Matrix3d reflection_matrix;
	reflection_matrix << 1, 0, 0, 0, 1, 0, 0, 0, 1;
	reflection_matrix = reflection_matrix - 2 * (sres.normal * sres.normal.transpose());
	double origin_plane_distance = (sres.center).dot(sres.normal);

	// reflect query set across initial plane		
	vertices *= reflection_matrix.transpose();
	vertices.rowwise() += 2.0 * origin_plane_distance * sres.normal.transpose();
}

template <typename MeshSampler>
class SymmetryDetector
{
public:
	SymmetryDetector(const ICPParams& icpparams = {}, const MeshSampler& meshsampler = {}) :
		m_meshsampler(meshsampler),
		m_icp_params(icpparams)
	{}

	SymmetryResult findMainSymmetryPlane(const Mesh& mesh, const Eigen::Vector3d& initial_plane_normal)
	{
		auto t0 = std::chrono::high_resolution_clock::now();
		// --- sample points on mesh ---
		std::cout << "--- Symmetry detection: prefiltering input mesh...\n";
		Eigen::MatrixXd Vt;
		Eigen::MatrixXd Nt;
		auto t1 = std::chrono::high_resolution_clock::now();
		m_meshsampler.sampleMeshPoints(mesh, Vt, Nt);
		auto t_prefiltering = std::chrono::high_resolution_clock::now() - t1;

		// --- find symmetry plane ---
		Eigen::Vector3d center_of_mass{ Vt.colwise().mean().transpose() };
		Eigen::Vector3d plane_normal{ initial_plane_normal };

		// icp instance
		ICPAligner icp(Vt);
		// query set copy
		Eigen::MatrixXd Vq(Vt);
		Eigen::MatrixXd Nq(Nt);

		// construct initial plane and reflection matrix
		Eigen::Matrix3d reflection_matrix;
		reflection_matrix << 1, 0, 0, 0, 1, 0, 0, 0, 1;
		reflection_matrix = reflection_matrix - 2 * (plane_normal * plane_normal.transpose());
		double origin_plane_distance = (center_of_mass).dot(plane_normal);
		std::cout << center_of_mass;
		// reflect query set across initial plane		
		Vq *= reflection_matrix.transpose();
		Vq.rowwise() += 2.0 * origin_plane_distance * plane_normal.transpose();

		Nq *= reflection_matrix.transpose();
		Nq.rowwise().normalize();

		// align target and reflected query set
		Eigen::Matrix3d optimal_rotation;
		Eigen::Vector3d optimal_translation;
		//icp.align(optimal_rotation, optimal_translation, Vq, Vt, Nt, Nq, 50.0, 0.2, 1e-2, 1e-4, 100);
		std::cout << "--- Symmetry detection: calculating optimal rigid transform...\n";
		t1 = std::chrono::high_resolution_clock::now();
		icp.align(optimal_rotation, optimal_translation, Vq, Vt, Nt, Nq, m_icp_params);
		auto t_align = std::chrono::high_resolution_clock::now() - t1;

		std::cout << "--- Symmetry detection: eigendecomposition of optimally rotated reflection matrix...\n";
		t1 = std::chrono::high_resolution_clock::now();
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
		

		std::cout << "--- Symmetry detection: calculating result reflection plane...\n";
		Eigen::Vector3d newplanepoint = 0.5 * (optimal_rotation * (2 * origin_plane_distance * plane_normal) + optimal_translation);
		Eigen::Vector3d newnormal = es.eigenvectors().col(smallesteigenidx)(Eigen::all, 0).real();

		auto t_result = std::chrono::high_resolution_clock::now() - t1;

		reflection_matrix << 1, 0, 0, 0, 1, 0, 0, 0, 1;
		reflection_matrix = reflection_matrix - 2 * (newnormal * newnormal.transpose());
		origin_plane_distance = (-newplanepoint).dot(newnormal);

		Eigen::MatrixXd Vq_new(Vt);
		// reflect query set across initial plane		
		Vq_new *= reflection_matrix.transpose();
		Vq_new.rowwise() += 2.0 * origin_plane_distance * newnormal.transpose();

		// return result
		std::cout << "--- Symmetry detection: finished.\n";

		auto t_total = std::chrono::high_resolution_clock::now() - t0;

		return { newplanepoint,
			newnormal,
			optimal_translation,
			optimal_rotation,
			Vt,
			Vq_new,
			static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(t_prefiltering).count()) * 1e-6,
			static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(t_align).count()) * 1e-6,
			static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(t_result).count()) * 1e-6,
			static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(t_total).count()) * 1e-6,
			0.0
		};
	}

private:
	MeshSampler m_meshsampler;
	ICPParams m_icp_params;
};

#endif