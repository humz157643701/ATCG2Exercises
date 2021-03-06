#ifndef _MESH_SALIENCY_H_
#define _MESH_SALIENCY_H_
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <Octree.h>
#include <cmath>
#include <mesh.h>

enum class ScaleType
{
	DOUBLE_SIGMA_EVERY_SCALE,
	LINEAR_INCREASE
};

void calculateMeshSaliency(const Mesh& mesh, double scale_base, std::size_t start_scale, std::size_t end_scale, Eigen::VectorXd& vertex_saliency, ScaleType scale_type = ScaleType::LINEAR_INCREASE);

namespace MeshSamplers
{
	struct MeshSaliencySampler
	{
		MeshSaliencySampler(double scale_base, std::size_t start_scale, std::size_t end_scale, bool normalize = false, double max_salient_points_fraction = 0.1, ScaleType scale_type = ScaleType::LINEAR_INCREASE, bool visualize = false) :
			m_scale_base(scale_base),
			m_start_scale(start_scale),
			m_end_scale(end_scale),
			m_normalize(normalize),
			m_max_sal_point_fraction(max_salient_points_fraction),
			m_scale_type(scale_type),
			m_visualize(visualize)
		{}

		void sampleMeshPoints(
			const Mesh& mesh,
			Eigen::MatrixXd& sampled_points,
			Eigen::MatrixXd& sampled_normals
		);		

		void sampleMeshPoints(
			const Mesh& mesh,
			Eigen::MatrixXd& sampled_points,
			Eigen::MatrixXd& sampled_normals,
			Eigen::VectorXd& mesh_saliency
		);

	private:
		double m_scale_base;
		std::size_t m_start_scale;
		std::size_t m_end_scale;
		bool m_normalize;
		double m_max_sal_point_fraction;
		ScaleType m_scale_type;
		bool m_visualize;
	};
}

#endif