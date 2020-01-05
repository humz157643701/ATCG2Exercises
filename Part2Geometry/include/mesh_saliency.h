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

void calculateMeshSaliency(const Mesh& mesh, double scale_base, std::size_t start_scale, std::size_t end_scale, Eigen::VectorXd& vertex_saliency);

namespace MeshSamplers
{
	struct MeshSaliencySampler
	{
		MeshSaliencySampler(double scale_base, std::size_t start_scale, std::size_t end_scale) :
			m_scale_base(scale_base),
			m_start_scale(start_scale),
			m_end_scale(end_scale)
		{}

		void sampleMeshPoints(
			const Mesh& mesh,
			Eigen::MatrixXd& sampled_points,
			Eigen::MatrixXd& sampled_normals
		)
		{
			
		}

	private:
		double m_scale_base;
		std::size_t m_start_scale;
		std::size_t m_end_scale;
	};
}

#endif