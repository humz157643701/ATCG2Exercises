#ifndef _MESH_SALIENCY_H_
#define _MESH_SALIENCY_H_
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <Octree.h>

namespace MeshSamplers
{
	// default sampler

	struct MeshSaliencySampler
	{
		MeshSaliencySampler(double neighborhood_radius, double scale_base, std::size_t start_scale, std::size_t end_scale) :
			m_scale_base(scale_base),
			m_start_scale(start_scale),
			m_end_scale(end_scale)
		{}

		void sampleMeshPoints(
			const Eigen::MatrixXd& vertices,
			const Eigen::MatrixXd& normals,
			const Eigen::MatrixXi& faces,
			Eigen::MatrixXd& sampled_points,
			Eigen::MatrixXd& sampled_normals
		)
		{
			// first compute mean curvature
			Eigen::SparseMatrix<double> L, M, Minv;
			igl::cotmatrix(vertices, faces, L);
			igl::massmatrix(vertices, faces, igl::MASSMATRIX_TYPE_VORONOI, M);
			igl::invert_diag(M, Minv);
			Eigen::MatrixXd mean_curvature_normals = -Minv * (L * vertices);
			Eigen::VectorXd mean_curvatures = mean_curvature_normals.rowwise().norm();

			// build octree for neighbourhood search
			std::vector<Vec3> octree_data(vertices.rows());
			for (Eigen::DenseIndex i = 0; i < vertices.rows(); ++i)
				octree_data.push_back(Vec3{ static_cast<float>(vertices(i, 0)), static_cast<float>(vertices(i, 1)), static_cast<float>(vertices(i, 2)) });
			Octree octr(octree_data);
			octr.build();

			std::size_t numscales = static_cast<size_t>(m_end_scale - m_start_scale + 1);
			Eigen::MatrixXd vertex_saliencies(vertices.rows(), static_cast<Eigen::DenseIndex>(numscales));
			
			// calculate aabb of vertices
			Eigen::Vector3d aabb_min = vertices.colwise().minCoeff().transpose();
			Eigen::Vector3d aabb_max = vertices.colwise().maxCoeff().transpose();
			double aabb_diag = (aabb_max - aabb_min).norm();
			double epsilon = aabb_diag * m_scale_base;

			// iterate through all scales and compute per-scale vertex saliencies
			Eigen::VectorXd G_fine(vertices.rows());
			Eigen::VectorXd G_coarse(vertices.rows());
			for (std::size_t scale = m_start_scale; scale <= m_end_scale; ++scale)
			{
				// calculate G_fine; gaussian weighted average of mean curvature with sdev scale * epsilon

				// calculate G_coarse; gaussian weighted average of mean curvature with sdev 2 * scale * epsilon

				// calculate saliency for this scale (G_fine - G_coarse)
			}

			// sum up saliencies for all scales (using non maximum suppression)

			// search for vertices with high saliency
		}

	private:
		double m_scale_base;
		std::size_t m_start_scale;
		std::size_t m_end_scale;
	};
}

#endif