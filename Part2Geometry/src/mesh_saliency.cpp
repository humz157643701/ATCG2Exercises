#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <Octree.h>
#include <cmath>
#include <mesh.h>


void calculateMeshSaliency(const Mesh& mesh, double scale_base, std::size_t start_scale, std::size_t end_scale, Eigen::VectorXd& vertex_saliency)
{
	vertex_saliency.resize(mesh.vertices().rows());
	vertex_saliency.setZero();
	// first compute mean curvature
	Eigen::SparseMatrix<double> L, M, Minv;
	igl::cotmatrix(mesh.vertices(), mesh.faces(), L);
	igl::massmatrix(mesh.vertices(), mesh.faces(), igl::MASSMATRIX_TYPE_VORONOI, M);
	igl::invert_diag(M, Minv);
	Eigen::MatrixXd mean_curvature_normals = -Minv * (L * mesh.vertices());
	Eigen::VectorXd mean_curvatures = mean_curvature_normals.rowwise().norm();

	std::size_t numscales = static_cast<size_t>(end_scale - start_scale + 1);
	// holds saliency scores for each vertex and each scale
	Eigen::MatrixXd vertex_saliencies(mesh.vertices().rows(), static_cast<Eigen::DenseIndex>(numscales));

	// calculate aabb of vertices
	Eigen::Vector3d aabb_min = mesh.vertices().colwise().minCoeff().transpose();
	Eigen::Vector3d aabb_max = mesh.vertices().colwise().maxCoeff().transpose();
	double aabb_diag = (aabb_max - aabb_min).norm();
	double epsilon = aabb_diag * scale_base;

	// iterate through all scales and compute per-scale vertex saliencies
	Eigen::VectorXd G_fine(mesh.vertices().rows());
	Eigen::VectorXd G_coarse(mesh.vertices().rows());
	for (std::size_t scale = start_scale; scale <= end_scale; ++scale)
	{
		double cur_sigma = static_cast<double>(scale) * epsilon;
		double cur_sigma2 = 2.0 * cur_sigma;
		for (Eigen::DenseIndex v = 0; v < mesh.vertices().rows(); ++v)
		{
			// search for neighbours using octree
			auto res_fine = mesh.octree().query_radius(Vec3(mesh.vertices()(v, 0), mesh.vertices()(v, 1), mesh.vertices()(v, 2)), 2.0 * cur_sigma);
			auto res_coarse = mesh.octree().query_radius(Vec3(mesh.vertices()(v, 0), mesh.vertices()(v, 1), mesh.vertices()(v, 2)), 2.0 * cur_sigma2);
			// calculate G_fine; gaussian weighted average of mean curvature with sdev scale * epsilon
			double wfine = 0.0;
			double gfine = 0.0;
			for (size_t i = 0; i < res_fine.idx_dist_pair.size(); ++i)
			{
				Eigen::RowVector3d distvec = mesh.vertices()(res_fine.idx_dist_pair[i].first, Eigen::all) - mesh.vertices()(v, Eigen::all);
				double d2 = distvec.dot(distvec);
				double w = std::exp(-d2 / (2.0 * cur_sigma * cur_sigma));
				gfine += mean_curvatures(res_fine.idx_dist_pair[i].first) * w;
				wfine += w;
			}
			if (wfine != 0.0)
				G_fine(v) = gfine / wfine;
			else
				G_fine(v) = 0.0;
			// calculate G_coarse; gaussian weighted average of mean curvature with sdev 2 * scale * epsilon
			double wcoarse = 0.0;
			double gcoarse = 0.0;
			for (size_t i = 0; i < res_coarse.idx_dist_pair.size(); ++i)
			{
				Eigen::RowVector3d distvec = mesh.vertices()(res_coarse.idx_dist_pair[i].first, Eigen::all) - mesh.vertices()(v, Eigen::all);
				double d2 = distvec.dot(distvec);
				double w = std::exp(-d2 / (2.0 * cur_sigma2 * cur_sigma2));
				gcoarse += mean_curvatures(res_coarse.idx_dist_pair[i].first) * w;
				wcoarse += w;
			}
			if (wcoarse != 0.0)
				G_coarse(v) = gcoarse / wcoarse;
			else
				G_coarse(v) = 0.0;
		}
		// calculate vertex saliencies for this scale (DoG's)
		vertex_saliencies(Eigen::all, scale - start_scale) = (G_fine - G_coarse).cwiseAbs();
	}

	// vertex saliencies now represents unnormalized DoG scale space of mean curvature

	// normalize scale saliencies
	for (std::size_t scale = start_scale; scale <= end_scale; ++scale)
	{
		// for each scale calculate mean local maxima and global maximum
		// global maximum
		Eigen::DenseIndex global_maximum_idx;
		double global_maximum = vertex_saliencies(Eigen::all, static_cast<Eigen::DenseIndex>(scale - start_scale)).maxCoeff(&global_maximum_idx);

		// local maxima: do simple non maximum suppression based on one ring neighbourhood
		double mean_local_maximum = 0.0;
		size_t num_local_maxima = 0;

		for (Eigen::DenseIndex v = 0; v < vertex_saliencies.rows(); ++v)
		{
			// if center vertex v is greater than neighbourhood, it is a local maximum
			if (v != global_maximum_idx)
			{
				bool is_local_maximum = true;
				for (std::size_t i = 0; i < mesh.adjacency_list()[static_cast<std::size_t>(v)].size(); ++i)
				{
					if (vertex_saliencies(mesh.adjacency_list()[static_cast<std::size_t>(v)][i], static_cast<Eigen::DenseIndex>(scale - start_scale)) > vertex_saliencies(v, static_cast<Eigen::DenseIndex>(scale - start_scale)))
					{
						is_local_maximum = false;
						break;
					}
				}
				if (is_local_maximum)
				{
					mean_local_maximum += vertex_saliencies(v, static_cast<Eigen::DenseIndex>(scale - start_scale));
					num_local_maxima++;
				}
			}
		}

		if (num_local_maxima > 0)
			mean_local_maximum /= static_cast<double>(num_local_maxima);
		else
			mean_local_maximum = 0.0;

		// normalize this scale's saliency
		vertex_saliencies(Eigen::all, static_cast<Eigen::DenseIndex>(scale - start_scale)) *= ((global_maximum - mean_local_maximum) * (global_maximum - mean_local_maximum));
	}
	// sum up saliencies for all scales
	Eigen::VectorXd vertex_saliency = vertex_saliencies.rowwise().sum();
}

