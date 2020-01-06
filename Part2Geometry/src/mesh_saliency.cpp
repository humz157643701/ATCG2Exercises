#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <cmath>
#include <mesh.h>
#include <iostream>
#include <algorithm>
#include <mesh_saliency.h>



void calculateMeshSaliency(const Mesh& mesh, double scale_base, std::size_t start_scale, std::size_t end_scale, Eigen::VectorXd& vertex_saliency, bool normalize)
{
	vertex_saliency.resize(mesh.vertices().rows());
	vertex_saliency.setZero();
	// first compute mean curvature
	std::cout << "Calculating mean curvature...\n";
	Eigen::SparseMatrix<double> L, M, Minv;
	igl::cotmatrix(mesh.vertices(), mesh.faces(), L);
	igl::massmatrix(mesh.vertices(), mesh.faces(), igl::MASSMATRIX_TYPE_VORONOI, M);
	igl::invert_diag(M, Minv);
	Eigen::MatrixXd mean_curvature_normals = -Minv * (L * mesh.vertices());
	Eigen::VectorXd mean_curvatures = mean_curvature_normals.rowwise().norm();

	//DEBUG
	//vertex_saliency = mean_curvatures;
	//return;

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
	std::cout << "Calculating per-scale saliencies...\n";
	std::vector<std::pair<long long, double>> rad_search_res_fine;
	std::vector<std::pair<long long, double>> rad_search_res_coarse;
	for (std::size_t scale = start_scale; scale <= end_scale; ++scale)
	{
		std::cout << "Scale " << std::to_string(scale) << "\n";
		double cur_sigma = static_cast<double>(scale) * epsilon;
		double cur_sigma2 = 2.0 * cur_sigma;
		
		for (Eigen::DenseIndex v = 0; v < mesh.vertices().rows(); ++v)
		{
			double query_point[] = { mesh.vertices()(v, 0), mesh.vertices()(v, 1), mesh.vertices()(v, 2) };
			// calculate G_fine; gaussian weighted average of mean curvature with sdev scale * epsilon			
			rad_search_res_fine.clear();			
			mesh.kdtree().index->radiusSearch(&query_point[0], (2.0 * cur_sigma) * (2.0 * cur_sigma), rad_search_res_fine, nanoflann::SearchParams(0, 0.0, false));			
			double wfine = 0.0;
			double gfine = 0.0;
			for (size_t i = 0; i < rad_search_res_fine.size(); ++i)
			{
				double w = std::exp(-rad_search_res_fine[i].second / (2.0 * cur_sigma * cur_sigma));
				gfine += mean_curvatures(rad_search_res_fine[i].first) * w;
				wfine += w;
			}
			if (wfine != 0.0)
				G_fine(v) = gfine / wfine;
			else
				G_fine(v) = 0.0;
			
			// calculate G_coarse; gaussian weighted average of mean curvature with sdev 2 * scale * epsilon
			rad_search_res_coarse.clear();
			mesh.kdtree().index->radiusSearch(&query_point[0], (2.0 * cur_sigma2) * (2.0 * cur_sigma2), rad_search_res_coarse, nanoflann::SearchParams(0, 0.0, false));
			double wcoarse = 0.0;
			double gcoarse = 0.0;
			for (size_t i = 0; i < rad_search_res_coarse.size(); ++i)
			{
				double w = std::exp(-rad_search_res_coarse[i].second / (2.0 * cur_sigma2 * cur_sigma2));
				gcoarse += mean_curvatures(rad_search_res_coarse[i].first) * w;
				wcoarse += w;
			}
			if (wcoarse != 0.0)
				G_coarse(v) = gcoarse / wcoarse;
			else
				G_coarse(v) = 0.0;
		}
		// calculate vertex saliencies for this scale (DoG's)
		vertex_saliencies(Eigen::all, scale - start_scale) = (G_fine - G_coarse).cwiseAbs();
		if (normalize)
			vertex_saliencies(Eigen::all, scale - start_scale) *= (cur_sigma * cur_sigma);
	}

	// vertex saliencies now represents unnormalized DoG scale space of mean curvature

	//DEBUG
	//vertex_saliency = vertex_saliencies(Eigen::all, 4);
	//return;

	// normalize scale saliencies
	std::cout << "Normalizing per-scale saliency...\n";
	for (std::size_t scale = start_scale; scale <= end_scale; ++scale)
	{
		std::cout << "Scale " << std::to_string(scale) << "\n";
		// for each scale calculate mean local maxima and global maximum
		// global maximum
		std::cout << "Calculating global maximum...\n";
		Eigen::DenseIndex global_maximum_idx;
		double global_maximum = vertex_saliencies(Eigen::all, static_cast<Eigen::DenseIndex>(scale - start_scale)).maxCoeff(&global_maximum_idx);

		// local maxima: do simple non maximum suppression based on one ring neighbourhood
		double mean_local_maximum = 0.0;
		size_t num_local_maxima = 0;

		std::cout << "Calculating mean local maximum...\n";
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
	std::cout << "Calculating final mesh saliency...\n";
	vertex_saliency = vertex_saliencies.rowwise().sum();
}

void MeshSamplers::MeshSaliencySampler::sampleMeshPoints(const Mesh & mesh, Eigen::MatrixXd & sampled_points, Eigen::MatrixXd & sampled_normals)
{
	Eigen::VectorXd mesh_saliency;
	calculateMeshSaliency(mesh, m_scale_base, m_start_scale, m_end_scale, mesh_saliency, m_normalize);

	// search local minima, sort by saliency and return the best 90% or so
	// local maxima: do simple non maximum suppression based on one ring neighbourhood
	size_t num_local_maxima = 0;
	std::vector<std::pair<Eigen::DenseIndex, double>> local_maxima(static_cast<size_t>(mesh.vertices().rows() / 10));

	for (Eigen::DenseIndex v = 0; v < mesh_saliency.rows(); ++v)
	{		
		bool is_local_maximum = true;
		for (std::size_t i = 0; i < mesh.adjacency_list()[static_cast<std::size_t>(v)].size(); ++i)
		{
			if (mesh_saliency(mesh.adjacency_list()[static_cast<std::size_t>(v)][i]) > mesh_saliency(v))
			{
				is_local_maximum = false;
				break;
			}
		}
		if (is_local_maximum)
		{
			local_maxima.push_back(std::make_pair(v, mesh_saliency(v)));
		}		
	}

	// sort local maxima by saliency score, descending order
	std::sort(local_maxima.begin(), local_maxima.end(), [](const std::pair<Eigen::DenseIndex, double>& a, const std::pair<Eigen::DenseIndex, double>& b) {
		return a.second > b.second;
	});

	std::size_t num_samples = std::min(local_maxima.size(), static_cast<size_t>(mesh.vertices().rows() * m_max_sal_point_fraction));

	sampled_points.resize(num_samples, 3);
	sampled_normals.resize(num_samples, 3);

	// return the best local minima points and normals
	for (Eigen::DenseIndex i = 0; i < static_cast<std::size_t>(num_samples); ++i)
	{
		sampled_points(i, Eigen::all) = mesh.vertices()(local_maxima[i].first, Eigen::all);
		sampled_normals(i, Eigen::all) = mesh.normals()(local_maxima[i].first, Eigen::all);
	}	
}

void MeshSamplers::MeshSaliencySampler::sampleMeshPoints(const Mesh & mesh, Eigen::MatrixXd & sampled_points, Eigen::MatrixXd & sampled_normals, Eigen::VectorXd & _mesh_saliency)
{
	Eigen::VectorXd mesh_saliency;
	calculateMeshSaliency(mesh, m_scale_base, m_start_scale, m_end_scale, mesh_saliency, m_normalize);

	// search local minima, sort by saliency and return the best 90% or so
	// local maxima: do simple non maximum suppression based on one ring neighbourhood
	size_t num_local_maxima = 0;
	std::vector<std::pair<Eigen::DenseIndex, double>> local_maxima(static_cast<size_t>(mesh.vertices().rows() / 10));

	for (Eigen::DenseIndex v = 0; v < mesh_saliency.rows(); ++v)
	{
		bool is_local_maximum = true;
		for (std::size_t i = 0; i < mesh.adjacency_list()[static_cast<std::size_t>(v)].size(); ++i)
		{
			if (mesh_saliency(mesh.adjacency_list()[static_cast<std::size_t>(v)][i]) > mesh_saliency(v))
			{
				is_local_maximum = false;
				break;
			}
		}
		if (is_local_maximum)
		{
			local_maxima.push_back(std::make_pair(v, mesh_saliency(v)));
		}
	}

	// sort local maxima by saliency score, descending order
	std::sort(local_maxima.begin(), local_maxima.end(), [](const std::pair<Eigen::DenseIndex, double>& a, const std::pair<Eigen::DenseIndex, double>& b) {
		return a.second > b.second;
	});

	std::size_t num_samples = std::min(local_maxima.size(), static_cast<size_t>(mesh.vertices().rows() * m_max_sal_point_fraction));

	sampled_points.resize(num_samples, 3);
	sampled_normals.resize(num_samples, 3);

	// return the best local minima points and normals
	for (Eigen::DenseIndex i = 0; i < static_cast<std::size_t>(num_samples); ++i)
	{
		sampled_points(i, Eigen::all) = mesh.vertices()(local_maxima[i].first, Eigen::all);
		sampled_normals(i, Eigen::all) = mesh.normals()(local_maxima[i].first, Eigen::all);
	}

	// return saliency scores for debug purposes
	_mesh_saliency.resize(mesh_saliency.rows());
	_mesh_saliency = mesh_saliency;
}
