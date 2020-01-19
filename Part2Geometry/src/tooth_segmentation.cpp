#include "..\include\tooth_segmentation.h"
#include <iostream>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/jet.h>
#include <vector>
#include <random>

void ToothSegmentation::segmentTeethFromMesh(const Mesh& mesh, const Eigen::Vector3d& mesh_up, const Eigen::Vector3d& mesh_right, std::vector<Mesh>& tooth_meshes, const ToothSegmentation::CuspDetectionParams& cuspd_params, bool visualize_steps)
{
	// transform mesh to canonical pose
	std::cout << "--- TOOTH SEGMENTATION ---\n";
	std::cout << "-- Transforming mesh into canonical pose...\n";
	Mesh working_mesh(mesh);
	Eigen::Matrix3d crot;
	crot.col(0) = mesh_right.normalized();
	crot.col(1) = mesh_up.normalized();
	crot.col(2) = mesh_right.cross(mesh_up).normalized();
	crot.transposeInPlace();
	Eigen::Vector3d ctrans = mesh.vertices().colwise().mean().transpose();

	working_mesh.vertices() *= crot.transpose();
	working_mesh.normals() *= crot.transpose();
	working_mesh.vertices().rowwise() += ctrans.transpose();

	if (visualize_steps)
	{
		igl::opengl::glfw::Viewer viewer;
		viewer.data().set_mesh(working_mesh.vertices(), working_mesh.faces());
		viewer.launch();
	}

	// compute cusp features
	std::cout << "-- Computing cusp features...\n";
	Eigen::MatrixXd cusps;
	computeCusps(working_mesh, cusps, cuspd_params, true);

	// fit plane

	// cut teeth

	// fit curve

	// assign features to teeth

	// harmonic field stuff
}

void ToothSegmentation::computeCusps(const Mesh& mesh, Eigen::MatrixXd& features, const ToothSegmentation::CuspDetectionParams& cuspd_params, bool visualize_steps)
{
	// calculate mean curvature
	std::cout << "- Calculating mean curvature...\n";
	Eigen::SparseMatrix<double> L, M, Minv;
	igl::cotmatrix(mesh.vertices(), mesh.faces(), L);
	igl::massmatrix(mesh.vertices(), mesh.faces(), igl::MASSMATRIX_TYPE_VORONOI, M);
	igl::invert_diag(M, Minv);
	Eigen::MatrixXd mean_curvature_normals = -Minv * (L * mesh.vertices());	

	// smooth the mesh to make curvature estimate less noisy
	Eigen::MatrixXd smoothed_vertices(mesh.vertices());
	Eigen::MatrixXd smoothed_normals(mesh.normals());
	if (cuspd_params.smoothing_steps > 0)
	{
		std::cout << "- Smooothing the mesh to make the curvature estimate less noisy...\n";		
		for (std::size_t i = 0; i < cuspd_params.smoothing_steps; ++i)
		{
			std::cout << "Iteration " << i << "\n";
			smoothed_vertices.array() -= mean_curvature_normals.array() * cuspd_params.smoothing_step_size;
			igl::cotmatrix(smoothed_vertices, mesh.faces(), L);
			igl::massmatrix(smoothed_vertices, mesh.faces(), igl::MASSMATRIX_TYPE_VORONOI, M);
			igl::invert_diag(M, Minv);
			mean_curvature_normals = -Minv * (L * smoothed_vertices);			
		}
	}

	// recalculate normals
	igl::per_vertex_normals(smoothed_vertices, mesh.faces(), smoothed_normals);

	// vertices below min_feature_height are not considered
	double aabbminheight = mesh.vertices()(Eigen::all, 1).minCoeff();
	double aabbheight = mesh.vertices()(Eigen::all, 1).maxCoeff() - aabbminheight;
	Eigen::VectorXd active_vertex_map = (((mesh.vertices()(Eigen::all, 1).array() - aabbminheight) / aabbheight) - cuspd_params.min_feature_height).cwiseSign().cwiseMax(0.0).matrix();
	Eigen::VectorXi active_indices(static_cast<Eigen::DenseIndex>(active_vertex_map.sum() + 0.5));
	Eigen::DenseIndex aiidx = 0;
	for (Eigen::DenseIndex i = 0; i < active_vertex_map.rows(); ++i)
		if (active_vertex_map(i) > 0.0)
			active_indices[aiidx++] = i;

	// calculate signed mean curvatures
	Eigen::VectorXd mean_curvatures = ((mean_curvature_normals.array() * smoothed_normals.array()).matrix().rowwise().sum().array()).matrix();

	// remove outliers
	double mean_mean_curvature = mean_curvatures(active_indices.array()).mean();
	Eigen::VectorXd mcdev = (mean_mean_curvature - mean_curvatures.array()).matrix();
	double mcsdev = std::sqrt(mcdev.dot(mcdev) / static_cast<double>(mcdev.rows()));

	// remove values with zscore > max_zscore
	Eigen::VectorXd mc_zscore = ((mcdev.array() - mean_mean_curvature) / mcsdev).cwiseAbs().matrix();
	for (Eigen::DenseIndex i = 0; i < mean_curvatures.rows(); ++i)
		if (mc_zscore(i) > cuspd_params.max_zscore)
			mean_curvatures(i) = mean_mean_curvature;
	mean_mean_curvature = mean_curvatures(active_indices.array()).mean();
	for (Eigen::DenseIndex i = 0; i < mean_curvatures.rows(); ++i)
		if (mc_zscore(i) > cuspd_params.max_zscore)
			mean_curvatures(i) = mean_mean_curvature;


	std::cout << "Mean curvature standard deviation: " << mcsdev << "\n";
	std::cout << "Mean curvature mean: " << mean_mean_curvature << "\n";
	std::cout << "Mean curvature max zscore: " << mc_zscore.maxCoeff() << "\n";
	std::cout << "Mean curvature min zscore: " << mc_zscore.minCoeff() << "\n";

	if (visualize_steps)
	{
		std::cout << "maxc " << mean_curvatures.maxCoeff() << "\n minc " << mean_curvatures.minCoeff() << "\n";
		Eigen::MatrixXd C;
		igl::opengl::glfw::Viewer viewer;
		igl::jet(((mean_curvatures.array() - mean_curvatures(active_indices.array()).minCoeff()) / (mean_curvatures(active_indices.array()).maxCoeff() - mean_curvatures(active_indices.array()).minCoeff())).matrix().cwiseMax(0.0), true, C);
		//igl::jet(mc_zscore, true, C);
		viewer.data().set_mesh(mesh.vertices(), mesh.faces());
		viewer.data().set_colors(C);
		viewer.launch();
	}	

	// compute weights
	std::cout << "- Calculating feature weights...\n";
	// weighted average between negative mean curvature and y height. the result is multiplied by 0 if y height < cuspd_params.min_feature_height
	double active_min_mean_curvature = mean_curvatures(active_indices.array()).minCoeff();
	double active_max_mean_curvature = mean_curvatures(active_indices.array()).maxCoeff();
	Eigen::VectorXd curv_weights = (((mean_curvatures.array() - active_min_mean_curvature) / (active_max_mean_curvature - active_min_mean_curvature)) * active_vertex_map.array()).matrix();// +(1.0 - active_vertex_map.array()) * active_min_mean_curvature).matrix();
	double minheight = aabbminheight + aabbheight * cuspd_params.min_feature_height;
	Eigen::VectorXd height_weights = (((mesh.vertices()(Eigen::all, 1).array() - minheight) / (mesh.vertices()(Eigen::all, 1).array().maxCoeff() - minheight)) * active_vertex_map.array()).matrix();// +(1.0 - active_vertex_map.array()) * minheight).matrix();
	// pow the weights
	curv_weights = curv_weights.cwiseMax(0.0).array().pow(cuspd_params.curve_exp).matrix();
	height_weights = height_weights.cwiseMax(0.0).array().pow(cuspd_params.height_exp).matrix();
	// renormalize
	curv_weights = (curv_weights.array() - curv_weights(active_indices).minCoeff()) / (curv_weights(active_indices).maxCoeff() - curv_weights(active_indices).minCoeff());
	height_weights = (height_weights.array() - height_weights(active_indices).minCoeff()) / (height_weights(active_indices).maxCoeff() - height_weights(active_indices).minCoeff());
	// combine weights
	Eigen::VectorXd weights = ((1.0 - cuspd_params.alpha) * curv_weights.array().cwiseMax(0.0).pow(cuspd_params.curve_exp) + cuspd_params.alpha * height_weights.array().cwiseMax(0.0).pow(cuspd_params.height_exp)).matrix();

	// renormalize weights
	weights = (weights.array() - weights.minCoeff()) / (weights.maxCoeff() - weights.minCoeff());

	/*for (Eigen::DenseIndex i = 0; i < weights.rows(); ++i)
	{
		std::cout << weights(i) << "\n";
	}*/

	if (visualize_steps)
	{
		std::cout << "maxw " << weights.maxCoeff() << "\n minw " << weights.minCoeff() << "\n";
		Eigen::MatrixXd C;
		igl::jet(weights, false, C);
		igl::opengl::glfw::Viewer viewer;
		viewer.data().set_mesh(mesh.vertices(), mesh.faces());
		viewer.data().set_colors(C);
		viewer.launch();
	}

	std::cout << "- Searching local maxima...\n";

	size_t num_local_maxima = 0;
	std::vector<std::pair<Eigen::DenseIndex, double>> local_maxima;

	for (Eigen::DenseIndex v = 0; v < active_indices.rows(); ++v)
	{
		Eigen::DenseIndex i = active_indices(v);
		bool is_local_maximum = true;
		for (std::size_t j = 0; j < mesh.adjacency_list()[static_cast<std::size_t>(i)].size(); ++j)
		{
			if (weights(mesh.adjacency_list()[static_cast<std::size_t>(i)][j]) > weights(i))
			{
				is_local_maximum = false;
				break;
			}
		}
		if (is_local_maximum)
		{
			local_maxima.push_back(std::make_pair(i, weights(i)));
		}
	}	

	std::cout << "Local maxima found: " << local_maxima.size() << "\n";

	// sort local maxima by weight
	std::sort(local_maxima.begin(), local_maxima.end(), [](const std::pair<Eigen::DenseIndex, double>& a, const std::pair<Eigen::DenseIndex, double>& b) {
		return a.second > b.second;
	});

	// extract a fraction of the best local maxima for mean shift
	Eigen::DenseIndex num_filtered_maxima = std::min(std::max(static_cast<Eigen::DenseIndex>(cuspd_params.ms_frac * static_cast<double>(local_maxima.size()) + 0.5), static_cast<Eigen::DenseIndex>(1)), static_cast<Eigen::DenseIndex>(local_maxima.size()));
	Eigen::MatrixXd particles(num_filtered_maxima, 3);
	for (Eigen::DenseIndex p = 0; p < num_filtered_maxima; ++p)
	{
		particles.row(p) = mesh.vertices().row(local_maxima[p].first);
	}

	std::cout << "Local maxima considered for mean shift: " << particles.rows() << "\n";

	// to make indexing faster make a compact copy of the active vertices and weights
	Eigen::MatrixXd active_vertices(mesh.vertices()(active_indices, Eigen::all));
	Eigen::VectorXd active_weights(weights(active_indices));
	Eigen::VectorXd inverse_active_weights((1.0 - active_weights.array()).matrix());

	// build kd-tree for mean shift
	kdtree_t kdtree(3, active_vertices);
	kdtree.index->buildIndex();
	const nanoflann::SearchParams radsearchparam(32, 0.0, false);

	// mean shift
	// window size
	Eigen::Vector3d aabb_min = active_vertices.colwise().minCoeff().transpose();
	Eigen::Vector3d aabb_max = active_vertices.colwise().maxCoeff().transpose();
	double aabb_diag = (aabb_max - aabb_min).norm();
	double h = aabb_diag * cuspd_params.ms_window_size;
	double search_rad = (2.0 * h) * (2.0 * h);
	double mswvar = 2.0 * h * h;

	// accumulates total amount of shift
	double total_shift;
	double querypt[3];

	// shift vectors
	Eigen::MatrixXd shift_vectors(particles.rows(), 3);

	std::vector<std::pair<long long, double>> rad_search_res;
	Eigen::RowVector3d mean;
	double normalizer;	
	Eigen::DenseIndex nbidx;
	for (std::size_t t = 0; t < cuspd_params.ms_max_iterations; ++t)
	{
		/*if (visualize_steps)
		{
			Eigen::MatrixXd C;
			igl::jet(weights, false, C);
			igl::opengl::glfw::Viewer viewer;
			viewer.data().set_mesh(mesh.vertices(), mesh.faces());
			viewer.data().set_colors(C);
			viewer.data().point_size = 2.0;
			viewer.data().set_points(particles, Eigen::RowVector3d(1.0, 1.0, 1.0));
			viewer.launch();
		}*/
		std::cout << "Mean-shift iteration " << t + 1 << "\n";

		shift_vectors.setZero();
		for (Eigen::DenseIndex p = 0; p < particles.rows(); ++p)
		{
			rad_search_res.clear();
			querypt[0] = particles(p, 0);
			querypt[1] = particles(p, 1);
			querypt[2] = particles(p, 2);

			kdtree.index->radiusSearch(querypt, search_rad, rad_search_res, radsearchparam);

			mean.setZero();
			//normalizer = 0.0;
			normalizer = std::numeric_limits<double>::lowest();
			for (std::size_t r = 0; r < rad_search_res.size(); ++r)
			{
				nbidx = rad_search_res[r].first;
				///*auto dvec = (active_vertices(nbidx, Eigen::all) - particles(p, Eigen::all));
				//double d2 = dvec.dot(dvec);*/
				//double w = active_weights(nbidx) / static_cast<double>(rad_search_res.size());
				//mean += active_vertices(nbidx, Eigen::all) * w;
				//normalizer += w;
				if (active_weights(nbidx) > normalizer)
				{
					mean = active_vertices(nbidx, Eigen::all);
					normalizer = active_weights(nbidx);
				}
			}
			/*if(normalizer > 0.0)
				shift_vectors(p, Eigen::all) = (mean / normalizer) - particles(p, Eigen::all);
			else
				shift_vectors(p, Eigen::all).setZero();*/
			shift_vectors(p, Eigen::all) = mean - particles(p, Eigen::all);
		}

		// shift all particles
		particles.array() += shift_vectors.array();
		total_shift = shift_vectors.rowwise().norm().sum();

		if (total_shift < cuspd_params.ms_min_total_shift)
		{
			break;
		}
	}

	/*if (visualize_steps)
	{
		Eigen::MatrixXd C;
		igl::jet(weights, false, C);
		igl::opengl::glfw::Viewer viewer;
		viewer.data().set_mesh(mesh.vertices(), mesh.faces());
		viewer.data().set_colors(C);
		viewer.data().point_size = 5.0;
		viewer.data().set_points(particles, Eigen::RowVector3d(1.0, 1.0, 1.0));
		viewer.launch();
	}*/

	std::vector<bool> duplmap;
	for (Eigen::DenseIndex i = 0; i < particles.rows(); ++i)
	{
		duplmap.push_back(false);
	}

	for (Eigen::DenseIndex i = 0; i < particles.rows(); ++i)
	{
		for (Eigen::DenseIndex j = 0; j < particles.rows(); ++j)
		{
			if (i != j && !duplmap[i])
			{
				if (((particles(i, Eigen::all) - particles(j, Eigen::all)).norm() < (cuspd_params.ms_ft_collapse_dist * aabb_diag)) && (!duplmap[j]))
				{
					duplmap[j] = true;
				}
			}
		}
	}

	Eigen::DenseIndex numfeatures = 0;
	for (std::size_t i = 0; i < duplmap.size(); ++i)
	{
		if (!duplmap[i])
			numfeatures++;
	}

	features.resize(numfeatures, 3);
	Eigen::DenseIndex ftct = 0;
	for (Eigen::DenseIndex i = 0; i < static_cast<Eigen::DenseIndex>(duplmap.size()); ++i)
	{
		if (!duplmap[i])
		{
			features(ftct++, Eigen::all) = particles(i, Eigen::all);
		}
	}

	std::cout << "Number of features found: " << features.rows() << "\n";

	if (visualize_steps)
	{
		Eigen::MatrixXd C;
		igl::jet(weights, false, C);
		igl::opengl::glfw::Viewer viewer;
		viewer.data().set_mesh(mesh.vertices(), mesh.faces());
		viewer.data().set_colors(C);
		viewer.data().point_size = 7.0;
		viewer.data().set_points(features, Eigen::RowVector3d(1.0, 1.0, 1.0));
		viewer.launch();
	}
}
