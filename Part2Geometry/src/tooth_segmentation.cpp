#include "..\include\tooth_segmentation.h"
#include <iostream>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/jet.h>
#include <vector>
#include <random>
#include <Eigen/Sparse>
#include <cmath>
//#include <Eigen/SparseQR>
#include <persistence1d.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include<Eigen/IterativeLinearSolvers>

void ToothSegmentation::segmentTeethFromMesh(const Mesh& mesh, const Eigen::Vector3d& approximate_mesh_up, const Eigen::Vector3d& mesh_right, std::vector<Mesh>& tooth_meshes, const ToothSegmentation::CuspDetectionParams& cuspd_params, const HarmonicFieldParams& hf_params, const MeanCurvatureParams& mc_params, const ToothMeshExtractionParams& tme_params, bool visualize_steps)
{
	std::cout << "--- TOOTH SEGMENTATION ---\n";
	Mesh working_mesh(mesh);

	// try to find a better up vector
	std::cout << "-- Estimating up vector...\n";
	Eigen::Vector3d mesh_up = estimateUpVector(working_mesh, approximate_mesh_up);
	std::cout << "Estimated up vector: (" << mesh_up(0) << " " << mesh_up(1) << " " << mesh_up(2) << ")\n";

	if (visualize_steps)
	{
		igl::opengl::glfw::Viewer viewer;
		viewer.data().set_mesh(working_mesh.vertices(), working_mesh.faces());
		viewer.launch();
	}

	// compute mean curvature estimate
	std::cout << "-- Computing mean curvature...\n";
	Eigen::VectorXd mean_curvature(working_mesh.vertices().rows());
	computeMeanCurvature(working_mesh, mean_curvature, mc_params, visualize_steps);

	// compute cusp features
	std::cout << "-- Computing cusp features...\n";
	Eigen::VectorXi cusps;
	computeCusps(working_mesh, mesh_up, mean_curvature, cusps, cuspd_params, visualize_steps);

	// second iteration
	mesh_up = estimateUpVector(working_mesh.vertices()(cusps.array(), Eigen::all), approximate_mesh_up);
	std::cout << "Estimated up vector: (" << mesh_up(0) << " " << mesh_up(1) << " " << mesh_up(2) << ")\n";
	computeCusps(working_mesh, mesh_up, mean_curvature, cusps, cuspd_params, visualize_steps);

	// gingiva cut
	Eigen::VectorXd vertex_heights(working_mesh.vertices().rows());
	for (Eigen::Index i = 0; i < working_mesh.vertices().rows(); ++i)
		vertex_heights(i) = working_mesh.vertices()(i, Eigen::all).dot(mesh_up.transpose());

	double aabbminheight = vertex_heights.minCoeff();
	double aabbheight = vertex_heights.maxCoeff() - aabbminheight;

	std::cout << "-- Performing gingiva cut...\n";

	Eigen::VectorXi cut_indices;
	Eigen::VectorXi inverse_index_map; // maps indices from cut mesh to indices from old mesh
	Eigen::VectorXi index_map; // original mesh indices -> cut mesh indices
	cutMesh(working_mesh, cut_indices, inverse_index_map, index_map, mesh_up, (aabbminheight + aabbheight * cuspd_params.min_feature_height) * mesh_up);
	// remap mean curvature to cut mesh
	Eigen::VectorXd cut_mean_curvature(mean_curvature(inverse_index_map.array()));

	/////////////////////////////////////////////////////////
	/// SPOKE FEATURE STUFF AND AUTOMATIC GINGIVA CUTTING ///
	/////////////////////////////////////////////////////////
	auto planeresult = fitPlane(cusps, working_mesh, index_map);
	//planeresult.first is normal, planeresult.second is point on plane
	// Gingiva cut
	working_mesh.recalculateKdTree();
	auto featuregroups = segmentFeatures(cusps, mean_curvature, working_mesh, index_map);

	// assign features to teeth
	std::vector<ToothFeature> tooth_features;
	//{
	//	{{index_map(228824), index_map(20453), index_map(0), index_map(0)}, 2},
	//	{{index_map(9335), index_map(158157), index_map(165224), index_map(163144)}, 4},
	//	{{index_map(3653), index_map(155934), index_map(0), index_map(0)}, 2},
	//	{{index_map(154708), index_map(465), index_map(0), index_map(0)}, 2},
	//	{{index_map(157325), index_map(4649), index_map(0), index_map(0)}, 2},
	//	{{index_map(156244), index_map(1727), index_map(0), index_map(0)}, 2},
	//	{{index_map(154437), index_map(154323), index_map(0), index_map(0)}, 2},
	//	{{index_map(154374), index_map(281), index_map(0), index_map(0)}, 2},
	//	{{index_map(6753), index_map(158652), index_map(0), index_map(0)}, 2},
	//	{{index_map(158728), index_map(160864), index_map(0), index_map(0)}, 2},
	//	{{index_map(2179), index_map(2224), index_map(0), index_map(0)}, 2},
	//	//{{index_map(5775), index_map(160101), index_map(0), index_map(0)}, 2},
	//	{{index_map(155337), index_map(4017), index_map(5812), index_map(8910)}, 4},
	//	{{index_map(16741), index_map(174894), index_map(160439), index_map(170818)}, 4}
	//};
	for (auto fg : featuregroups)
	{
		tooth_features.push_back({ {}, Eigen::DenseIndex(fg.size()) });
		for (auto f : fg)
		{
			tooth_features.back().featurePointIndices.push_back(f);
		}
	}
	
	Eigen::Index tidx = 0;
	for (std::size_t i = 0; i < tooth_features.size(); ++i)
	{
		tidx += tooth_features[i].numFeaturePoints;
	}
	Eigen::VectorXi toothftidcs(tidx);
	tidx = 0;
	for (std::size_t i = 0; i < tooth_features.size(); ++i)
	{
		for (std::size_t t = 0; t < tooth_features[i].numFeaturePoints; ++t)
		{
			toothftidcs[tidx++] =  tooth_features[i].featurePointIndices[t];
		}
	}

	if (visualize_steps)
	{
		igl::opengl::glfw::Viewer viewer;
		viewer.data().set_mesh(working_mesh.vertices(), working_mesh.faces());
		viewer.data().set_points(working_mesh.vertices()(toothftidcs.array(), Eigen::all), Eigen::RowVector3d(1.0, 1.0, 1.0));
		viewer.data().point_size = 5.0;
		viewer.launch();
	}

	//// harmonic field stuff
	std::cout << "Solving harmonic field...\n";
	Eigen::VectorXd harmonic_field;
	calculateHarmonicField(working_mesh, cut_mean_curvature, tooth_features, cut_indices, harmonic_field, hf_params, visualize_steps);

	if (visualize_steps)
	{
		igl::opengl::glfw::Viewer viewer;
		viewer.data().set_mesh(working_mesh.vertices(), working_mesh.faces());
		viewer.data().set_points(working_mesh.vertices()(toothftidcs.array(), Eigen::all), Eigen::RowVector3d(1.0, 1.0, 1.0));
		viewer.data().point_size = 5.0;
		Eigen::MatrixXd C(working_mesh.vertices().rows(), 3);
		igl::jet(harmonic_field, true, C);
		viewer.data().set_colors(C);
		viewer.launch();
	}

	// extract tooth meshes and return
	tooth_meshes = extractToothMeshes(working_mesh, harmonic_field, tooth_features, tme_params, visualize_steps);
}

void ToothSegmentation::computeMeanCurvature(const Mesh & mesh, Eigen::VectorXd & mean_curvature, const MeanCurvatureParams & mc_params, bool visualize_steps)
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
	if (mc_params.smoothing_steps > 0)
	{
		std::cout << "- Smooothing the mesh to make the curvature estimate less noisy...\n";
		for (std::size_t i = 0; i < mc_params.smoothing_steps; ++i)
		{
			std::cout << "Iteration " << i << "\n";
			smoothed_vertices.array() -= mean_curvature_normals.array() * mc_params.smoothing_step_size;
			igl::cotmatrix(smoothed_vertices, mesh.faces(), L);
			igl::massmatrix(smoothed_vertices, mesh.faces(), igl::MASSMATRIX_TYPE_VORONOI, M);
			igl::invert_diag(M, Minv);
			mean_curvature_normals = -Minv * (L * smoothed_vertices);
		}
	}

	// recalculate normals
	igl::per_vertex_normals(smoothed_vertices, mesh.faces(), smoothed_normals);

	// calculate signed mean curvatures
	Eigen::VectorXd mean_curvatures = ((mean_curvature_normals.array() * smoothed_normals.array()).matrix().rowwise().sum().array()).matrix();
	mean_curvature.resize(mean_curvatures.rows());
	mean_curvature = mean_curvatures;

	// remove outliers
	double mean_mean_curvature = mean_curvatures.mean();
	Eigen::VectorXd mcdev = (mean_mean_curvature - mean_curvatures.array()).matrix();
	double mcsdev = std::sqrt(mcdev.dot(mcdev) / static_cast<double>(mcdev.rows()));

	// remove values with zscore > max_zscore
	Eigen::VectorXd mc_zscore = ((mcdev.array() - mean_mean_curvature) / mcsdev).cwiseAbs().matrix();
	for (Eigen::DenseIndex i = 0; i < mean_curvatures.rows(); ++i)
		if (mc_zscore(i) > mc_params.max_zscore)
			mean_curvatures(i) = mean_mean_curvature;
	mean_mean_curvature = mean_curvatures.mean();
	for (Eigen::DenseIndex i = 0; i < mean_curvatures.rows(); ++i)
		if (mc_zscore(i) > mc_params.max_zscore)
			mean_curvatures(i) = mean_mean_curvature;

	mean_curvature = mean_curvatures;
	std::cout << "Mean curvature standard deviation: " << mcsdev << "\n";
	std::cout << "Mean curvature mean: " << mean_mean_curvature << "\n";
	std::cout << "Mean curvature max zscore: " << mc_zscore.maxCoeff() << "\n";
	std::cout << "Mean curvature min zscore: " << mc_zscore.minCoeff() << "\n";

	if (visualize_steps)
	{
		Eigen::MatrixXd C;
		igl::opengl::glfw::Viewer viewer;
		igl::jet(((mean_curvatures.array() - mean_curvature.minCoeff()) / (mean_curvatures.maxCoeff() - mean_curvatures.minCoeff())).matrix().cwiseMax(0.0), true, C);
		//igl::jet(mc_zscore, true, C);
		viewer.data().set_mesh(mesh.vertices(), mesh.faces());
		viewer.data().set_colors(C);
		viewer.launch();
	}
}

void ToothSegmentation::computeCusps(const Mesh & mesh, const Eigen::Vector3d & mesh_up, const Eigen::VectorXd & mean_curvature, Eigen::VectorXi & features, const CuspDetectionParams & cuspd_params, bool visualize_steps)
{
	// vertices below min_feature_height are not considered
	Eigen::VectorXd vertex_heights(mesh.vertices().rows());
	for (Eigen::Index i = 0; i < mesh.vertices().rows(); ++i)
		vertex_heights(i) = mesh.vertices()(i, Eigen::all).dot(mesh_up.transpose());

	double aabbminheight = vertex_heights.minCoeff();
	double aabbheight = vertex_heights.maxCoeff() - aabbminheight;
	Eigen::VectorXd active_vertex_map = (((vertex_heights.array() - aabbminheight) / aabbheight) - cuspd_params.min_feature_height).cwiseSign().cwiseMax(0.0).matrix();
	Eigen::VectorXi active_indices(static_cast<Eigen::DenseIndex>(active_vertex_map.sum() + 0.5));
	Eigen::DenseIndex aiidx = 0;
	for (Eigen::DenseIndex i = 0; i < active_vertex_map.rows(); ++i)
		if (active_vertex_map(i) > 0.0)
			active_indices[aiidx++] = i;

	// compute weights
	std::cout << "- Calculating feature weights...\n";
	// weighted average between negative mean curvature and y height. the result is multiplied by 0 if y height < cuspd_params.min_feature_height
	double active_min_mean_curvature = mean_curvature(active_indices.array()).minCoeff();
	double active_max_mean_curvature = mean_curvature(active_indices.array()).maxCoeff();
	Eigen::VectorXd curv_weights = (((mean_curvature.array() - active_min_mean_curvature) / (active_max_mean_curvature - active_min_mean_curvature)) * active_vertex_map.array()).matrix();
	double minheight = aabbminheight + aabbheight * cuspd_params.min_feature_height;
	Eigen::VectorXd height_weights = (((vertex_heights.array() - minheight) / (vertex_heights.maxCoeff() - minheight)) * active_vertex_map.array()).matrix();
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

	// non-maximum suppression
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

	// extract a fraction of the best local maxima for optimum shift
	Eigen::DenseIndex num_filtered_maxima = std::min(std::max(static_cast<Eigen::DenseIndex>(cuspd_params.os_frac * static_cast<double>(local_maxima.size()) + 0.5), static_cast<Eigen::DenseIndex>(1)), static_cast<Eigen::DenseIndex>(local_maxima.size()));
	Eigen::MatrixXd particles(num_filtered_maxima, 3);
	for (Eigen::DenseIndex p = 0; p < num_filtered_maxima; ++p)
	{
		particles.row(p) = mesh.vertices().row(local_maxima[p].first);
	}

	std::cout << "Local maxima considered for optimum shift: " << particles.rows() << "\n";

	// to make indexing faster make a compact copy of the active vertices and weights
	Eigen::MatrixXd active_vertices(mesh.vertices()(active_indices, Eigen::all));
	Eigen::VectorXd active_weights(weights(active_indices));
	Eigen::VectorXd inverse_active_weights((1.0 - active_weights.array()).matrix());

	// build kd-tree for optimum shift
	kdtree_t kdtree(3, active_vertices);
	kdtree.index->buildIndex();
	const nanoflann::SearchParams radsearchparam(32, 0.0, false);

	// optimum shift
	// window size
	Eigen::Vector3d aabb_min = active_vertices.colwise().minCoeff().transpose();
	Eigen::Vector3d aabb_max = active_vertices.colwise().maxCoeff().transpose();
	double aabb_diag = (aabb_max - aabb_min).norm();
	double h = aabb_diag * cuspd_params.os_window_size;
	double search_rad = (2.0 * h) * (2.0 * h);
	//double mswvar = 2.0 * h * h;

	// accumulates total amount of shift
	double total_shift;
	double querypt[3];

	// shift vectors
	Eigen::MatrixXd shift_vectors(particles.rows(), 3);

	std::vector<std::pair<long long, double>> rad_search_res;
	Eigen::RowVector3d mean;
	double normalizer;
	Eigen::DenseIndex nbidx;
	for (std::size_t t = 0; t < cuspd_params.os_max_iterations; ++t)
	{
		std::cout << "Optimum-shift iteration " << t + 1 << "\n";

		shift_vectors.setZero();
		for (Eigen::DenseIndex p = 0; p < particles.rows(); ++p)
		{
			rad_search_res.clear();
			querypt[0] = particles(p, 0);
			querypt[1] = particles(p, 1);
			querypt[2] = particles(p, 2);

			kdtree.index->radiusSearch(querypt, search_rad, rad_search_res, radsearchparam);

			mean.setZero();
			normalizer = std::numeric_limits<double>::lowest();
			for (std::size_t r = 0; r < rad_search_res.size(); ++r)
			{
				nbidx = rad_search_res[r].first;
				if (active_weights(nbidx) > normalizer)
				{
					mean = active_vertices(nbidx, Eigen::all);
					normalizer = active_weights(nbidx);
				}
			}
			shift_vectors(p, Eigen::all) = mean - particles(p, Eigen::all);
		}

		// shift all particles
		particles.array() += shift_vectors.array();
		total_shift = shift_vectors.rowwise().norm().sum();

		if (total_shift < cuspd_params.os_min_total_shift)
		{
			break;
		}
	}

	// collapse features with distance < threshold
	std::cout << "Merging features...\n";
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
				if (((particles(i, Eigen::all) - particles(j, Eigen::all)).norm() < (cuspd_params.ft_collapse_dist * aabb_diag)) && (!duplmap[j]))
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

	features.resize(numfeatures);
	Eigen::DenseIndex ftct = 0;
	for (Eigen::DenseIndex i = 0; i < static_cast<Eigen::DenseIndex>(duplmap.size()); ++i)
	{
		if (!duplmap[i])
		{
			double qp[] = { particles(i, 0), particles(i, 1), particles(i, 2) };
			Eigen::Index closestidx;
			double dist;
			mesh.kdtree().query(qp, 1, &closestidx, &dist);
			features(ftct++) = closestidx;
		}
	}

	// remove features with low local neighborhood (essentially removes too small features)
	std::cout << "Removing small features...\n";
	double sftr_window_size = aabb_diag * cuspd_params.small_ft_window_size;
	search_rad = sftr_window_size * sftr_window_size;
	double ft_mean = weights(features.array()).mean();
	std::vector<Eigen::DenseIndex> final_features;
	//Eigen::VectorXd sizeweight(features.rows());
	for (Eigen::DenseIndex i = 0; i < features.rows(); ++i)
	{
		rad_search_res.clear();
		querypt[0] = mesh.vertices()(features(i), 0);
		querypt[1] = mesh.vertices()(features(i), 1);
		querypt[2] = mesh.vertices()(features(i), 2);

		mesh.kdtree().index->radiusSearch(querypt, search_rad, rad_search_res, radsearchparam);

		double mean_nb_score = 0.0;
		for (size_t s = 0; s < rad_search_res.size(); ++s)
		{
			if (rad_search_res[s].first != features(i))
				mean_nb_score += weights(rad_search_res[s].first);
		}
		mean_nb_score /= std::max(rad_search_res.size() - 1, 1ull);

		// if mean neighborhood score is less than center feature score - threshold, remove this feature
		if (ft_mean - mean_nb_score <= cuspd_params.small_ft_threshold)
			final_features.push_back(features(i));

		//sizeweight(i) = weights(features(i)) - mean_nb_score;
	}
	features.resize(final_features.size());
	for (size_t i = 0; i < final_features.size(); ++i)
		features(i) = final_features[i];

	
	std::cout << "Number of features found: " << features.rows() << "\n";

	if (visualize_steps)
	{
		Eigen::MatrixXd C;
		//Eigen::MatrixXd CF;
		igl::jet(weights, false, C);
		//igl::jet(sizeweight, true, CF);
		igl::opengl::glfw::Viewer viewer;
		viewer.data().set_mesh(mesh.vertices(), mesh.faces());
		viewer.data().set_colors(C);
		viewer.data().point_size = 7.0;
		viewer.data().set_points(mesh.vertices()(features.array(), Eigen::all), Eigen::RowVector3d{ 1.0, 1.0, 1.0 });
		viewer.launch();
	}
}

void ToothSegmentation::calculateHarmonicField(const Mesh& mesh, const Eigen::VectorXd& mean_curvature, const std::vector<ToothFeature>& toothFeatures, const Eigen::VectorXi& cutIndices, Eigen::VectorXd& harmonic_field, const HarmonicFieldParams& hf_params, bool visualize_steps)
{
	std::cout << "Building laplacian matrix...\n";
	// sort indices by constraint type
	// tooth features
	Eigen::Index num_tooth_features = 0;
	for (const auto& t : toothFeatures)
		num_tooth_features += t.numFeaturePoints;

	// calculate laplacian matrix
	Eigen::SparseMatrix<double> L(mesh.vertices().rows(), mesh.vertices().rows());
	L.setZero();
	Eigen::VectorXd b(mesh.vertices().rows());	
	b.setZero();

	std::vector<Eigen::Triplet<double>> Ltripls;

	double max_curvature = mean_curvature.maxCoeff();
	double min_curvature = mean_curvature.minCoeff();
	for (Eigen::Index i = 0; i < mesh.vertices().rows(); ++i)
	{
		double iweight = 0.0;
		for (Eigen::Index a = 0; a < mesh.adjacency_list()[i].size(); ++a)
		{
			Eigen::Index j = mesh.adjacency_list()[i][a];			
			double cotij = calcCotanWeight(i, j, mesh) * calcCurvatureWeight(i, j, mean_curvature, hf_params);
			Ltripls.push_back(Eigen::Triplet<double>(i, j, -cotij));
			iweight += cotij;
		}
		Ltripls.push_back(Eigen::Triplet<double>(i, i, iweight));
	}
	std::cout << "L non-zero entries: " << Ltripls.size() << "\n";
	std::cout << "Num vertices: " << mesh.vertices().rows() << "\n";

	L.setFromTriplets(Ltripls.begin(), Ltripls.end());

	// mass matrix to account for triangulation
	Eigen::SparseMatrix<double> M, Minv;
	igl::massmatrix(mesh.vertices(), mesh.faces(), igl::MASSMATRIX_TYPE_VORONOI, M);
	igl::invert_diag(M, Minv);
	L = Minv * L;

	std::cout << "Applying constraints...\n";
	// tooth features
	std::cout << "Tooth features...\n";
	for (Eigen::Index t = 0; t < toothFeatures.size(); ++t)
	{
		for (Eigen::Index i = 0; i < toothFeatures[t].numFeaturePoints; ++i)
		{
			L.row(toothFeatures[t].featurePointIndices[i]) *= 0.0;
			L.coeffRef(toothFeatures[t].featurePointIndices[i], toothFeatures[t].featurePointIndices[i]) = 1.0;
			b(toothFeatures[t].featurePointIndices[i]) = ((t % 2 == 0) ? 0.0 : hf_params.w);
		}
	}

	// cut boundary vertices
	std::cout << "Cut boundary...\n";
	for (Eigen::Index i = 0; i < cutIndices.rows(); ++i)
	{
		L.row(cutIndices(i)) *= 0.0;
		L.coeffRef(cutIndices(i), cutIndices(i)) = 1.0;
		b(cutIndices(i)) = 0.5 * hf_params.w;
	}

	// solve sparse system
	std::cout << "Solving linear system...\n";
	//Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
	solver.compute(L);
	if (solver.info() != Eigen::Success)
	{
		std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SPARSE QR DECOMPOSTION FAILED !!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
		harmonic_field.resize(mesh.vertices().rows());
		harmonic_field.setZero();
	}
	harmonic_field = solver.solve(b);
}

void ToothSegmentation::cutMesh(Mesh& mesh, Eigen::VectorXi& cut_indices, Eigen::VectorXi& ivrs_index_map, Eigen::VectorXi& _index_map, const Eigen::Vector3d& normal, const Eigen::Vector3d& plane_point)
{
	std::vector<Eigen::RowVector3d> newvertices;
	std::vector<Eigen::RowVector3i> newfaces;
	std::vector<Eigen::DenseIndex> cutindices_old;

	Eigen::VectorXi index_map(mesh.vertices().rows());
	index_map.setConstant(-1);
	
	Eigen::VectorXd vertex_plane_distances = (mesh.vertices().rowwise() - plane_point.transpose()) * normal;

	for (Eigen::DenseIndex f = 0; f < mesh.faces().rows(); ++f)
	{
		if (vertex_plane_distances(mesh.faces()(f, 0)) > 0.0 && vertex_plane_distances(mesh.faces()(f, 1)) > 0.0 && vertex_plane_distances(mesh.faces()(f, 2)) > 0.0)
		{
			Eigen::RowVector3i newface;
			for (Eigen::DenseIndex i = 0; i < 3; ++i)
			{
				if (index_map(mesh.faces()(f, i)) == -1) // if vertex is not yet in new set, add it and put the correponding index into index map
				{
					newvertices.push_back(mesh.vertices()(mesh.faces()(f, i), Eigen::all));
					index_map(mesh.faces()(f, i)) = newvertices.size() - 1;
					newface(i) = newvertices.size() - 1;
				}
				else // append existing new index to face
				{
					newface(i) = index_map(mesh.faces()(f, i));
				}
			}
			newfaces.push_back(newface);
		}
		else
		{
			if (vertex_plane_distances(mesh.faces()(f, 0)) > 0.0)
			{
				cutindices_old.push_back(mesh.faces()(f, 0));
			}
			if (vertex_plane_distances(mesh.faces()(f, 1)) > 0.0)
			{
				cutindices_old.push_back(mesh.faces()(f, 1));
			}
			if (vertex_plane_distances(mesh.faces()(f, 2)) > 0.0)
			{
				cutindices_old.push_back(mesh.faces()(f, 2));
			}
		}		
	}

	Eigen::MatrixXd Vnew(newvertices.size(), 3);
	for (std::size_t i = 0; i < newvertices.size(); ++i)
		Vnew(i, Eigen::all) = newvertices[i];

	Eigen::MatrixXi Fnew(newfaces.size(), 3);
	for (std::size_t i = 0; i < newfaces.size(); ++i)
		Fnew(i, Eigen::all) = newfaces[i];

	Eigen::MatrixXd Nnew;
	igl::per_vertex_normals(Vnew, Fnew, Nnew);

	mesh = Mesh(Vnew, Nnew, Fnew);

	Eigen::DenseIndex num_cut_indices = 0;
	for (Eigen::DenseIndex i = 0; i < cutindices_old.size(); ++i)
		if (index_map(cutindices_old[i]) != -1)
			num_cut_indices++;

	cut_indices.resize(num_cut_indices);
	num_cut_indices = 0;
	for (Eigen::DenseIndex i = 0; i < cutindices_old.size(); ++i)
		if (index_map(cutindices_old[i]) != -1)
			cut_indices[num_cut_indices++] = index_map(cutindices_old[i]);

	// construct inverse index map (for handling old arrays defined over the mesh)
	Eigen::DenseIndex ivrsidxct = 0;
	ivrs_index_map.resize(Vnew.rows());
	for (Eigen::DenseIndex i = 0; i < index_map.rows(); ++i)
		if(index_map(i) != -1)
			ivrs_index_map(index_map(i)) = i;

	_index_map = index_map;
}

double ToothSegmentation::calcCotanWeight(const Eigen::Index & i, const Eigen::Index & j, const Mesh & mesh)
{
	bool is_adjacent = false;
	for (std::size_t n = 0; n < mesh.adjacency_list()[i].size(); ++n)
	{
		if (mesh.adjacency_list()[i][n] == j)
		{
			is_adjacent = true;
		}
	}

	if (is_adjacent)
	{
		Eigen::Vector3d vi = mesh.vertices()(i, Eigen::all);
		Eigen::Vector3d vj = mesh.vertices()(j, Eigen::all);
		double res = 0.0;
		for (Eigen::Index t = 0; t < mesh.triangle_list()[i].size(); ++t)
		{
			for (Eigen::Index k = 0; k < 3; ++k)
			{
				if (mesh.faces()(mesh.triangle_list()[i][t], k) == j)
				{
					for (Eigen::Index l = 0; l < 3; ++l)
					{
						if (mesh.faces()(mesh.triangle_list()[i][t], l) != i && mesh.faces()(mesh.triangle_list()[i][t], l) != j)
						{
							Eigen::Vector3d vl = mesh.vertices()(mesh.faces()(mesh.triangle_list()[i][t], l), Eigen::all);
							res += ((vi - vl).normalized().dot((vj - vl).normalized())) / ((vi - vl).normalized().cross((vj - vl).normalized())).norm();
						}
					}
				}
			}
		}
		return res / 2.0;
	}
	else
	{
		return 0.0;
	}
}

double ToothSegmentation::calcCurvatureWeight(const Eigen::Index & i, const Eigen::Index & j, const Eigen::VectorXd & mean_curvature, const HarmonicFieldParams & hf_params)
{
	double nmci = std::abs(std::min(mean_curvature(i), 0.0));
	double nmcj = std::abs(std::min(mean_curvature(j), 0.0));
	if (nmci >= hf_params.concavity_threshold || nmcj >= hf_params.concavity_threshold)
		return hf_params.gamma_low;
	else
		return hf_params.gamma_high;
}

Eigen::Vector3d ToothSegmentation::estimateUpVector(const Mesh & mesh, const Eigen::Vector3d & approximate_up)
{
	Eigen::MatrixXd zero_mean_data(mesh.vertices().rowwise() - mesh.vertices().colwise().mean());
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(zero_mean_data, Eigen::ComputeThinU | Eigen::ComputeThinV);
	return (svd.matrixV().col(2) * (svd.matrixV().col(2).dot(approximate_up) >= 0.0 ? 1.0 : -1.0)).normalized();
}

Eigen::Vector3d ToothSegmentation::estimateUpVector(const Eigen::MatrixXd& points, const Eigen::Vector3d & approximate_up)
{
	Eigen::MatrixXd zero_mean_data(points.rowwise() - points.colwise().mean());
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(zero_mean_data, Eigen::ComputeThinU | Eigen::ComputeThinV);
	return (svd.matrixV().col(2) * (svd.matrixV().col(2).dot(approximate_up) >= 0.0 ? 1.0 : -1.0)).normalized();
}
std::pair<Eigen::Vector3d, Eigen::Vector3d> ToothSegmentation::fitPlane(const Eigen::VectorXi& featureindices, const Mesh& mesh, Eigen::VectorXi idmap)
{
	Eigen::MatrixXd features(featureindices.rows(), 3);// = mesh.vertices();
	for (size_t idx = 0; idx < featureindices.size(); ++idx)
	{
		std::cout << featureindices(idx) << std::endl;
		features.row(idx) = mesh.vertices().row(idmap(featureindices(idx)));
	}
	Eigen::Vector3d center = features.colwise().mean();
	auto normed = Eigen::MatrixXd(features);
	normed.rowwise() -= center.transpose();
	auto svd = normed.bdcSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
	std::cout << "V: \n" << svd.matrixV().transpose() << std::endl;
	return { svd.matrixV().transpose().bottomRows(1).row(0) , center };
	
}

std::vector<std::vector<size_t>> ToothSegmentation::segmentFeatures(const Eigen::VectorXi& featureindices, Eigen::VectorXd meancurvature, const Mesh& mesh, Eigen::VectorXi idmap)
{
	auto Vs = mesh.vertices();
	auto Fs = mesh.faces();
	auto Ns = mesh.normals();
	Vs *= Eigen::MatrixXd(Eigen::AngleAxis<double>(M_PI, Eigen::Vector3d{ 0,0,1 }).toRotationMatrix());
	Mesh newm = Mesh(Vs, Ns, Fs);
	newm.recalculateKdTree();
	//############ Calculating curve, 3rd degreee poly
	Eigen::MatrixXd M(4, 4);
	M.fill(0);
	Eigen::MatrixXd polied(featureindices.rows(), 4);
	//Construct beta
	Eigen::VectorXd b_3(4);
	auto wf = [](double x, double xm) {return 1 / (std::abs(x - xm) + 1); };
	double x_middle = 0;
	double z_max = -999999;
	Eigen::MatrixXd features(featureindices.rows(), 3);// = mesh.vertices();
	for (size_t idx=0; idx<featureindices.size();++idx)
	{
		features.row(idx) = newm.vertices().row(idmap(featureindices(idx)));
	}
	for (auto& row : features.rowwise())
	{
		if (row(2) > z_max)
		{
			z_max = row(2);
			x_middle = row(0);
		}

	}
	for (size_t x = 0; x < polied.rows(); ++x)
	{
		polied.row(x) << std::pow(features(x, 0), 3), std::pow(features(x, 0), 2), std::pow(features(x, 0), 1), 1; //std::pow(features(x, 0), 4),
		M += Eigen::MatrixXd(polied.row(x).transpose() * polied.row(x));
		//std::cout << "polied row: " << polied.row(x).transpose() * polied.row(x) << std::endl;
		//polied.row(x) *= wf(features(x, 0), x_middle);
	}

	//M = polied.colwise().sum().transpose() * polied.colwise().sum();
	std::cout << "polied: " << polied << std::endl;
	std::cout << "M: \n" << M << std::endl;
	for (size_t x = 0; x < 4; ++x)
	{
		double v = 0;
		for (auto& xyz : features.rowwise())
		{
			v += std::pow(xyz(0), 3 - x) * xyz(2); //wf(xyz(0), x_middle) *
		}

		b_3.row(x) << v;
	}
	Eigen::MatrixXd curveparams = M.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_3);

	//################ Calculate spokes
	double min_x, max_x;
	min_x = features.colwise().minCoeff()(0);
	max_x = features.colwise().maxCoeff()(0);
	std::cout << "min, max x: " << min_x << " | " << max_x << std::endl;
	std::vector<std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>> spokes{};
	int num_spokes = 100;
	spokes.reserve(100);
	std::vector<Eigen::Vector3d> curvepoints;
	double fpymax = features.colwise().maxCoeff()(1);
	double verticesminy = newm.vertices().colwise().minCoeff()(1);
	std::cout << "verticesminy: " << verticesminy << std::endl;
	double stepsize = std::abs(verticesminy - newm.vertices().colwise().maxCoeff()(1)) / 200.0;
	std::cout << "fpymax: " << fpymax << std::endl;
	double curvelength = (max_x - min_x);
	double spokeinterdistance = curvelength / (14*10.0);
	for (double x = min_x - spokeinterdistance; x <= max_x + spokeinterdistance; x += spokeinterdistance) //+(max_x-min_x)/25
	{

		double z = curveparams(0, 0) * std::pow(x, 3) + curveparams(1, 0) * std::pow(x, 2) + curveparams(2, 0) * std::pow(x, 1) + curveparams(3, 0) * std::pow(x, 0);// +curveparams(4, 0);
		double dz = curveparams(0, 0) * std::pow(x, 2) * 3 + curveparams(1, 0) * std::pow(x, 1) * 2 + curveparams(2, 0) * std::pow(x, 0) * 1;// +curveparams(3, 0);
		curvepoints.push_back({ x,fpymax,z });
		//z = dz*z+b
		Eigen::Vector2d normal{ 1, -(1 / dz) };
		normal.normalize();
		spokes.push_back({});
		for (int theta = -3; theta <= 3; ++theta)
		{
			double t = theta * 5 * (M_PI / 180.0);
			Eigen::Matrix2d rot;
			rot << std::cos(t), -std::sin(t), std::sin(t), std::cos(t);
			//std::cout << "norm: " << normal << std::endl;
			auto rotnorm = rot * normal;
			//std::cout << "rotated: " << rotnorm << std::endl;
			spokes.back().push_back({ { rotnorm(0), 0, rotnorm(1) }, { x, fpymax, z } });
		}
	}

	// For every spoke, create samples on spoke and find find closest vertice to plane spun up by thhe normal vector and up vector OR vertices in radius of point, store y value.
	// Pick highest y value along spoke -> among all spokes pick the one with the lowest max
	std::vector<std::vector<double>> depths;
	depths.reserve(100);
	std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> finalspokes;
	std::vector<double> finalspokecurvatures;

	std::cout << "stepsize: " << stepsize << std::endl;
	for (size_t outer = 0; outer < spokes.size(); outer += 1)
	{
		double minspokey = std::numeric_limits<double>::max();
		long minspokeIdx = -1;
		Eigen::Index minspokepointid;
		double meancurvaturealongspoke;
		for (size_t inner = 0; inner < spokes[outer].size(); inner++)
		{
			double maxy = std::numeric_limits<double>::lowest();
			Eigen::Vector3d loc = spokes[outer][inner].second;
			std::vector<double> curves;
			for (double x = -1.5; x <= 1.5; x += 0.1)
			{
				Eigen::Vector3d pos = loc + (spokes[outer][inner].first * x * 5);
				//Let's go down step by step in negative y direction until we find a close vertice
				std::vector<std::pair<Eigen::Index, double>> srchres;
				for (double y = pos(1); y >= verticesminy; y -= stepsize)
				{
					auto neighbor = newm.kdtree().index->radiusSearch(pos.data(), stepsize * stepsize, srchres, {});
					if (srchres.size() > 0)
					{
						auto r = newm.vertices()(srchres.front().first, 1);
						curves.push_back(meancurvature[srchres.front().first] );
						//std::cout << " srchres size: " << srchres.size() << "y: " << pos << " id: "<<spokes[outer][inner].second<< std::endl;

						if (r > maxy)// && meancurvature[srchres.front().first] <= 0)
						{
							maxy = r;
							spokes[outer][inner].second(1) = r;
							minspokepointid = srchres.front().first;
							//curvepoints.push_back(mesh.vertices().row(srchres.front().first));
						//	curvepoints.push_back(pos);
						}
						break;
					}

					pos(1) -= stepsize;
				}
			}
			if (maxy < minspokey)
			{
				minspokey = maxy;
				minspokeIdx = inner;
				for (auto& x : curves)
				{
					meancurvaturealongspoke += x;
				}
				meancurvaturealongspoke /= curves.size();
			}
		}
		if (minspokeIdx != -1 )
		{
			finalspokes.push_back({ spokes[outer][minspokeIdx] });
			finalspokecurvatures.push_back(meancurvaturealongspoke);
			std::cout << "Spoke point curvature: " << meancurvaturealongspoke << std::endl;
		}
		else
		{
			std::cout << "spoke unsolved: " << outer << std::endl;
		}
	}
	//Find minima
	std::vector<float> spokeheights;
	spokeheights.reserve(finalspokes.size());
	for (auto& row : finalspokes)
	{
		spokeheights.push_back(row.second(1));
	}
	p1d::Persistence1D pers1d;
	pers1d.RunPersistence(spokeheights);
	std::vector<p1d::TPairedExtrema> extrema;
	pers1d.GetPairedExtrema(extrema, 0.7);
	std::sort(extrema.begin(), extrema.end(), [](p1d::TPairedExtrema& a, p1d::TPairedExtrema& b) {return a.MinIndex < b.MinIndex; });
	
	Eigen::MatrixXd spokepoints(30 * finalspokes.size() + 1, 3);
	spokepoints.fill(0);
	for (size_t outer = 0; outer < finalspokes.size(); outer += 1)
	{
		int x1 = 0;
		for (double x = -1; x <= 1; x += 0.1)
		{
			Eigen::Vector3d pos = finalspokes[outer].second + (finalspokes[outer].first * x*5); //finalspokes[outer].second + (finalspokes[outer].first * x*5); //
			//spokepoints.row(outer * 20 + x1) = pos;
			x1++;
		}
	}
	//auto extr2(extrema);
	//extr2.clear();
	//for (auto& x : extrema)
	//{
	//	std::vector<std::pair<Eigen::Index, double>> srchres;
	//	auto neighbor = mesh.kdtree().query()
	//	if()
	//}
	// Remove spokes which have no feature points in between one another
	// FOr every 2 spokes, check on what side each feature point is, in relation to those 2 spokes. 
	double memeavg = 0;
	std::vector<size_t> extremaindices;
	extremaindices.push_back(0);
	for (size_t s = 0; s < extrema.size(); s += 1)
	{
		std::cout << "Persistence: " << extrema[s].Persistence << std::endl;
		//if (finalspokecurvatures[extrema[s].MinIndex] > 0)
		//{
		//	continue;
		//}
		extremaindices.push_back(extrema[s].MinIndex);
	}
	extremaindices.push_back(finalspokes.size() - 1);

	for (auto& id : extremaindices)
	{
		std::cout << "spoke curvature: " << finalspokecurvatures[id] << std::endl;
		memeavg += finalspokecurvatures[id];
		int x1 = 0;
		for (double x = -1.5; x <= 1.5; x += 0.1)
		{
			Eigen::Vector3d pos = finalspokes[id].second + (finalspokes[id].first * x * 5); //finalspokes[outer].second + (finalspokes[outer].first * x*5); //
			curvepoints.push_back(pos);
			x1++;
		}
	}



	std::vector<std::vector<size_t>> featuregroups;

	for (size_t s = 0; s < extremaindices.size() - 1; s += 1)
	{
		Eigen::Vector3d spoke1_leveled, spoke2_leveled;
		Eigen::Vector3d comSpokes = (finalspokes[extremaindices[s]].second + finalspokes[extremaindices[s+1]].second) / 2.0;
		double featuregroupdistThreshold = (comSpokes - features.colwise().mean().transpose()).norm();
		int xp = 0;
		std::vector<Eigen::Vector3d> spokps;
		std::vector<Eigen::Vector3d> grpp;
		std::vector<std::pair<Eigen::Vector3d, size_t>> group1;
		std::vector<std::pair<Eigen::Vector3d, size_t>> group2;

		for (size_t idx = 0; idx<features.rows();++idx)
		{
			Eigen::Vector3d feature_leveled;
			auto fp = features.row(idx);
			feature_leveled << fp(0), 0, fp(2);

			spoke1_leveled << finalspokes[extremaindices[s]].second(0), 0, finalspokes[extremaindices[s]].second(2);
			spoke2_leveled << finalspokes[extremaindices[s+1]].second(0), 0, finalspokes[extremaindices[s+1]].second(2);
			Eigen::Vector3d pos = feature_leveled - spoke1_leveled; //finalspokes[outer].second + (finalspokes[outer].first * x*5); //
			Eigen::Vector3d pos2 = feature_leveled - spoke2_leveled;
			Eigen::Vector3d norm1, norm2;
			norm1 = finalspokes[extremaindices[s]].first;
			norm2 = finalspokes[extremaindices[s + 1]].first;
			
			if (std::signbit(norm1.cross(pos)(1)) != std::signbit(norm2.cross(pos2)(1)))
			{
				group1.push_back({ fp.transpose() + Eigen::Vector3d{ 0, 1, 0 }, idx});
				grpp.push_back(group1.back().first);
			}
			else
			{
				group2.push_back({ fp.transpose() + Eigen::Vector3d{ 0, 1, 0 }, idx });
			}
		}


		featuregroups.push_back({});
		if (group1.size() > group2.size())
		{
			group1 = group2;
		}
		for (auto& fp : group1)
		{
			double distance = (comSpokes - fp.first).norm();
			if (distance > featuregroupdistThreshold)
			{
				continue;
			}
			featuregroups.back().push_back(idmap(featureindices(fp.second)));
		}

		int x1 = 0;
		for (double x = -1.5; x <= 1.5; x += 0.1)
		{
			Eigen::Vector3d pos = finalspokes[extremaindices[s]].second + (finalspokes[extremaindices[s]].first * x * 5); //finalspokes[outer].second + (finalspokes[outer].first * x*5); //
			spokps.push_back(pos);
			pos = finalspokes[extremaindices[s+1]].second + (finalspokes[extremaindices[s+1]].first * x * 5); //finalspokes[outer].second + (finalspokes[outer].first * x*5); //
			spokps.push_back(pos);
			x1++;
		}

		Eigen::MatrixXd pb1(spokps.size(), 3);
		Eigen::MatrixXd pb2(grpp.size(), 3);
		for (int x = 0; x < spokps.size(); ++x)
		{
			pb1.row(x) = spokps[x];
		}
		for (int x = 0; x < grpp.size(); ++x)
		{
			pb2.row(x) = grpp[x];
		}

		//igl::opengl::glfw::Viewer viewer2;
		//viewer2.data().set_mesh(newm.vertices(), newm.faces());
		//viewer2.data().set_colors(Eigen::RowVector3d(1, 1, 0.5));// Eigen::RowVector3d(1, 0, 1));
		//viewer2.data().add_points(features, Eigen::RowVector3d(1, 0, 0));
		//viewer2.data().add_points(pb1, Eigen::RowVector3d(0, 1, 1));
		//viewer2.data().add_points(pb2, Eigen::RowVector3d(0, 1, 0));

		//viewer2.data().point_size = 10;
		//viewer2.launch();


	}
	std::cout << "spoke curvature avg: " << memeavg/double(extrema.size()) << std::endl;
	std::vector<Eigen::MatrixXd> results;
	results.push_back(Eigen::MatrixXd{ curvepoints.size(), 3 });
	int m = spokes.size() * spokes[0].size();
	for (size_t x = 0; x < curvepoints.size(); x++)
	{
		results.back().row(x) = curvepoints[x];
	}
	results.push_back(spokepoints);
	igl::opengl::glfw::Viewer viewer2;
	viewer2.data().set_mesh(newm.vertices(), newm.faces());
	viewer2.data().set_colors(Eigen::RowVector3d(1, 1, 0.5));// Eigen::RowVector3d(1, 0, 1));
	viewer2.data().add_points(features, Eigen::RowVector3d(1, 0, 0));
	viewer2.data().add_points(results.front(), Eigen::RowVector3d(0, 1, 1));
	viewer2.data().add_points(spokepoints, Eigen::RowVector3d(0, 1, 0));

	viewer2.data().point_size = 10;
	viewer2.launch();
	return featuregroups;
}

std::vector<Mesh> ToothSegmentation::extractToothMeshes(const Mesh & mesh, const Eigen::VectorXd & harmonic_field, const std::vector<ToothSegmentation::ToothFeature>& teeth, const ToothSegmentation::ToothMeshExtractionParams& tme_params, bool visualize_steps)
{
	// idea: for every tooth, start a dfs (more efficient with std::vector) or bfs at every tooth feature point
	// and mark vertices belonging to the tooth and store their indices. Then assemble the mesh for the tooth.
	std::vector<Mesh> tooth_meshes;
	// 1 if vertex belongs to a tooth, 0 otherwise
	Eigen::VectorXi tooth_map(mesh.vertices().rows());
	tooth_map.setZero();

	// list of tooth vertex indices
	std::vector<Eigen::Index> tooth_indices;

	std::vector<Eigen::Index> stack;

	for (std::size_t t = 0; t < teeth.size(); ++t)
	{
		tooth_map.setZero();
		tooth_indices.clear();
		stack.clear();
		for (std::size_t f = 0; f < teeth[t].numFeaturePoints; ++f)
		{
			// do dfs starting from this feature
			stack.clear();
			stack.push_back(teeth[t].featurePointIndices[f]);
			while (!stack.empty())
			{
				// pop a node
				auto cidx = stack.back();
				stack.pop_back();
				if (tooth_map(cidx) == 0 && (t % 2 == 0 ? harmonic_field(cidx) < tme_params.even_tooth_threshold : harmonic_field(cidx) > tme_params.odd_tooth_threshold))
				{
					// mark node as tooth
					tooth_map(cidx) = 1;
					// add index to tooth index list
					tooth_indices.push_back(cidx);
					// push children onto stack
					for (const auto& a : mesh.adjacency_list()[cidx])
					{
						stack.push_back(a);
					}
				}
			}
		}

		if (visualize_steps)
		{
			igl::opengl::glfw::Viewer viewer;
			viewer.data().set_mesh(mesh.vertices(), mesh.faces());
			Eigen::MatrixXd C(mesh.vertices().rows(), 3);
			igl::jet(tooth_map, true, C);
			viewer.data().set_colors(C);
			viewer.launch();
		}

		// extract faces and vertices for this tooth

		// create tooth mesh
	}
	return tooth_meshes;
}