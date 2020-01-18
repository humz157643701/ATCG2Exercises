#include "..\include\tooth_segmentation.h"
#include <iostream>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/jet.h>

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
	Eigen::VectorXi cuspIndices;
	computeCusps(working_mesh, cuspIndices, cuspd_params, true);

	// fit plane

	// cut teeth

	// fit curve

	// assign features to teeth

	// harmonic field stuff
}

void ToothSegmentation::computeCusps(const Mesh& mesh, Eigen::VectorXi& featureIndices, const ToothSegmentation::CuspDetectionParams& cuspd_params, bool visualize_steps)
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
	Eigen::VectorXd weights = ((1.0 - cuspd_params.alpha) * curv_weights.array().pow(cuspd_params.curve_exp) + cuspd_params.alpha * height_weights.array().pow(cuspd_params.height_exp)).matrix();

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


}
