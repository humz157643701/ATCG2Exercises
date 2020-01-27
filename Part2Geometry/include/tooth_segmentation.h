#ifndef _TOOTH_SEGMENTATION_H_
#define _TOOTH_SEGMENTATION_H_
#include <mesh.h>
#include <Eigen/Dense>
#include <memory>

class ToothSegmentation
{
public:
	// ALGORITHM OVERVIEW
	// 1. find cusps features
	// 2. fit plane through found cusps
	// 3. slice mesh below teeth
	// 4. identify per-tooth features and place boundary conditions
	// 5. solve for harmonic field
	// 6. cut out individual tooth meshes

	struct CuspDetectionParams
	{
		// weighting parameter: curvature <-- [0, 1] --> mesh z value
		double alpha;
		double curve_exp;
		double height_exp;
		// fraction of bounding box height e [0,1] from jaw
		double min_feature_height;
		// fraction of local maxima to mean shift
		double os_frac;
		// optimum shift bandwidth param
		double os_window_size; // [0, 1] fraction of bounding box diagonal
		// min change
		double os_min_total_shift;
		// max iterations
		std::size_t os_max_iterations;
		// feature collapse distance
		double ft_collapse_dist;
		// small feature removal window size
		double small_ft_window_size;
		// small feature removal threshold
		double small_ft_threshold;
	};

	struct ToothFeature
	{
		std::vector<Eigen::DenseIndex> featurePointIndices;
		Eigen::DenseIndex numFeaturePoints;
	};

	struct HarmonicFieldParams
	{
		double w;
		double gamma_high;
		double gamma_low;
		double concavity_threshold;
	};

	struct MeanCurvatureParams
	{
		// smoothing step size
		double smoothing_step_size;
		// number of smoothing steps
		double smoothing_steps;
		// outlier detection
		double max_zscore;
	};

	struct ToothMeshExtractionParams
	{
		// even tooth threshold
		double even_tooth_threshold;
		// odd tooth threshold
		double odd_tooth_threshold;
	};

	static void segmentTeethFromMesh(const Mesh& mesh,
		const Eigen::Vector3d& approximate_mesh_up,
		const Eigen::Vector3d& mesh_right,
		std::vector<Mesh>& tooth_meshes,
		const CuspDetectionParams& cuspd_params = {0.4, 1000, 0.5},
		const HarmonicFieldParams& hf_params = { 1.0, 1.0 },
		const MeanCurvatureParams& mc_params = {0.00025, 50, 2.0},
		const ToothMeshExtractionParams& tme_params = {0.3, 0.6},
		bool visualize_steps = false);

private:
	static void computeMeanCurvature(const Mesh& mesh,
		Eigen::VectorXd& mean_curvature,
		const MeanCurvatureParams& mc_params,
		bool visualize_steps = false);
	static void computeCusps(const Mesh& mesh,
		const Eigen::Vector3d& mesh_up,
		const Eigen::VectorXd& mean_curvature,
		Eigen::VectorXi& features,
		const CuspDetectionParams& cuspd_params,
		bool visualize_steps = false);
	static void calculateHarmonicField(const Mesh& mesh,
		const Eigen::VectorXd& mean_curvature,
		const std::vector<ToothFeature>& toothFeatures,
		const Eigen::VectorXi& cut_indices,
		Eigen::VectorXd& harmonic_field,
		const HarmonicFieldParams& hf_params,
		bool visualize_steps = false);
	static void cutMesh(Mesh& mesh,
		Eigen::VectorXi& cut_indices,
		Eigen::VectorXi& inv_index_map,
		Eigen::VectorXi& index_map,
		const Eigen::Vector3d& normal,
		const Eigen::Vector3d& plane_point);
	static double calcCotanWeight(const Eigen::Index& i,
		const Eigen::Index& j,
		const Mesh& mesh);
	static double calcCurvatureWeight(const Eigen::Index& i,
		const Eigen::Index& j,
		const Eigen::VectorXd& mean_curvature,
		const HarmonicFieldParams& hf_params);
	static Eigen::Vector3d estimateUpVector(const Mesh& mesh, const Eigen::Vector3d& approximate_up);
	static Eigen::Vector3d estimateUpVector(const Eigen::MatrixXd& points, const Eigen::Vector3d& approximate_up);

	static std::pair<Eigen::Vector3d, Eigen::Vector3d> fitPlane(const Eigen::VectorXi& featureindices, const Mesh& mesh, Eigen::VectorXi idmap);

	static std::vector<std::vector<size_t>> segmentFeatures(const Eigen::VectorXi& featureindices, Eigen::VectorXd meancurvature, const Mesh& mesh, Eigen::VectorXi idmap);

	static std::vector<Mesh> extractToothMeshes(const Mesh& mesh, const Eigen::VectorXd& harmonic_field, const std::vector<ToothFeature>& teeth, const ToothMeshExtractionParams& tme_params, bool visualize_steps = false);
};


#endif