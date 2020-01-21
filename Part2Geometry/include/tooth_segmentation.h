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
		// smoothing step size
		double smoothing_step_size;
		// number of smoothing steps
		double smoothing_steps;
		// outlier detection
		double max_zscore;
		// fraction of local maxima to mean shift
		double ms_frac;
		// mean shift bandwidth param
		double ms_window_size; // [0, 1] fraction of bounding box diagonal
		// min change
		double ms_min_total_shift;
		// max iterations
		std::size_t ms_max_iterations;
		// feature collapse distance
		double ms_ft_collapse_dist;
	};

	struct ToothFeature
	{
		Eigen::DenseIndex featurePointIndices[4];
		Eigen::DenseIndex numFeaturePoints;
	};

	static void segmentTeethFromMesh(const Mesh& mesh,
		const Eigen::Vector3d& mesh_up,
		const Eigen::Vector3d& mesh_right,
		std::vector<Mesh>& tooth_meshes,
		const CuspDetectionParams& cuspd_params = {0.4, 1000, 0.5},
		double harmonic_field_w = 1000.0,
		bool visualize_steps = false
		);

private:
	static void computeCusps(const Mesh& mesh,
		Eigen::VectorXi& features,
		Eigen::VectorXd& mean_curvature,
		const CuspDetectionParams& cuspd_params,
		bool visualize_steps = false);
	static void calculateHarmonicField(const Mesh& mesh,
		const Eigen::VectorXd& mean_curvature,
		const std::vector<ToothFeature>& toothFeatures,
		const Eigen::VectorXi& cut_indices,
		Eigen::VectorXd& harmonic_field,
		double w = 1000.0,
		bool visualize_steps = false);
	static void cutMesh(Mesh& mesh,
		Eigen::VectorXi& cut_indices,
		Eigen::VectorXi& inv_index_map,
		Eigen::VectorXi& index_map,
		const Eigen::Vector3d& normal,
		const Eigen::Vector3d& plane_point);
	static double calcCotanWeight(const Eigen::Index& i, const Eigen::Index& j, const Mesh& mesh);

};


#endif