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

	static void segmentTeethFromMesh(const Mesh& mesh,
		const Eigen::Vector3d& mesh_up,
		const Eigen::Vector3d& mesh_right,
		std::vector<Mesh>& tooth_meshes,
		const CuspDetectionParams& cuspd_params = {0.4, 1000, 0.5},
		bool visualize_steps = false
		);

private:
	static void computeCusps(const Mesh& mesh, Eigen::MatrixXd& features, const CuspDetectionParams& cuspd_params, bool visualize_steps = false);
};


#endif