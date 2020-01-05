#ifndef _ICP_H_
#define _ICP_H_
#include <nanoflann.hpp>
#include <memory>
#include <Eigen/Dense>
#include <limits>

struct ICPParams
{
	double max_distance = std::numeric_limits<double>::max();
	double min_normal_cos_theta = -1.0;
	double min_err = 1e-5;
	double min_err_change = 1e-8;
	std::size_t max_iterations = 200;
};

class ICPAligner
{
	using kdtree_t = nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixXd, 3, nanoflann::metric_L2>;
public:
	ICPAligner(const Eigen::MatrixXd& target_points);

	double align(Eigen::Matrix3d& optimal_rotation,
		Eigen::Vector3d& optimal_translation,
		Eigen::MatrixXd& query_points,
		const Eigen::MatrixXd & target_points,
		const Eigen::MatrixXd& target_normals,
		const Eigen::MatrixXd& query_normals,
		double max_distance = std::numeric_limits<double>::max(),
		double min_normal_cos_theta = -1.0,
		double min_err = 1e-5,
		double min_err_change = 1e-8,
		size_t max_iterations = 200);
	double align(Eigen::Matrix3d& optimal_rotation,
		Eigen::Vector3d& optimal_translation,
		Eigen::MatrixXd& query_points,
		const Eigen::MatrixXd & target_points,
		const Eigen::MatrixXd& target_normals,
		const Eigen::MatrixXd& query_normals,
		const ICPParams& params);		
	void setTargetPoints(const Eigen::MatrixXd& target_points);
	static void applyRigidTransform(Eigen::MatrixXd& points, const Eigen::Matrix3d& optimal_rotation, const Eigen::Vector3d& optimal_translation, bool reorthonormalize_rotation = false);
private:
	void buildKDTree(const Eigen::MatrixXd& points);
	double calcMAE(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b);
	double calcMAE(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b, const Eigen::MatrixXi& correspondences);
	void calcWeights(Eigen::VectorXd& weights, const Eigen::MatrixXd& distances);
	void optimalRigidTransform(const Eigen::MatrixXd& query_points, const Eigen::MatrixXd& target_points, const Eigen::MatrixXi& correspondences, const Eigen::VectorXd& weights, Eigen::Matrix3d& optimal_rotation, Eigen::Vector3d& optimal_translation);
	void genCorrespondences(
		Eigen::MatrixXi& correspondences,
		Eigen::MatrixXd& distances,
		const Eigen::MatrixXd& query_points,
		const Eigen::MatrixXd& target_points,
		const Eigen::MatrixXd& target_normals,
		const Eigen::MatrixXd& query_normals,
		double max_distance = std::numeric_limits<double>::max(),
		double min_normal_cos_theta = -1.0);
	std::unique_ptr<kdtree_t> m_target_kdtree;
};

#endif