#ifndef _ICP_H_
#define _ICP_H_
#include <Octree.h>
#include <memory>
#include <Eigen/Dense>
#include <limits>

class ICPAligner
{
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
		double min_err = 1e-5);
	void setTargetPoints(const Eigen::MatrixXd& target_points);
private:
	void buildOctree(const std::vector<Vec3>& points);
	double calcSquaredError(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b);
	double calcSquaredError(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b, const Eigen::MatrixXi& correspondences);
	void genCorrespondences(
		Eigen::MatrixXi& correspondences,
		Eigen::MatrixXd& distances,
		const Eigen::MatrixXd& query_points,
		const Eigen::MatrixXd& target_points,
		const Eigen::MatrixXd& target_normals,
		const Eigen::MatrixXd& query_normals,
		double max_distance = std::numeric_limits<double>::max(),
		double min_normal_cos_theta = -1.0);
	std::unique_ptr<Octree> m_target_octree;
};

#endif