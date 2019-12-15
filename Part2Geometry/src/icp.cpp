#include <icp.h>
#include <Eigen/Core>

ICPAligner::ICPAligner(const Eigen::MatrixXd & target_points) :
	m_target_octree(nullptr)
{
	std::vector<Vec3> tpoints;
	tpoints.reserve(static_cast<std::size_t>(target_points.rows()));
	for (Eigen::DenseIndex i = 0; i < target_points.rows(); ++i)
		tpoints.push_back(Vec3{ target_points(i, 0), target_points(i, 1), target_points(i, 2) });
	buildOctree(tpoints);
}

double ICPAligner::align(Eigen::Matrix3d & optimal_rotation, Eigen::Vector3d & optimal_translation, Eigen::MatrixXd & query_points, const Eigen::MatrixXd & target_points, const Eigen::MatrixXd & target_normals, const Eigen::MatrixXd & query_normals, double max_distance, double min_normal_cos_theta, double min_err)
{
	optimal_rotation.setIdentity();
	optimal_translation.setZero();

	
}

void ICPAligner::setTargetPoints(const Eigen::MatrixXd & target_points)
{
	std::vector<Vec3> tpoints;
	tpoints.reserve(static_cast<std::size_t>(target_points.rows()));
	for (Eigen::DenseIndex i = 0; i < target_points.rows(); ++i)
		tpoints.push_back(Vec3{ target_points(i, 0), target_points(i, 1), target_points(i, 2) });
	buildOctree(tpoints);
}

void ICPAligner::buildOctree(const std::vector<Vec3>& points)
{
	m_target_octree.reset(new Octree(points));
	m_target_octree->build(5);
}

double ICPAligner::calcSquaredError(const Eigen::MatrixXd & a, const Eigen::MatrixXd & b)
{
	return (a.array() * b.array()).rowwise().sum().sum();
}

double ICPAligner::calcSquaredError(const Eigen::MatrixXd & a, const Eigen::MatrixXd & b, const Eigen::MatrixXi & correspondences)
{
	auto indicesa = correspondences.col(0).array();
	auto indicesb = correspondences.col(1).array();

	return (a(indicesa, Eigen::all) * b.array()).rowwise().sum().sum();
}

void ICPAligner::genCorrespondences(Eigen::MatrixXi& correspondences,
	Eigen::MatrixXd& distances,
	const Eigen::MatrixXd & query_points,
	const Eigen::MatrixXd & target_points, 
	const Eigen::MatrixXd & target_normals, 
	const Eigen::MatrixXd & query_normals, 
	double max_distance, 
	double min_normal_cos_theta)
{
	static std::vector<Eigen::RowVector2i> corr_cache;
	static std::vector<Eigen::RowVector2d> dist_cache;
	corr_cache.clear();
	dist_cache.clear();

	for (Eigen::DenseIndex q = 0; q < query_points.rows(); ++q)
	{
		auto res = m_target_octree->query_knn(Vec3{ query_points(q, 0), query_points(q, 1), query_points(q, 2) }, 1);
		auto nnidx = static_cast<Eigen::DenseIndex>(res.idx_dist_pair[0].first);
		double distance = static_cast<double>(res.idx_dist_pair[0].second);
		double costheta = query_normals.row(q).dot(target_normals.row(nnidx));
		if (distance <= max_distance && costheta >= min_normal_cos_theta)
		{
			corr_cache.push_back(Eigen::RowVector2i{ q, nnidx });
			dist_cache.push_back(Eigen::Vector2d{ distance, costheta });
		}
	}

	for (std::size_t i = 0; i < corr_cache.size(); ++i)
	{
		correspondences.row(static_cast<Eigen::DenseIndex>(i)) = corr_cache[i];
		distances.row(static_cast<Eigen::DenseIndex>(i)) = dist_cache[i];
	}
}
