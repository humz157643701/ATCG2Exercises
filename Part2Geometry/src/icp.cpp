#include <icp.h>
#include <Eigen/Core>
#include <iostream>
#include <igl/opengl/glfw/Viewer.h>

ICPAligner::ICPAligner(const Eigen::MatrixXd & target_points) :
	m_target_octree(nullptr),
	m_octree_data()
{
	m_octree_data.reserve(static_cast<std::size_t>(target_points.rows()));
	for (Eigen::DenseIndex i = 0; i < target_points.rows(); ++i)
		m_octree_data.push_back(Vec3{ static_cast<float>(target_points(i, 0)), static_cast<float>(target_points(i, 1)), static_cast<float>(target_points(i, 2) )});
	buildOctree(m_octree_data);
}

double ICPAligner::align(Eigen::Matrix3d & optimal_rotation, Eigen::Vector3d & optimal_translation, Eigen::MatrixXd & query_points, const Eigen::MatrixXd & target_points, const Eigen::MatrixXd & target_normals, const Eigen::MatrixXd & query_normals, double max_distance, double min_normal_cos_theta, double min_err)
{
	optimal_rotation.setIdentity();
	optimal_translation.setZero();

	Eigen::MatrixXi correspondences;
	Eigen::MatrixXd distances;

	double error = std::numeric_limits<double>::max();
	std::size_t itct = 0;
	genCorrespondences(correspondences, distances, query_points, target_points, target_normals, query_normals, max_distance, min_normal_cos_theta);
	while (true)
	{		
		igl::opengl::glfw::Viewer viewr;
		viewr.data().clear();
		viewr.data().point_size = 1.0;
		viewr.data().add_points(target_points(correspondences.col(1).array(), Eigen::all) * 0.05, Eigen::RowVector3d(1.0, 1.0, 1.0));
		viewr.data().add_points(query_points * 0.05, Eigen::RowVector3d(1.0, 0.0, 0.0));
		viewr.launch();
		genCorrespondences(correspondences, distances, query_points, target_points, target_normals, query_normals, max_distance, min_normal_cos_theta);
		double newerror = calcMeanSquaredError(query_points, target_points, correspondences);
		if (std::abs(newerror - error) < min_err)
			return newerror;
		error = newerror;
		std::cout << "ICP Iteration: " << itct << ", Current Error: " << error << "\n";		
		optimalRigidTransform(query_points, target_points, correspondences, distances, optimal_rotation, optimal_translation);
		applyRigidTransform(query_points, optimal_rotation, optimal_translation);	
		itct++;
		
	}
}

void ICPAligner::setTargetPoints(const Eigen::MatrixXd & target_points)
{
	m_octree_data.clear();
	for (Eigen::DenseIndex i = 0; i < target_points.rows(); ++i)
		m_octree_data.push_back(Vec3{ static_cast<float>(target_points(i, 0)), static_cast<float>(target_points(i, 1)), static_cast<float>(target_points(i, 2) )});
	buildOctree(m_octree_data);
}

void ICPAligner::applyRigidTransform(Eigen::MatrixXd & points, const Eigen::Matrix3d & optimal_rotation, const Eigen::Vector3d & optimal_translation)
{
	points = points * optimal_rotation.transpose();
	points.rowwise() += optimal_translation.transpose();
}

void ICPAligner::buildOctree(const std::vector<Vec3>& points)
{
	m_target_octree.reset(new Octree(points));
	m_target_octree->build(1);
}

double ICPAligner::calcMeanSquaredError(const Eigen::MatrixXd & a, const Eigen::MatrixXd & b)
{
	return (a.array() * b.array()).rowwise().sum().sum() / static_cast<double>(a.rows());
}

double ICPAligner::calcMeanSquaredError(const Eigen::MatrixXd & a, const Eigen::MatrixXd & b, const Eigen::MatrixXi & correspondences)
{
	return (a(correspondences.col(0).array(), Eigen::all).array() * b(correspondences.col(1).array(), Eigen::all).array()).rowwise().sum().sum() / static_cast<double>(correspondences.rows());
}

void ICPAligner::calcWeights(Eigen::VectorXd & weights, const Eigen::MatrixXd & distances)
{
	weights.resize(distances.rows(), 1);
	weights = distances.exp
}

void ICPAligner::optimalRigidTransform(const Eigen::MatrixXd & query_points, const Eigen::MatrixXd & target_points, const Eigen::MatrixXi & correspondences, const Eigen::VectorXd & weights, Eigen::Matrix3d& optimal_rotation, Eigen::Vector3d& optimal_translation)
{
	if (correspondences.rows() >= 2)
	{
		auto Q = query_points(correspondences.col(0).array(), Eigen::all);
		auto T = target_points(correspondences.col(1).array(), Eigen::all);
		Eigen::RowVector3d Q_mean = Q.colwise().mean();
		Eigen::RowVector3d T_mean = T.colwise().mean();
		Q.rowwise() -= Q_mean;
		T.rowwise() -= T_mean;
		Eigen::Matrix3d H;
		H.setZero();
		for (Eigen::DenseIndex i = 0; i < Q.rows(); ++i)
		{
			H += Q(i, Eigen::all) * T(i, Eigen::all).transpose() * weights(i);
		}
		Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Eigen::Matrix3d F = Eigen::MatrixXd::Identity(3, 3);
		F(2, 2) = (svd.matrixV() * svd.matrixU().transpose()).determinant();
		Eigen::Matrix3d R = svd.matrixV() * F * svd.matrixU().transpose();
		Eigen::RowVector3d t = T_mean - (Q_mean * R.transpose());
		optimal_rotation = R;
		optimal_translation = t.transpose();
	}
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
		auto res = m_target_octree->query_knn(Vec3{ static_cast<float>(query_points(q, 0)), static_cast<float>(query_points(q, 1)), static_cast<float>(query_points(q, 2)) }, 1);
		auto nnidx = static_cast<Eigen::DenseIndex>(res.idx_dist_pair[0].first);
		double distance = static_cast<double>(res.idx_dist_pair[0].second);
		double costheta = query_normals.row(q).dot(target_normals.row(nnidx));
		if (distance <= max_distance && costheta >= min_normal_cos_theta)
		{
			corr_cache.push_back(Eigen::RowVector2i{ q, nnidx });
			dist_cache.push_back(Eigen::Vector2d{ distance, costheta });
		}
	}

	correspondences.resize(corr_cache.size(), 2);
	distances.resize(dist_cache.size(), 2);

	for (std::size_t i = 0; i < corr_cache.size(); ++i)
	{
		correspondences.row(static_cast<Eigen::DenseIndex>(i)) = corr_cache[i];
		distances.row(static_cast<Eigen::DenseIndex>(i)) = dist_cache[i];
	}
}
