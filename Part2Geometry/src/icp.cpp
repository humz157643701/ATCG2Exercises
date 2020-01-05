#include <icp.h>
#include <Eigen/Core>
#include <iostream>
#include <igl/opengl/glfw/Viewer.h>
#include <exception>
#include <stdexcept>

ICPAligner::ICPAligner(const Eigen::MatrixXd & target_points) :
	m_target_kdtree(nullptr)
{
	buildKDTree(target_points);
}

double ICPAligner::align(Eigen::Matrix3d & optimal_rotation,
	Eigen::Vector3d & optimal_translation,
	Eigen::MatrixXd & query_points,
	const Eigen::MatrixXd & target_points,
	const Eigen::MatrixXd & target_normals,
	const Eigen::MatrixXd & query_normals,
	double max_distance,
	double min_normal_cos_theta,
	double min_err,
	double min_err_change,
	size_t max_iterations)
{
	optimal_rotation.setIdentity();
	optimal_translation.setZero();
	Eigen::Matrix3d optimal_rotation_delta;
	Eigen::Vector3d optimal_translation_delta;
	optimal_rotation_delta.setIdentity();
	optimal_translation_delta.setZero();

	Eigen::MatrixXi correspondences;
	Eigen::MatrixXd distances;
	Eigen::VectorXd weights;

	double error = std::numeric_limits<double>::max();
	std::size_t itct = 0;
	while (true)
	{
		std::cout << "\nICP: Iteration " << itct << "\n";
		std::cout << "ICP: Calculating correspondences...\n";
		genCorrespondences(correspondences, distances, query_points, target_points, target_normals, query_normals, max_distance, min_normal_cos_theta);
		std::cout << "ICP: " << correspondences.rows() <<  " correspondences found.\n";
		if(correspondences.rows() == 0)
			throw std::logic_error("ICP: No correspondences found.\n");
		std::cout << "ICP: Calculating error...\n";
		double newerror = calcMAE(query_points, target_points, correspondences);
		std::cout << "ICP: Current Error: " << newerror << "\n";
		if (std::abs(newerror - error) < min_err_change || newerror < min_err)
		{
			std::cout << "ICP: Registration finished!\n";
			return newerror;
		}
		else if (itct > max_iterations)
		{
			std::cout << "ICP: Registration failed.\n";
			return newerror;
		}
		error = newerror;
		
		std::cout << "ICP: Aligning correspondences...\n";
		optimalRigidTransform(query_points, target_points, correspondences, weights, optimal_rotation_delta, optimal_translation_delta);
		std::cout << "ICP: Transforming query set...\n";
		applyRigidTransform(query_points, optimal_rotation_delta, optimal_translation_delta);
		// update global transformation
		optimal_rotation = optimal_rotation_delta * optimal_rotation;
		optimal_translation = optimal_rotation_delta * optimal_translation + optimal_translation_delta;
		itct++;		
	}
}

double ICPAligner::align(Eigen::Matrix3d & optimal_rotation, Eigen::Vector3d & optimal_translation, Eigen::MatrixXd & query_points, const Eigen::MatrixXd & target_points, const Eigen::MatrixXd & target_normals, const Eigen::MatrixXd & query_normals, const ICPParams & params)
{
	return align(optimal_rotation, optimal_translation, query_points, target_points, target_normals, query_normals, params.max_distance, params.min_normal_cos_theta, params.min_err, params.min_err_change, params.max_iterations);
}

void ICPAligner::setTargetPoints(const Eigen::MatrixXd & target_points)
{
	buildKDTree(target_points);
}

void ICPAligner::applyRigidTransform(Eigen::MatrixXd & points, const Eigen::Matrix3d & optimal_rotation, const Eigen::Vector3d & optimal_translation, bool reorthonormalize_rotation)
{
	//// re-orthonormalize rotation matrix
	//Eigen::Matrix3d R(optimal_rotation);
	//if (reorthonormalize_rotation)
	//{
	//	Eigen::JacobiSVD<Eigen::Matrix3d> svd(optimal_rotation, Eigen::ComputeFullU | Eigen::ComputeFullV);
	//	Eigen::Matrix3d F = Eigen::MatrixXd::Identity(3, 3);
	//	F(2, 2) = (svd.matrixU() * svd.matrixV().transpose()).determinant();
	//	R = svd.matrixU() * F * svd.matrixV().transpose();
	//}

	points *= optimal_rotation.transpose();
	points.rowwise() += optimal_translation.transpose();
}

void ICPAligner::buildKDTree(const Eigen::MatrixXd& points)
{
	m_target_kdtree.reset(new kdtree_t(3, points));
	m_target_kdtree->index->buildIndex();
}

double ICPAligner::calcMAE(const Eigen::MatrixXd & a, const Eigen::MatrixXd & b)
{
	return (a - b).rowwise().norm().mean();
}

double ICPAligner::calcMAE(const Eigen::MatrixXd & a, const Eigen::MatrixXd & b, const Eigen::MatrixXi & correspondences)
{
	return (a(correspondences.col(0).array(), Eigen::all) - b(correspondences.col(1).array(), Eigen::all)).rowwise().norm().mean();
}

void ICPAligner::calcWeights(Eigen::VectorXd & weights, const Eigen::MatrixXd & distances)
{
	weights.resize(distances.rows(), 1);
	weights = (-(distances.col(0).array() * distances.col(0).array())).exp() * (distances.col(1).array()).max(0.0);
}

void ICPAligner::optimalRigidTransform(const Eigen::MatrixXd & query_points, const Eigen::MatrixXd & target_points, const Eigen::MatrixXi & correspondences, const Eigen::VectorXd & weights, Eigen::Matrix3d& optimal_rotation, Eigen::Vector3d& optimal_translation)
{
	if (correspondences.rows() >= 2)
	{
		auto Q = query_points(correspondences.col(0).array(), Eigen::all);
		auto T = target_points(correspondences.col(1).array(), Eigen::all);
		Eigen::RowVector3d Q_mean = Q.colwise().mean();
		Eigen::RowVector3d T_mean = T.colwise().mean();
		Eigen::Matrix3d H = (Q.rowwise() - Q_mean).transpose() * (T.rowwise() - T_mean);
		Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Eigen::Matrix3d F = Eigen::MatrixXd::Identity(3, 3);
		F(2, 2) = (svd.matrixV() * svd.matrixU().transpose()).determinant();
		Eigen::Matrix3d R = svd.matrixV() /* F*/ * svd.matrixU().transpose();
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
	static std::vector<Eigen::RowVector2i> corr_cache(query_points.rows());
	//static std::vector<Eigen::RowVector2d> dist_cache(query_points.rows());
	corr_cache.clear();
	//dist_cache.clear();

	for (Eigen::DenseIndex q = 0; q < query_points.rows(); ++q)
	{
		double qp[] = { query_points(q, 0), query_points(q, 1), query_points(q, 2) };
		Eigen::DenseIndex nnidx;
		double distance;
		m_target_kdtree->query(&qp[0], 1, &nnidx, &distance);
		distance = std::sqrt(distance);
		double costheta = query_normals.row(q).dot(target_normals.row(nnidx));
		if (distance <= max_distance && costheta >= min_normal_cos_theta)
		{
			corr_cache.push_back(Eigen::RowVector2i{ q, nnidx });
			//dist_cache.push_back(Eigen::Vector2d{ distance, costheta });
		}
	}

	correspondences.resize(corr_cache.size(), 2);
	//distances.resize(dist_cache.size(), 2);

	for (std::size_t i = 0; i < corr_cache.size(); ++i)
	{
		correspondences.row(static_cast<Eigen::DenseIndex>(i)) = corr_cache[i];
		//distances.row(static_cast<Eigen::DenseIndex>(i)) = dist_cache[i];
	}
}
