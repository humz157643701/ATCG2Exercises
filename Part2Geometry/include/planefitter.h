#ifndef _PLANEFITTER_H_
#define _PLANEFITTER_H_
#include <Eigen/Dense>
#include <iostream>
class PlaneFitter
{
public:
	static std::pair<Eigen::Vector3d, Eigen::Vector3d> fitPlane(const Eigen::MatrixXd& features)
	{
		Eigen::Vector3d center = features.colwise().mean();
		auto normed = Eigen::MatrixXd(features);
		normed.rowwise()-=center.transpose();
		auto svd = normed.bdcSvd(Eigen::ComputeFullU|Eigen::ComputeFullV);
		std::cout << "V: \n" << svd.matrixV().transpose() << std::endl;
		return { svd.matrixV().transpose().bottomRows(1).row(0) , center };
	}
};

#endif