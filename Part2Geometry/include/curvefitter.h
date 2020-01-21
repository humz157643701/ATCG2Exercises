#ifndef _CURVEFITTER_H_
#define _CURVEFITTER_H_
#include <Eigen/Dense>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <mesh.h>
class CurveFitter
{
public:
	static Eigen::MatrixXd fitCurve(const Eigen::MatrixXd& features)
	{

		Eigen::MatrixXd M(5, 5);
		M.fill(0);
		Eigen::MatrixXd polied(features.rows(), 5);
		//Construct beta
		Eigen::VectorXd b_3(5);
		auto wf = [](double x, double xm) {return 1 / (std::abs(x - xm) + 1); };
		double x_middle = 0;
		double z_max = -999999;

		for (auto& row : features.rowwise())
		{
			if (row(2) > z_max)
			{
				z_max = row(2);
				x_middle = row(0);
			}

		}
		for (size_t x = 0; x < polied.rows(); ++x)
		{
			polied.row(x) << std::pow(features(x, 0), 4), std::pow(features(x, 0), 3), std::pow(features(x, 0), 2), std::pow(features(x, 0), 1), 1;
			M += Eigen::MatrixXd(polied.row(x).transpose() * polied.row(x));
			//std::cout << "polied row: " << polied.row(x).transpose() * polied.row(x) << std::endl;
			//polied.row(x) *= wf(features(x, 0), x_middle);
		}

		//M = polied.colwise().sum().transpose() * polied.colwise().sum();
		std::cout << "polied: " <<  polied << std::endl;
		std::cout << "M: \n" << M << std::endl;
		for (size_t x = 0; x < 5; ++x)
		{
			double v = 0;
			for (auto& xyz : features.rowwise())
			{
				v +=  std::pow(xyz(0), 4 - x) * xyz(2); //wf(xyz(0), x_middle) *
			}

			b_3.row(x) << v;
		}
		Eigen::MatrixXd svd = M.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b_3);
		std::cout << "svd sol:\n" << svd << std::endl;
		std::cout << "svd reprod:\n" << M*svd << std::endl;

		return svd;
		//return { svd.matrixV().bottomRows(1).row(0) , center };
		//ÜLTRA
	}
	static std::vector<Eigen::MatrixXd> segmentFeatures(const Eigen::MatrixXd& features, const Eigen::MatrixXd curveparams, double curvey, const Mesh& mesh)
	{
		double min_x, max_x;
 		min_x = features.colwise().minCoeff()(0);
		max_x = features.colwise().maxCoeff()(0);
		std::cout << "min, max x: " << min_x << " | " << max_x << std::endl;
		std::vector<std::vector<std::pair<Eigen::Vector3d,Eigen::Vector3d>>> spokes{};
		std::cout << "curve: " << curveparams << std::endl;
		int num_spokes = 100;
		spokes.reserve(100);
		std::vector<Eigen::Vector3d> curvepoints;
		double fpymax = features.colwise().maxCoeff()(1);
		std::cout << "fpymax: " << fpymax << std::endl;
		for (double x = min_x; x <= max_x; x += (max_x - min_x) / 100.0)
		{
			
			double z = curveparams(0, 0) * std::pow(x, 4) + curveparams(1, 0) * std::pow(x, 3) + curveparams(2, 0) * std::pow(x, 2) + curveparams(3, 0) *std::pow(x, 1) + curveparams(4, 0);
			double dz = curveparams(0, 0) * std::pow(x, 3) * 4 + curveparams(1, 0) * std::pow(x, 2) * 3 + curveparams(2, 0) *std::pow(x, 1) * 2 + curveparams(3, 0);
			curvepoints.push_back({ x,fpymax,z });
			//z = dz*z+b
			Eigen::Vector2d normal{ 1, -(1/dz) };
			normal.normalize();
			spokes.push_back({});
			for (int theta = -3; theta <= 3; ++theta)
			{
				double t = theta * 10 * (M_PI / 180.0);
				Eigen::Matrix2d rot;
				rot << std::cos(t), -std::sin(t), std::sin(t), std::cos(t);
				//std::cout << "norm: " << normal << std::endl;
				auto rotnorm = rot* normal;
				//std::cout << "rotated: " << rotnorm << std::endl;
				spokes.back().push_back({ { rotnorm(0), 0, rotnorm(1) }, { x, fpymax, z } });
			}
		}

		// For every spoke, create samples on spoke and find find closest vertice to plane spun up by thhe normal vector and up vector OR vertices in radius of point, store y value.
		// Pick highest y value along spoke -> among all spokes pick the one with the lowest max
		std::vector<std::vector<double>> depths;
		depths.reserve(100);
		std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> finalspokes;
		double verticesminy = mesh.vertices().colwise().minCoeff()(1);
		std::cout << "verticesminy: " << fpymax << std::endl;
		double stepsize = std::abs(verticesminy - curvey) / 500.0;
		std::cout << "stepsize: " << fpymax << std::endl;
		for (size_t outer = 50; outer < spokes.size(); outer+=1)
		{
			double minspokey = std::numeric_limits<double>::max();
			size_t minspokeIdx = 0;
			for (size_t inner = 0; inner < spokes[outer].size(); inner++)
			{
				double maxy = std::numeric_limits<double>::lowest();
				for (double x = -1; x <= 1; x += 0.1)
				{
					Eigen::Vector3d pos = spokes[outer][inner].second + (spokes[outer][inner].first * x * 5);
					//Let's go down step by step in negative y direction until we find a close vertice
					std::vector<std::pair<Eigen::Index, double>> srchres;
					for (double y = pos(1); y >= verticesminy; y -= stepsize)
					{
						auto neighbor = mesh.kdtree().index->radiusSearch(pos.data(), stepsize*stepsize, srchres, {});
						if (srchres.size() > 0)
						{
							std::cout << " srchres size: " << srchres.size() << "y: " << y<<  std::endl;
							auto r = mesh.vertices()(srchres.front().first, 1);
							if (r > maxy)
							{
								maxy = r;
							}
							break;
						}
					}
				}
				if (maxy < minspokey)
				{
					minspokey = maxy;
					minspokeIdx = inner;
				}
			}
			finalspokes.push_back({ spokes[outer][minspokeIdx] });
		}
		Eigen::MatrixXd spokepoints(20*finalspokes.size()+1, 3);
		for (size_t outer = 0; outer < finalspokes.size(); outer += 5)
		{
			int x1 = 0;
			for (double x = -1; x <= 1; x += 0.1)
			{
				Eigen::Vector3d pos = finalspokes[outer].second + (finalspokes[outer].first * x*5);
				spokepoints.row(outer * 20 + x1) = pos;
				x1++;
			}
		}
		// Remove spokes which have no feature points in between one another


		std::vector<Eigen::MatrixXd> results;
		results.push_back(Eigen::MatrixXd{ curvepoints.size(), 3 });
		int m = spokes.size() * spokes[0].size();
		for (size_t x = 0; x < curvepoints.size(); x++)
		{
			results.back().row(x) = curvepoints[x];
		}
		results.push_back(spokepoints);
		return results;
	}
};

#endif