#include <iostream>
#include<igl/opengl/glfw/Viewer.h>
#include<igl/readOBJ.h>
#include <icp.h>
#include <Eigen/Geometry>



int main(int argc, char *argv[])
{
	try
	{

		//Eigen::MatrixXd planeV11(4, 3);
		//Eigen::MatrixXi planeF11(2, 3);
		//Eigen::MatrixXd planeC1(planeV11.rows(), 3);
		//planeF11 << 0, 1, 2, 3, 2, 1;
		//Eigen::Vector3d meme(1, 0, 0);
		//// define plane points
		//Eigen::Vector3d perp11 = Eigen::Vector3d(meme(2), meme(1), meme(0)).cross(meme);
		//std::cout << perp11 << "\n";
		//Eigen::Vector3d perp22 = perp11.cross(meme);
		//planeV11(0, Eigen::all) = perp11 * 100 + perp22 * 100;
		//planeV11(1, Eigen::all) = perp11 * -100 + perp22 * 100;
		//planeV11(2, Eigen::all) = perp11 * 100 + perp22 * -100;
		//planeV11(3, Eigen::all) = perp11 * -100 + perp22 * -100;
		////planeV1 << 10.0, 0.0, 0.0, -10.0, 0.0, 0.0, 10.0, -10.0, 0.0, -10.0, -10.0, 0.0;
		//std::cout << planeF11 << "\n";
		//std::cout << planeV11;
		//planeC1 << Eigen::RowVector3d(1.0, 1.0, 1.0).replicate(planeV11.rows(), 1);


		//// Plot the mesh
		//igl::opengl::glfw::Viewer viewer22;
		//viewer22.data().set_mesh(planeV11, planeF11);
		//viewer22.data().set_colors(planeC1);
		//viewer22.launch();


		Eigen::MatrixXd V1;
		Eigen::MatrixXi F1;
		std::cout << "Loading meshes...\n";
		igl::readOBJ("assets/models/RD-01/16021_OnyxCeph3_Export_OK-A.obj", V1, F1);
		Eigen::MatrixXd V2(V1);
		Eigen::MatrixXi F2(F1);
		//ICPAligner::applyRigidTransform(V2, Eigen::Matrix3d(Eigen::AngleAxisd((5.0 / 180.0) * EIGEN_PI, Eigen::Vector3d(1.0, 2.0, 0.0))), Eigen::Vector3d(2.0, 0.0, 0.0));

		Eigen::MatrixXd N1;
		Eigen::MatrixXd N2;

		igl::per_vertex_normals(V1, F1, N1);
		igl::per_vertex_normals(V2, F2, N2);

		Eigen::MatrixXd C1(V1.rows(), 3);
		Eigen::MatrixXd C2(V2.rows(), 3);
		C1 << Eigen::RowVector3d(1.0, 1.0, 1.0).replicate(V1.rows(), 1);
		C2 << Eigen::RowVector3d(1.0, 0.0, 0.0).replicate(V2.rows(), 1);

		Eigen::Matrix3d rot;
		Eigen::Vector3d trans;

		Eigen::MatrixXd V(V1.rows() + V2.rows(), V1.cols());
		Eigen::MatrixXi F(F1.rows() + F2.rows(), F1.cols());
		Eigen::MatrixXd C(C1.rows() + C2.rows(), C1.cols());

		V << V1, V2;
		F << F1, (F2.array() + V1.rows());
		C << C1, C2;

		

		// Plot the mesh
		igl::opengl::glfw::Viewer viewer;
		//viewer.data().set_mesh(V, F);
		//viewer.data().set_colors(C);
		//viewer.launch();

		Eigen::Vector3d centerofmass{ 0,0,0 };
		Eigen::Vector3d planenormal{ 1,0,0 };
		double d_originplanedistance = 0;
		Eigen::Matrix3d reflectionmatrix;
		reflectionmatrix << 1, 0, 0, 0, 1, 0, 0, 0, 1;
		std::cout << reflectionmatrix << "\n";
		reflectionmatrix = reflectionmatrix - 2 * (planenormal * planenormal.transpose());
		std::cout << reflectionmatrix << "\n";

		for (size_t x = 0; x < V1.rows(); ++x)
		{
			centerofmass(0) += V1(x, 0);
			centerofmass(1) += V1(x, 1);
			centerofmass(2) += V1(x, 2);
		}
		centerofmass /= V1.rows();
		d_originplanedistance = (Eigen::Vector3d(0, 0, 0) - centerofmass).dot(planenormal);
		std::cout << centerofmass << "\n";
		std::cout << d_originplanedistance;

		for (size_t x = 0; x < V2.rows(); ++x)
		{
			Eigen::Vector3d refl = reflectionmatrix * V2(x, Eigen::all).transpose() + 2 * d_originplanedistance * planenormal;

			
			V2(x, Eigen::all) = refl;
		}

		for (size_t x = 0; x < N2.rows(); ++x)
		{
			Eigen::Vector3d refl = reflectionmatrix * N2(x, Eigen::all).transpose() + 2 * d_originplanedistance * planenormal;


			N2(x, Eigen::all) = refl;
		}

		ICPAligner icp(V1);
		// align transforms the query set (V2) and returns optimal translation and rotation
		icp.align(rot, trans, V2, V1, N1, N2, 50.0, 0.2, 1e-2, 1e-4, 500);

		auto reflectedrot = reflectionmatrix * rot.transpose();
		Eigen::EigenSolver<Eigen::MatrixXd> es(reflectedrot);
		int smallesteigenidx = 0;
		double smallesteigendiff = 1;
		for (size_t x = 0; x < es.eigenvalues().rows(); ++x)
		{
			auto diff = es.eigenvalues()(x, 0).real() + 1;
			if (diff < smallesteigendiff)
			{
				smallesteigenidx = x;
				smallesteigendiff = diff;
			}
		}
		std::cout << es.eigenvalues()(smallesteigenidx).real() << " | " << smallesteigendiff << "\n";
		std::cout << es.eigenvectors().col(smallesteigenidx)(Eigen::all, 0).real() << "\n";

		// New symm plane point
		Eigen::Vector3d newplanepoint = 0.5 * (rot * (2 * d_originplanedistance * planenormal) + trans);
		Eigen::Vector3d newnormal = es.eigenvectors().col(smallesteigenidx)(Eigen::all, 0).real();
		std::cout << newplanepoint << "\n";
		// reset V2 and test if global transformation works
		// V2 = Eigen::MatrixXd(V1);
		//ICPAligner::applyRigidTransform(V2, Eigen::Matrix3d(Eigen::AngleAxisd((5.0 / 180.0) * EIGEN_PI, Eigen::Vector3d(1.0, 2.0, 0.0))), Eigen::Vector3d(2.0, 0.0, 0.0));
		// apply rotation and translation calculated using icp
		//ICPAligner::applyRigidTransform(V2, reflectedrot, trans, true);
		std::cout << "Optimal rotation:\n" << rot << "\n";
		std::cout << "Optimal translation:\n" << trans << "\n";

		Eigen::MatrixXd planeV1(4, 3);
		Eigen::MatrixXi planeF1(2, 3);
		Eigen::MatrixXd planeC(planeV1.rows(), 3);
		planeF1 << 0, 1, 2, 3, 2, 1;

		// define plane points
		Eigen::Vector3d perp1 = Eigen::Vector3d(newnormal(2), newnormal(0), newnormal(1)).cross(newnormal);
		Eigen::Vector3d perp2 = perp1.cross(newnormal);
		planeV1(0, Eigen::all) = perp1 * 100 + perp2 * 100;
		planeV1(1, Eigen::all) = perp1 * -100 + perp2 * 100;
		planeV1(2, Eigen::all) = perp1 * 100 + perp2 * -100;
		planeV1(3, Eigen::all) = perp1 * -100 + perp2 * -100;
		//planeV1 << 10.0, 0.0, 0.0, -10.0, 0.0, 0.0, 10.0, -10.0, 0.0, -10.0, -10.0, 0.0;
		std::cout << planeF1 << "\n";
		std::cout << planeV1;

		V = Eigen::MatrixXd(V1.rows() + V2.rows() + planeV1.rows(), V1.cols());
		F = Eigen::MatrixXi(F1.rows() + F2.rows() + planeF1.rows(), F1.cols());
		C = Eigen::MatrixXd(C1.rows() + C2.rows() + planeC.rows(), C1.cols());

		V << V1, V2, planeV1;
		F << F1, (F2.array() + V1.rows()), (planeF1.array() + V1.rows()+V2.rows());
		planeC << Eigen::RowVector3d(0.0, 1.0, 1.0).replicate(4, 1);
		C << C1, C2, planeC;
		
		// Plot the mesh
		igl::opengl::glfw::Viewer viewer2;
		viewer2.data().set_mesh(V, F);
		viewer2.data().set_colors(C);
		viewer2.launch();
	}
	catch (const std::exception& ex)
	{
		std::cerr << "Error: " << ex.what() << "\n";
	}
}