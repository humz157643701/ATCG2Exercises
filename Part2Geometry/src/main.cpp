#include <iostream>
#include<igl/opengl/glfw/Viewer.h>
#include<igl/readOBJ.h>
#include <icp.h>
#include <Eigen/Geometry>



int main(int argc, char *argv[])
{
	try
	{
		Eigen::MatrixXd V1;
		Eigen::MatrixXi F1;
		std::cout << "Loading meshes...\n";
		igl::readOBJ("assets/models/RD-01/16021_OnyxCeph3_Export_OK-A.obj", V1, F1);
		Eigen::MatrixXd V2(V1);
		Eigen::MatrixXi F2(F1);
		ICPAligner::applyRigidTransform(V2, Eigen::Matrix3d(Eigen::AngleAxisd((5.0 / 180.0) * EIGEN_PI, Eigen::Vector3d(1.0, 2.0, 0.0))), Eigen::Vector3d(2.0, 0.0, 0.0));

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
		viewer.data().set_mesh(V, F);
		viewer.data().set_colors(C);
		viewer.launch();

		ICPAligner icp(V1);
		icp.align(rot, trans, V2, V1, N1, N2, 10.0, 0.75, 1e-2, 200);

		V << V1, V2;

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