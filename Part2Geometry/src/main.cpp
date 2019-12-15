#include <iostream>
#include<igl/opengl/glfw/Viewer.h>
#include<igl/readOBJ.h>
#include <icp.h>
#include <Eigen/Geometry>



int main(int argc, char *argv[])
{
	Eigen::MatrixXd V1;
	Eigen::MatrixXi F1;
	igl::readOBJ("assets/models/RD-01/16021_OnyxCeph3_Export_OK-A.obj", V1, F1);
	Eigen::MatrixXd V2(V1);
	Eigen::MatrixXi F2(F1);
	ICPAligner::applyRigidTransform(V2, Eigen::Matrix3d(Eigen::AngleAxisd((30.0 / 180.0) * EIGEN_PI, Eigen::Vector3d(0.0, 1.0, 0.0))), Eigen::Vector3d(8.0, 0.0, 0.0));

	Eigen::MatrixXd N1;
	Eigen::MatrixXd N2;

	igl::per_vertex_normals(V1, F1, N1);
	igl::per_vertex_normals(V2, F2, N2);

	Eigen::Matrix3d rot;
	Eigen::Vector3d trans;

	ICPAligner icp(V1);
	icp.align(rot, trans, V2, V1, N1, N2, 3.0, 0.75, 1e-5);

	Eigen::MatrixXd V(V1.rows() + V2.rows(), V1.cols());
	Eigen::MatrixXi F(F1.rows() + F2.rows(), F1.cols());

	V << V1, V2;
	F << F1, (F2.array() + V1.rows());
	
	// Plot the mesh
	/*igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(V, F);
	viewer.launch();*/
}