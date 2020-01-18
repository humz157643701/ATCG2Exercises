#include <iostream>
#include<igl/opengl/glfw/Viewer.h>
#include<igl/readOBJ.h>
#include <Eigen/Geometry>
#include <chrono>
#include <mesh.h>
#include <igl/jet.h>
#include <fstream>
#include <tooth_segmentation.h>

int main(int argc, char* argv[])
{
	try
	{
		Eigen::MatrixXd V1;
		Eigen::MatrixXi F1;
		Eigen::MatrixXd N1;
		std::cout << "--- Loading meshes...\n";
		std::string model = "assets/models/RD-01/16021_OnyxCeph3_Export_OK-A.obj";
		igl::readOBJ(model, V1, F1);
		igl::per_vertex_normals(V1, F1, N1);
		Mesh original(V1, N1, F1);

		std::vector<Mesh> teeth;
		ToothSegmentation::segmentTeethFromMesh(original,
			Eigen::Vector3d{ 0.0, -1.0, 0.0 },
			Eigen::Vector3d{ -1.0, 0.0, 0.0 },
			teeth,
			ToothSegmentation::CuspDetectionParams{
				0.4,
				0.75,
				0.9,
				1000,
				0.5,
				0.0002,
				100,
				2.5
			},
			true
			);
	}
	catch (const std::exception & ex)
	{
		std::cerr << "Error: " << ex.what() << "\n";
	}
}