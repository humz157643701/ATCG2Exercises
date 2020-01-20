#include <iostream>
#include<igl/opengl/glfw/Viewer.h>
#include<igl/readOBJ.h>
#include <Eigen/Geometry>
#include <chrono>
#include <mesh.h>
#include <igl/jet.h>
#include <fstream>
#include <tooth_segmentation.h>
#include <cmath>

int main(int argc, char* argv[])
{
	try
	{
		Eigen::MatrixXd V1;
		Eigen::MatrixXi F1;
		Eigen::MatrixXd N1;
		std::cout << "--- Loading meshes...\n";
		std::string model = "assets/models/RD-04/14645_OnyxCeph3_Export_OK-A.obj";
		igl::readOBJ(model, V1, F1);
		igl::per_vertex_normals(V1, F1, N1);
		Mesh original(V1, N1, F1);

		std::vector<Mesh> teeth;
		ToothSegmentation::segmentTeethFromMesh(original,
			Eigen::Vector3d{ 0.0, -1.0, 0.0 },
			Eigen::Vector3d{ -1.0, 0.0, 0.0 },
			teeth,
			ToothSegmentation::CuspDetectionParams{
				0.5, // curvature <-> height
				0.8, // curvature exponent
				0.9, // height exponent
				0.5, // min feature height
				0.0003, // smoothing step size
				50, // smoothing steps
				2.0, // max zscore
				0.5, // fraction of local maxima to be considered for mean shift
				0.012, // mean shift window size
				1e-4, // min total shift before convergence
				1000, // max mean shift iterations
				0.03 // feature collapse distance
			},
			true
			);
	}
	catch (const std::exception & ex)
	{
		std::cerr << "Error: " << ex.what() << "\n";
	}
}