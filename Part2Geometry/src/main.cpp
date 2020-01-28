#include <iostream>
#include<igl/opengl/glfw/Viewer.h>
#include<igl/readOBJ.h>
#include <Eigen/Geometry>
#include <chrono>
#include <mesh.h>
#include <igl/jet.h>
#include <fstream>
#include <planefitter.h>
#include <curvefitter.h>
#include <tooth_segmentation.h>
#define _USE_MATH_DEFINES
#include <math.h>
struct Plane
{
	Eigen::MatrixXd Vs;
	Eigen::MatrixXi Fs;
	Eigen::MatrixXd Cs;
};

Plane createPlane(Eigen::Vector3d &normal, Eigen::Vector3d& trans, Eigen::RowVector3d color)
{
	normal.normalize();
	Eigen::MatrixXd planeV1(4, 3);
	Eigen::MatrixXi planeF1(2, 3);
	Eigen::MatrixXd planeC(planeV1.rows(), 3);
	planeC << color.replicate(4, 1);
	planeF1 << 0, 1, 2, 3, 2, 1;

	// define plane points
	Eigen::Vector3d perp1 = Eigen::Vector3d(normal(2), normal(0), normal(1)).cross(normal);
	Eigen::Vector3d perp2 = perp1.cross(normal);
	planeV1(0, Eigen::all) = perp1 * 20 + perp2 * 20;
	planeV1(1, Eigen::all) = perp1 * -20 + perp2 * 20;
	planeV1(2, Eigen::all) = perp1 * 20 + perp2 * -20;
	planeV1(3, Eigen::all) = perp1 * -20 + perp2 * -20;
	planeV1.rowwise() += (trans.dot(normal) * normal).transpose();
	return { planeV1, planeF1, planeC };
}

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
		// build mesh
		std::cout << "--- Building mesh data structure...\n";
		Mesh mesh(V1, N1, F1);



		std::vector<Mesh> teeth;
		ToothSegmentation::segmentTeethFromMesh(mesh,
			Eigen::Vector3d{ 0.0, -1.0, 0.0 },
			Eigen::Vector3d{ -1.0, 0.0, 0.0 },
			teeth,
			ToothSegmentation::CuspDetectionParams{
				0.6, // curvature <-> height
				0.85, // curvature exponent
				0.8, // height exponent
				0.5, // min feature height
				0.75, // fraction of local maxima to be considered for mean shift
				0.012, // mean shift window size
				1e-4, // min total shift before convergence
				1000, // max mean shift iterations
				0.003, // feature collapse distance
				0.01, // small feature detection window size
				0.2 // small feature detection threshold
			},
			ToothSegmentation::HarmonicFieldParams{
				1.0,	// w
				1.0,	// cvtr weight high
				0.001, // cvtr weight low
				0.1		// neg. cvtr threshold
			},
			ToothSegmentation::MeanCurvatureParams{
				0.00025, // smoothing step size
				50, // smoothing steps
				2.0 // max zscore
			},
			ToothSegmentation::ToothMeshExtractionParams{
				0.25, // even tooth threshold
				0.75 // odd tooth treshold
			},
			true
		);


	}
	catch (const std::exception& ex)
	{
		std::cerr << "Error: " << ex.what() << "\n";
	}
}