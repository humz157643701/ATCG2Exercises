#include <iostream>
#include<igl/opengl/glfw/Viewer.h>
#include<igl/readOBJ.h>
#include<igl/readOFF.h>
#include <symmetry.h>
#include <Eigen/Geometry>
#include <meshsamplers.h>
#include <chrono>
#include <mesh.h>
#include <igl/jet.h>
#include <fstream>
#include <planefitter.h>
#include <curvefitter.h>
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

int main(int argc, char *argv[])
{
	try
	{
		Eigen::MatrixXd V1;
		Eigen::MatrixXi F1;
		Eigen::MatrixXd N1;
		std::cout << "--- Loading meshes...\n";
		std::string model = "assets/models/RD-01/16021_OnyxCeph3_Export_OK-A.obj";
		std::string vxm = "assets/models/RD-01/16021_OnyxCeph3_Export_OK-A.obj_256.obj";
		//std::string model = "assets/models/head.obj";
		//std::string vxm = "assets/models/head.obj_256.obj";
		igl::readOBJ(model, V1, F1);
		//igl::readOBJ("assets/models/head.obj", V1, F1);
		//igl::readOBJ("assets/models/dino.obj", V1, F1);
		//igl::readOFF("assets/models/bumpy.off", V1, F1);
		igl::per_vertex_normals(V1, F1, N1);
		// build mesh
		std::cout << "--- Building mesh data structure...\n";
	
		V1 *= Eigen::MatrixXd(Eigen::AngleAxis<double>(M_PI, Eigen::Vector3d{ 0,0,1 }).toRotationMatrix());
		Mesh mesh(V1, N1, F1);
		std::vector<size_t> fpoints{ 228824,
									20453,

									9335  ,
									158157,
									165224,
									163144,

									3653  ,
									155934,

									154708,
									465	  ,

									157325,
									4649  ,

									156244,
									1727  ,

									154437,
									154323,

									154374,
									281	  ,

									6753  ,
									158652,

									158728,
									160864,

									2179  ,
									2224  ,

									5775  ,
									160101,

									155337,
									4017  ,
									5812  ,
									8910  ,

									16741 ,
									174894,
									160439,
									170818 };

		Eigen::MatrixXd features{ fpoints.size(),3 };
		for (int x = 0; x < fpoints.size(); ++x)
		{
			features.row(x) = mesh.vertices().row(fpoints[x]);// << std::cos(M_PI * (x / 10.0) - M_PI) * 5, 4.5 + ((std::rand() % 100) / 100.0), std::sin(M_PI * (x / 10.0)) * 5;
		}
		std::cout << features << std::endl;
		auto planeres = PlaneFitter::fitPlane(features);
		std::cout << "normal: " << planeres.first << std::endl;
		auto plane = createPlane(planeres.first, planeres.second, { 1, 1,0 });
		auto curveres = CurveFitter::fitCurve(features);
		std::cout << "curvemain: " << curveres << std::endl;
		auto res = CurveFitter::segmentFeatures(features, curveres, planeres.second(1),mesh);
		igl::opengl::glfw::Viewer viewer2;
		auto m1 = viewer2.append_mesh();
		auto m2 = viewer2.append_mesh();

		viewer2.data(m1).set_mesh(plane.Vs, plane.Fs);
		viewer2.data(m1).set_colors(plane.Cs);
		viewer2.data(m2).set_mesh(mesh.vertices(), mesh.faces());
		viewer2.data(m2).set_colors(plane.Cs);// Eigen::RowVector3d(1, 0, 1));
		viewer2.data().add_points(features, Eigen::RowVector3d(1, 0, 0));
		viewer2.data().add_points(*res.begin(), Eigen::RowVector3d(0, 1, 0));
		viewer2.data().add_points(res.back(), Eigen::RowVector3d(1, 1, 0));
		viewer2.data().point_size = 15;
		//for (size_t x = 0; x < res.size() / 20; x += 20)
		//{
			//viewer2.data().add_edges(res[0], res[1], Eigen::RowVector3d(0, 1, 0));
		//}
		viewer2.launch();

		//Task2
		{
		

		std::cout << "--- building voxel data structures...\n";
		//teeth factor 0.310623
		// head factor 0.001771
		MeshSamplers::IntegralInvariantSignaturesSampler integralsampler{ vxm,0.310623, 0.02, true };


		/////// MESH SALIENCY TEST STUFF /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		Eigen::Vector3d origvec = { 1,1,0.3 };

		SymmetryDetector<MeshSamplers::MeshSaliencySampler> symmd(ICPParams{ (mesh.vertices().colwise().maxCoeff() - mesh.vertices().colwise().minCoeff()).norm() * 0.3, 0.4, 1e-2, 1e-4, 150 }, MeshSamplers::MeshSaliencySampler(0.0010, 1, 5, false, 0.0085, ScaleType::DOUBLE_SIGMA_EVERY_SCALE, true));
		auto res = symmd.findMainSymmetryPlane(mesh, origvec.normalized());

		std::cout << "############################# MeshSaliencySampler Symmetry found! ###############################\n";
		std::cout << "Prefiltering: " << res.time_for_prefiltering << "s\n";
		std::cout << "ICP alignment: " << res.time_for_alignment << "s\n";
		std::cout << "Symmetry plane calculation: " << res.time_for_sym_calculation << "s\n";
		std::cout << "Time total: " << res.time_total << "s\n";

		{
			//Eigen::MatrixXd V2(V1);
			//Eigen::MatrixXi F2(F1);
			////ICPAligner::applyRigidTransform(V2, res.rot, res.trans);
			//double d_originplanedistance = 0;
			//Eigen::Matrix3d reflectionmatrix;
			//reflectionmatrix << 1, 0, 0, 0, 1, 0, 0, 0, 1;
			//reflectionmatrix = reflectionmatrix - 2 * (res.normal * res.normal.transpose());
			//d_originplanedistance = (Eigen::Vector3d(0, 0, 0) - res.center).dot(res.center);
			//
			//for (size_t x = 0; x < V2.rows(); ++x)
			//{
			//	Eigen::Vector3d refl = reflectionmatrix * V2(x, Eigen::all).transpose() + 2 * d_originplanedistance * res.normal;


			//	V2(x, Eigen::all) = refl;
			//}

			//
			//Eigen::MatrixXd C1(V1.rows(), 3);
			//Eigen::MatrixXd C2(V2.rows(), 3);
			//C1 << Eigen::RowVector3d(1.0, 1.0, 1.0).replicate(V1.rows(), 1);
			//C2 << Eigen::RowVector3d(1.0, 0.0, 0.0).replicate(V2.rows(), 1);

			//auto planeRes = createPlane(res.normal);

			//Eigen::MatrixXd V(V1.rows() + V2.rows() + planeRes.Vs.rows(), V1.cols());
			//Eigen::MatrixXi F(F1.rows() + F2.rows() + planeRes.Fs.rows(), F1.cols());
			//Eigen::MatrixXd C(C1.rows() + C2.rows() + planeRes.Cs.rows(), C1.cols());

			//V << V1, V2, planeRes.Vs;
			//F << F1, (F2.array() + V1.rows()), (planeRes.Fs.array() + V1.rows() + V2.rows());
			//C << C1, C2, planeRes.Cs;


			//auto res2 = calcSymmetryPlane(V1, F1, { 0,1,0 });

			//Eigen::MatrixXd V3(V1);
			//Eigen::MatrixXi F3(F1);
			////ICPAligner::applyRigidTransform(V3, res2.rot, res2.trans);
			//d_originplanedistance = 0;
			//reflectionmatrix << 1, 0, 0, 0, 1, 0, 0, 0, 1;
			//reflectionmatrix = reflectionmatrix - 2 * (res2.normal * res2.normal.transpose());
			//d_originplanedistance = (Eigen::Vector3d(0, 0, 0) - res2.center).dot(res2.center);

			//for (size_t x = 0; x < V3.rows(); ++x)
			//{
			//	Eigen::Vector3d refl = reflectionmatrix * V3(x, Eigen::all).transpose() + 2 * d_originplanedistance * res2.normal;


			//	V3(x, Eigen::all) = refl;
			//}
			//
			//Eigen::MatrixXd C3(V3.rows(), 3);
			//C3 << Eigen::RowVector3d(0.0, 1.0, 0.0).replicate(V3.rows(), 1);

			//auto planeRes2 = createPlane(res2.normal);

			//auto res3 = calcSymmetryPlane(V1, F1, { 0,0,1 });

			//Eigen::MatrixXd V4(V1);
			//Eigen::MatrixXi F4(F1);
			////ICPAligner::applyRigidTransform(V4, res3.rot, res3.trans);
			//d_originplanedistance = 0;
			//reflectionmatrix << 1, 0, 0, 0, 1, 0, 0, 0, 1;
			//reflectionmatrix = reflectionmatrix - 2 * (res3.normal * res3.normal.transpose());
			//d_originplanedistance = (Eigen::Vector3d(0, 0, 0) - res3.center).dot(res3.center);

			//for (size_t x = 0; x < V4.rows(); ++x)
			//{
			//	Eigen::Vector3d refl = reflectionmatrix * V4(x, Eigen::all).transpose() + 2 * d_originplanedistance * res3.normal;


			//	V4(x, Eigen::all) = refl;
			//}
			//Eigen::MatrixXd C4(V4.rows(), 3);
			//C4 << Eigen::RowVector3d(0.0, 0.0, 1.0).replicate(V4.rows(), 1);

			//auto planeRes3 = createPlane(res2.normal);
		}
		Eigen::MatrixXd C1(V1.rows(), 3);
		Eigen::MatrixXd C2(V1.rows(), 3);
		C1 << Eigen::RowVector3d(1.0, 1.0, 1.0).replicate(V1.rows(), 1);
		C2 << Eigen::RowVector3d(1.0, 0.0, 1.0).replicate(V1.rows(), 1);
		auto pl_res = createPlane(res.normal, res.center, { 0,1,0 });
		Eigen::Vector3d center_of_mass{ V1.colwise().mean().transpose() };
		auto pl_init = createPlane(origvec, center_of_mass, { 1,1,0 });

		// Plot the mesh
		igl::opengl::glfw::Viewer viewer2;
		//auto m3 = viewer2.append_mesh();
		auto m1 = viewer2.append_mesh();
		auto m2 = viewer2.append_mesh();
		auto m3 = viewer2.append_mesh();
		auto m4 = viewer2.append_mesh();

		viewer2.data(m1).set_mesh(V1, F1);
		viewer2.data(m1).set_colors(C1);
		Eigen::MatrixXd reflV1(V1);
		//std::cout << "#############################\n Symm Result -\n normal: \n" << res.normal << "\n center: \n" << res.center << "\n  " << "\n#############################\n";


		reflectAlongPlane(res, reflV1);

		viewer2.data(m2).set_mesh(reflV1, F1);
		viewer2.data(m2).set_colors(C2);
		viewer2.data(m3).set_mesh(pl_init.Vs, pl_init.Fs);
		viewer2.data(m3).set_colors(pl_init.Cs);
		viewer2.data(m4).set_mesh(pl_res.Vs, pl_res.Fs);
		viewer2.data(m4).set_colors(pl_res.Cs);
		viewer2.data().add_points(Eigen::RowVector3d(0, 0, 0), Eigen::RowVector3d(1, 0, 0));
		viewer2.launch();

		SymmetryDetector<MeshSamplers::IntegralInvariantSignaturesSampler> symmd2(ICPParams{ (mesh.vertices().colwise().maxCoeff() - mesh.vertices().colwise().minCoeff()).norm() * 0.3, 0.4, 1e-2, 1e-4, 150 }, integralsampler);//0.310623 });
		res = symmd2.findMainSymmetryPlane(mesh, origvec.normalized());

		std::cout << "############################# IntegralInvariantSignaturesSampler Symmetry found! ###############################\n";
		std::cout << "Prefiltering: " << res.time_for_prefiltering << "s\n";
		std::cout << "ICP alignment: " << res.time_for_alignment << "s\n";
		std::cout << "Symmetry plane calculation: " << res.time_for_sym_calculation << "s\n";
		std::cout << "Time total: " << res.time_total << "s\n";

		C1 = Eigen::MatrixXd(V1.rows(), 3);
		C2 = Eigen::MatrixXd(V1.rows(), 3);
		C1 << Eigen::RowVector3d(1.0, 1.0, 1.0).replicate(V1.rows(), 1);
		C2 << Eigen::RowVector3d(1.0, 0.0, 1.0).replicate(V1.rows(), 1);
		pl_res = createPlane(res.normal, res.center, { 0,1,0 });
		center_of_mass = Eigen::Vector3d{ V1.colwise().mean().transpose() };
		pl_init = createPlane(origvec, center_of_mass, { 1,1,0 });

		// Plot the mesh
		viewer2 = igl::opengl::glfw::Viewer();
		//auto m3 = viewer2.append_mesh();
		m1 = viewer2.append_mesh();
		m2 = viewer2.append_mesh();
		m3 = viewer2.append_mesh();
		m4 = viewer2.append_mesh();

		viewer2.data(m1).set_mesh(V1, F1);
		viewer2.data(m1).set_colors(C1);
		reflV1 = Eigen::MatrixXd(V1);

		reflectAlongPlane(res, reflV1);

		viewer2.data(m2).set_mesh(reflV1, F1);
		viewer2.data(m2).set_colors(C2);
		viewer2.data(m3).set_mesh(pl_init.Vs, pl_init.Fs);
		viewer2.data(m3).set_colors(pl_init.Cs);
		viewer2.data(m4).set_mesh(pl_res.Vs, pl_res.Fs);
		viewer2.data(m4).set_colors(pl_res.Cs);
		viewer2.data().add_points(Eigen::RowVector3d(0, 0, 0), Eigen::RowVector3d(1, 0, 0));
		viewer2.launch();

		SymmetryDetector<MeshSamplers::PassthroughSampler> symmd3(ICPParams{ (mesh.vertices().colwise().maxCoeff() - mesh.vertices().colwise().minCoeff()).norm() * 0.3, 0.4, 1e-2, 1e-4, 150 }, {});//0.310623 });
		res = symmd3.findMainSymmetryPlane(mesh, origvec.normalized());

		std::cout << "############################# PassThroughSampler Symmetry found! ###############################\n";
		std::cout << "Prefiltering: " << res.time_for_prefiltering << "s\n";
		std::cout << "ICP alignment: " << res.time_for_alignment << "s\n";
		std::cout << "Symmetry plane calculation: " << res.time_for_sym_calculation << "s\n";
		std::cout << "Time total: " << res.time_total << "s\n";

		C1 = Eigen::MatrixXd(V1.rows(), 3);
		C2 = Eigen::MatrixXd(V1.rows(), 3);
		C1 << Eigen::RowVector3d(1.0, 1.0, 1.0).replicate(V1.rows(), 1);
		C2 << Eigen::RowVector3d(1.0, 0.0, 1.0).replicate(V1.rows(), 1);
		pl_res = createPlane(res.normal, res.center, { 0,1,0 });
		center_of_mass = Eigen::Vector3d{ V1.colwise().mean().transpose() };
		pl_init = createPlane(origvec, center_of_mass, { 1,1,0 });

		// Plot the mesh
		viewer2 = igl::opengl::glfw::Viewer();
		//auto m3 = viewer2.append_mesh();
		m1 = viewer2.append_mesh();
		m2 = viewer2.append_mesh();
		m3 = viewer2.append_mesh();
		m4 = viewer2.append_mesh();

		viewer2.data(m1).set_mesh(V1, F1);
		viewer2.data(m1).set_colors(C1);
		reflV1 = Eigen::MatrixXd(V1);

		reflectAlongPlane(res, reflV1);

		viewer2.data(m2).set_mesh(reflV1, F1);
		viewer2.data(m2).set_colors(C2);
		viewer2.data(m3).set_mesh(pl_init.Vs, pl_init.Fs);
		viewer2.data(m3).set_colors(pl_init.Cs);
		viewer2.data(m4).set_mesh(pl_res.Vs, pl_res.Fs);
		viewer2.data(m4).set_colors(pl_res.Cs);
		viewer2.data().add_points(Eigen::RowVector3d(0, 0, 0), Eigen::RowVector3d(1, 0, 0));
		viewer2.launch();
	}
	}
	catch (const std::exception& ex)
	{
		std::cerr << "Error: " << ex.what() << "\n";
	}
}