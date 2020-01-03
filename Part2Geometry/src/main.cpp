#include <iostream>
#include<igl/opengl/glfw/Viewer.h>
#include<igl/readOBJ.h>
#include <symmetry.h>
#include <Eigen/Geometry>
#include <meshsamplers.h>

struct ResultData
{
	Eigen::Vector3d center;
	Eigen::Vector3d normal;
	Eigen::Vector3d trans;
	Eigen::Matrix3d rot;
};

struct Plane
{
	Eigen::MatrixXd Vs;
	Eigen::MatrixXi Fs;
	Eigen::MatrixXd Cs;
};

Plane createPlane(Eigen::Vector3d &normal)
{
	Eigen::MatrixXd planeV1(4, 3);
	Eigen::MatrixXi planeF1(2, 3);
	Eigen::MatrixXd planeC(planeV1.rows(), 3);
	planeC << Eigen::RowVector3d(0.0, 1.0, 1.0).replicate(4, 1);
	planeF1 << 0, 1, 2, 3, 2, 1;

	// define plane points
	Eigen::Vector3d perp1 = Eigen::Vector3d(normal(2), normal(0), normal(1)).cross(normal);
	Eigen::Vector3d perp2 = perp1.cross(normal);
	planeV1(0, Eigen::all) = perp1 * 50 + perp2 * 50;
	planeV1(1, Eigen::all) = perp1 * -50 + perp2 * 50;
	planeV1(2, Eigen::all) = perp1 * 50 + perp2 * -50;
	planeV1(3, Eigen::all) = perp1 * -50 + perp2 * -50;
	//planeV1 << 10.0, 0.0, 0.0, -10.0, 0.0, 0.0, 10.0, -10.0, 0.0, -10.0, -10.0, 0.0;
	std::cout << planeF1 << "\n";
	std::cout << planeV1;

	return { planeV1, planeF1, planeC };
}

//ResultData calcSymmetryPlane(Eigen::MatrixXd& V1, Eigen::MatrixXi &F1, Eigen::Vector3d initialPlaneNormal)
//{
//	Eigen::MatrixXd V2(V1);
//	Eigen::MatrixXi F2(F1);
//	//ICPAligner::applyRigidTransform(V2, Eigen::Matrix3d(Eigen::AngleAxisd((5.0 / 180.0) * EIGEN_PI, Eigen::Vector3d(1.0, 2.0, 0.0))), Eigen::Vector3d(2.0, 0.0, 0.0));
//
//	Eigen::MatrixXd N1;
//	Eigen::MatrixXd N2;
//
//	igl::per_vertex_normals(V1, F1, N1);
//	igl::per_vertex_normals(V2, F2, N2);
//
//	
//	Eigen::Matrix3d rot;
//	Eigen::Vector3d trans;
//
//	
//	//V << V1, V2;
//	//F << F1, (F2.array() + V1.rows());
//	//C << C1, C2;
//
//
//
//	// Plot the mesh
//	//igl::opengl::glfw::Viewer viewer;
//	//viewer.data().set_mesh(V, F);
//	//viewer.data().set_colors(C);
//	//viewer.launch();
//
//	Eigen::Vector3d centerofmass{ 0,0,0 };
//	Eigen::Vector3d planenormal(initialPlaneNormal);
//	double d_originplanedistance = 0;
//	Eigen::Matrix3d reflectionmatrix;
//	reflectionmatrix << 1, 0, 0, 0, 1, 0, 0, 0, 1;
//	std::cout << reflectionmatrix << "\n";
//	reflectionmatrix = reflectionmatrix - 2 * (planenormal * planenormal.transpose());
//	std::cout << reflectionmatrix << "\n";
//
//	for (size_t x = 0; x < V1.rows(); ++x)
//	{
//		centerofmass(0) += V1(x, 0);
//		centerofmass(1) += V1(x, 1);
//		centerofmass(2) += V1(x, 2);
//	}
//	centerofmass /= V1.rows();
//	d_originplanedistance = (Eigen::Vector3d(0, 0, 0) - centerofmass).dot(planenormal);
//	std::cout << centerofmass << "\n";
//	std::cout << d_originplanedistance;
//
//	for (size_t x = 0; x < V2.rows(); ++x)
//	{
//		Eigen::Vector3d refl = reflectionmatrix * V2(x, Eigen::all).transpose() + 2 * d_originplanedistance * planenormal;
//
//
//		V2(x, Eigen::all) = refl;
//	}
//
//	for (size_t x = 0; x < N2.rows(); ++x)
//	{
//		Eigen::Vector3d refl = reflectionmatrix * N2(x, Eigen::all).transpose() + 2 * d_originplanedistance * planenormal;
//
//
//		N2(x, Eigen::all) = refl;
//	}
//
//	ICPAligner icp(V1);
//	// align transforms the query set (V2) and returns optimal translation and rotation
//	icp.align(rot, trans, V2, V1, N1, N2, 50.0, 0.2, 1e-2, 1e-4, 35);
//
//	auto reflectedrot = reflectionmatrix * rot.transpose();
//	Eigen::EigenSolver<Eigen::MatrixXd> es(reflectedrot);
//	int smallesteigenidx = 0;
//	double smallesteigendiff = 1;
//	for (size_t x = 0; x < es.eigenvalues().rows(); ++x)
//	{
//		auto diff = es.eigenvalues()(x, 0).real() + 1;
//		if (diff < smallesteigendiff)
//		{
//			smallesteigenidx = x;
//			smallesteigendiff = diff;
//		}
//	}
//	std::cout << es.eigenvalues()(smallesteigenidx).real() << " | " << smallesteigendiff << "\n";
//	std::cout << es.eigenvectors().col(smallesteigenidx)(Eigen::all, 0).real() << "\n";
//
//	// New symm plane point
//	Eigen::Vector3d newplanepoint = 0.5 * (rot * (2 * d_originplanedistance * planenormal) + trans);
//	Eigen::Vector3d newnormal = es.eigenvectors().col(smallesteigenidx)(Eigen::all, 0).real();
//	std::cout << newplanepoint << "\n";
//	// reset V2 and test if global transformation works
//	// V2 = Eigen::MatrixXd(V1);
//	//ICPAligner::applyRigidTransform(V2, Eigen::Matrix3d(Eigen::AngleAxisd((5.0 / 180.0) * EIGEN_PI, Eigen::Vector3d(1.0, 2.0, 0.0))), Eigen::Vector3d(2.0, 0.0, 0.0));
//	// apply rotation and translation calculated using icp
//	//ICPAligner::applyRigidTransform(V2, reflectedrot, trans, true);
//	std::cout << "Optimal rotation:\n" << rot << "\n";
//	std::cout << "Optimal translation:\n" << trans << "\n";
//
//	return { newplanepoint, newnormal, trans, rot };
//}


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
		Eigen::MatrixXd N1;
		std::cout << "Loading meshes...\n";
		igl::readOBJ("assets/models/RD-01/16021_OnyxCeph3_Export_OK-A.obj", V1, F1);
		igl::per_vertex_normals(V1, F1, N1);
		Eigen::Vector3d origvec = { 1,1,0 };
		SymmetryDetector<MeshSamplers::PassthroughSampler> symmd(ICPParams{ 50.0, 0.2, 1e-2, 1e-4, 50 }, MeshSamplers::PassthroughSampler{});
		auto res = symmd.findMainSymmetryPlane(V1, N1, F1, origvec.normalized());

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


		//// Plot the mesh
		//igl::opengl::glfw::Viewer viewer2;
		////auto m3 = viewer2.append_mesh();
		//viewer2.data().set_mesh(V, F);
		//viewer2.data().set_colors(C);
		//viewer2.launch();
		//viewer2 = igl::opengl::glfw::Viewer();
		//auto m1 = viewer2.append_mesh();
		//auto m2 = viewer2.append_mesh();
		//auto m3 = viewer2.append_mesh();

		//viewer2.data(m1).set_mesh(V3, F3);
		//viewer2.data(m1).set_colors(C3);
		//viewer2.data(m2).set_mesh(planeRes2.Vs, planeRes2.Fs);
		//viewer2.data(m2).set_colors(planeRes2.Cs);
		//viewer2.data(m3).set_mesh(V1, F1);
		//viewer2.data(m3).set_colors(C1);
		//viewer2.launch();
		//viewer2 = igl::opengl::glfw::Viewer();
		//m1 = viewer2.append_mesh();
		//m2 = viewer2.append_mesh();
		//m3 = viewer2.append_mesh();

		//viewer2.data(m1).set_mesh(V4, F4);
		//viewer2.data(m1).set_colors(C4);
		//viewer2.data(m2).set_mesh(planeRes3.Vs, planeRes3.Fs);
		//viewer2.data(m2).set_colors(planeRes3.Cs);
		//viewer2.data(m3).set_mesh(V1, F1);
		//viewer2.data(m3).set_colors(C1);
		//viewer2.launch();

	}
	catch (const std::exception& ex)
	{
		std::cerr << "Error: " << ex.what() << "\n";
	}
}