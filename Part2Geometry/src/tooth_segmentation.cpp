#include "..\include\tooth_segmentation.h"

void ToothSegmentation::segmentTeethFromMesh(const Mesh& mesh, const Eigen::Vector3d& mesh_up, const Eigen::Vector3d& mesh_right, std::size_t expected_toothcount, std::vector<Mesh>& tooth_meshes)
{
	// transform mesh to canonical pose
	Mesh working_mesh(mesh);
	Eigen::Matrix3d crot;
	crot << mesh_right.normalized(), mesh_up.normalized(), mesh_right.cross(mesh_up).normalized();
	crot = crot.transpose();
	Eigen::Vector3d ctrans = mesh.vertices().colwise().mean().transpose();

	working_mesh.vertices() *= crot.transpose();
	working_mesh.vertices().rowwise() += ctrans.transpose();

	// compute cusp features
	Eigen::VectorXi cuspIndices;
	computeCusps(working_mesh, expected_toothcount, cuspIndices);

	// fit plane

	// cut teeth

	// fit curve

	// assign features to teeth

	// harmonic field stuff
}

void ToothSegmentation::computeCusps(const Mesh& mesh, std::size_t expected_toothcount, Eigen::VectorXi& featureIndices)
{
}
