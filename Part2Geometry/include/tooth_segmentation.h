#ifndef _TOOTH_SEGMENTATION_H_
#define _TOOTH_SEGMENTATION_H_
#include <mesh.h>
#include <Eigen/Dense>
#include <memory>

class ToothSegmentation
{
public:
	// ALGORITHM OVERVIEW
	// 1. find cusps features
	// 2. fit plane through found cusps
	// 3. slice mesh below teeth
	// 4. identify per-tooth features and place boundary conditions
	// 5. solve for harmonic field
	// 6. cut out individual tooth meshes

	static void segmentTeethFromMesh(const Mesh& mesh, const Eigen::Vector3d& mesh_up, const Eigen::Vector3d& mesh_right, std::size_t expected_toothcount, std::vector<Mesh>& tooth_meshes);

private:
	static void computeCusps(const Mesh& mesh, std::size_t expected_toothcount, Eigen::VectorXi& featureIndices);
};


#endif