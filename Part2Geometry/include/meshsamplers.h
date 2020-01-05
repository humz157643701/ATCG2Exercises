#ifndef _MESH_SAMPLERS_H_
#define _MESH_SAMPLERS_H_
#include <Eigen/Dense>
#include <vector>
#include <mesh.h>
namespace MeshSamplers
{
	// default sampler

	struct PassthroughSampler
	{
		static void sampleMeshPoints(
			const Mesh& mesh,
			Eigen::MatrixXd& sampled_points,
			Eigen::MatrixXd& sampled_normals
		)
		{
			sampled_points = mesh.vertices();
			sampled_normals = mesh.normals();
		}
	};
}

// advanced samplers

#include <integral_invariant_signatures.h>
#include <mesh_saliency.h>

#endif