#ifndef _MESH_SAMPLERS_H_
#define _MESH_SAMPLERS_H_
#include <Eigen/Dense>

namespace MeshSamplers
{
	// default sampler

	struct PassthroughSampler
	{
		static void sampleMeshPoints(
			const Eigen::MatrixXd& vertices,
			const Eigen::MatrixXd& normals,
			const Eigen::MatrixXi& faces,
			Eigen::MatrixXd& sampled_points,
			Eigen::MatrixXd& sampled_normals
		)
		{
			sampled_points = vertices;
			sampled_normals = normals;
		}
	};


	// advanced samplers

	//#include <mesh_saliency.h>
	//#include <integral_invariant_signatures.h>
	//#include <salient_features_for_partial_shape_matching_and_similarity.h>


}

#endif