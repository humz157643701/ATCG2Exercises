#ifndef _MESH_SALIENCY_H_
#define _MESH_SALIENCY_H_

namespace MeshSamplers
{
	// default sampler

	struct MeshSaliencySampler
	{
		static void sampleMeshPoints(
			const Eigen::MatrixXd& vertices,
			const Eigen::MatrixXd& normals,
			const Eigen::MatrixXi& faces,
			Eigen::MatrixXd& sampled_points,
			Eigen::MatrixXd& sampled_normals
		)
		{
			
		}
	};
}

#endif