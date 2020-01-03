#ifndef _INTEGRAL_INVARIANT_SIGNATURES_H_
#define _INTEGRAL_INVARIANT_SIGNATURES_H_
namespace MeshSamplers
{
	struct IntegralInvariantSignaturesSampler
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
}
#endif