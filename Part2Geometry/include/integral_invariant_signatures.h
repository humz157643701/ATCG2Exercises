#ifndef _INTEGRAL_INVARIANT_SIGNATURES_H_
#define _INTEGRAL_INVARIANT_SIGNATURES_H_
namespace MeshSamplers
{
	struct IntegralInvariantSignaturesSampler
	{
		std::string voxelObj;
		IntegralInvariantSignaturesSampler(std::string pathToVoxel): voxelObj(pathToVoxel)
		{

		}

		static Eigen::MatrixXd readVoxelOBJ(std::string path)
		{
			std::ifstream voxelrep(path);
			std::string buf;
			size_t n = 0;
			while (std::getline(voxelrep, buf))
				n++;
			Eigen::MatrixXd Vx1(n, 3);
			voxelrep.clear();
			voxelrep.seekg(0, std::ios::beg);

			for (size_t i = 0; i < n; ++i)
			{
				int waste;
				double x, y, z;
				std::getline(voxelrep, buf);
				std::istringstream iss(buf);
				std::vector<std::string> results(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
				Vx1.row(i) << std::stod(results[1]), std::stod(results[2]), std::stod(results[3]);
			}
			return Vx1;
		}

		void sampleMeshPoints(
			const Mesh& mesh,
			Eigen::MatrixXd& sampled_points,
			Eigen::MatrixXd& sampled_normals
		)
		{
			double cellsize = 0.05;
			double kernelsize = 0.15;
			auto Voxels = readVoxelOBJ(this->voxelObj);
			igl::opengl::glfw::Viewer view;
			view.data().add_points(Voxels, Eigen::RowVector3d(1, 0, 0));
			view.launch();

			sampled_points = mesh.vertices();
			sampled_normals = mesh.normals();
		}
	};
}
#endif