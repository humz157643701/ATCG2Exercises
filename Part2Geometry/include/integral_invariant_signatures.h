#ifndef _INTEGRAL_INVARIANT_SIGNATURES_H_
#define _INTEGRAL_INVARIANT_SIGNATURES_H_
namespace MeshSamplers
{
	struct IntegralInvariantSignaturesSampler
	{
		std::string voxelObj;
		double voxelScale;
		Eigen::MatrixXd voxels;
		double featureCountScale;
		bool visualize = false;
		IntegralInvariantSignaturesSampler(std::string pathToVoxel, double vxlScale, double featCntScale, bool visualizeEnable): voxelObj(pathToVoxel), voxelScale(vxlScale), featureCountScale(featCntScale), visualize(visualizeEnable)
		{
			this->voxels = readVoxelOBJ(this->voxelObj);
			//Voxels.rowwise() -= (Eigen::Vector3d{ .5, .5, .5 }).transpose();
			this->voxels = fillVoxelMatrix(this->voxels, 256);
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

		//Fill surface voxel matrix to make it dense
		static Eigen::MatrixXd fillVoxelMatrix(const Eigen::MatrixXd& src, int vxlDim)
		{
			std::unique_ptr<kdtree_t> m_kdtree;
			m_kdtree.reset(new kdtree_t(3, src));
			m_kdtree->index->buildIndex();
			std::vector<Eigen::Vector3d> additions;
			additions.reserve(vxlDim * vxlDim * vxlDim - (src.cols()*src.rows()));
			std::vector<Eigen::Vector3d> additionstmp;
			additionstmp.reserve(256);
			bool inside = false;
			bool lastinside = false;
			int cntr = 0;
			std::vector<Eigen::Index>   ret_index(1);
			std::vector<double> out_dist_sqr(1);
			out_dist_sqr[0] = 0.5;
			nanoflann::SearchParams params;
			for (size_t x = 0; x <= vxlDim; ++x)
			{
				for (size_t y = 0; y <= vxlDim; ++y)
				{
					for (size_t z = 0; z <= vxlDim; ++z)
					{
						double query_pt[3] = { x, y, z };
						
						m_kdtree->query(query_pt, 1, &ret_index[0], &out_dist_sqr[0]);
						if ( out_dist_sqr[0]<0.3&&!lastinside)
						{
							inside = !inside;
						}
						else if (out_dist_sqr[0] >= 0.3 && inside)
						{
							additionstmp.push_back({ float(x), float(y), float(z) });
							cntr++;
						}
						lastinside = out_dist_sqr[0] < 0.3;
						
					}
					if(!(inside && !lastinside))
						additions.insert(additions.end(), additionstmp.begin(), additionstmp.end());
					additionstmp.clear();
					inside = false;
					lastinside = false;
				}
			}
			std::cout << "num adds: " << cntr << std::endl;
			Eigen::MatrixXd newm(src.rows()+ additions.size(), 3);
			newm << src;
			for (size_t x = 0; x < additions.size(); ++x)
			{
				newm.row(x+src.rows()) = additions[x];
			}
			return newm;
		}

		void sampleMeshPoints(
			const Mesh& mesh,
			Eigen::MatrixXd& sampled_points,
			Eigen::MatrixXd& sampled_normals
		)
		{
			auto mins = mesh.vertices().colwise().minCoeff();

			auto& Voxels = this->voxels;
			//Align point cloud wish original mesh
			Voxels = Voxels * this->voxelScale;// *0.001771;
			Voxels.rowwise() += (mins- Voxels.colwise().minCoeff());


			std::unique_ptr<kdtree_t> m_kdtree;
			m_kdtree.reset(new kdtree_t(3, Voxels));
			m_kdtree->index->buildIndex();
			//Eigen::MatrixXd descriptors(mesh.vertices().rows(),1);

			//descriptor list with just value to vertex index
			std::vector<std::pair<size_t, Eigen::Index>> descr(mesh.vertices().rows());
			//temporary map storing feature value to feature count and vertex index
			std::map<size_t, std::pair<size_t, std::vector<Eigen::Index>>> descrHist;

			Eigen::MatrixXd cs(mesh.vertices().rows(),3);

			size_t maxFeatureVal = 0;
			for (size_t i = 0; i < mesh.vertices().rows();++i)
			{
				std::vector<std::pair<Eigen::Index, double>> ret;
				double pt[3]{ mesh.vertices()(i,0), mesh.vertices()(i,1), mesh.vertices()(i,2) };
				m_kdtree->index->radiusSearch(&pt[0], (this->voxelScale*5)*(this->voxelScale*5), ret, nanoflann::SearchParams(0, 0.0, false));
				
				//mesh color stuff
				descr[i] = { ret.size(), i };
				if (ret.size() > maxFeatureVal)
					maxFeatureVal = ret.size();

				//filling histogram stuff
				if (descrHist.count(ret.size()) == 0)
				{
					descrHist.insert({ ret.size(),{} });
				}
				descrHist[ret.size()].first += 1;
				descrHist[ret.size()].second.push_back(i);
			}
			//coloring, vector index is still identical to row index
			for (size_t i = 0; i < mesh.vertices().rows(); ++i)
			{
				double ratio = float(descr[i].first) / maxFeatureVal;
				cs.row(i) << Eigen::Vector3d(ratio, 0, 1 - ratio).transpose();
			}

			//Put map into vector to sort it
			std::vector<std::pair<size_t, std::pair<size_t, std::vector<Eigen::Index>>>> descrHist2;
			descrHist2.reserve(descrHist.size());
			for (auto kv : descrHist)
			{
				descrHist2.push_back({ kv.first, kv.second });
			}
			std::sort(descrHist2.begin(), descrHist2.end(), [](const std::pair<size_t, std::pair<size_t, std::vector<Eigen::Index>>>& a, const std::pair<size_t, std::pair<size_t, std::vector<Eigen::Index>>>& b) {
				return a.second.first> b.second.first;
			});

			//Go through first some % rarest features and add them to our output
			Eigen::MatrixXd memes(size_t(mesh.vertices().rows()*this->featureCountScale), 3);
			Eigen::MatrixXd memesNormal(size_t(mesh.vertices().rows() * this->featureCountScale), 3);
			size_t cntr = 0;
			for (auto d : descrHist2)
			{	
				for (auto i : d.second.second)
				{
					memes.row(cntr) << mesh.vertices().row(i);
					memesNormal.row(cntr++) << mesh.normals().row(i);
					if (cntr >= size_t(mesh.vertices().rows() * this->featureCountScale))
						break;
				}
				if (cntr >= size_t(mesh.vertices().rows() * this->featureCountScale))
					break;
			}

			if (this->visualize)
			{
				igl::opengl::glfw::Viewer view;
				view.data().set_mesh(mesh.vertices(), mesh.faces());
				view.data().set_colors(cs);
				view.data().add_points(Eigen::RowVector3d(0, 0, 0), Eigen::RowVector3d(1, 1, 0));
				view.data().add_points(memes, Eigen::RowVector3d(0, 1, 0));
				view.data().point_size = 5.0;
				view.launch();
			}
			sampled_points.resize(memes.rows(), memes.cols());
			sampled_points << memes;
			sampled_normals.resize(memes.rows(), memes.cols());
			sampled_normals << memesNormal;
			//sampled_points = mesh.vertices();
			//sampled_normals = mesh.normals();
		}
	};
}
#endif