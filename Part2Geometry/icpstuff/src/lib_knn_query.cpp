#include "lib_knn_query.h"
#include <nanoflann.hpp>
#include <cstdint>
#include <vector>

/**
 * \brief Find the k nearest neighbors for each point in query_points in target_points using a KDTree
 *
 * what it does:
 *  -build a kd_tree for the target_points
 *  -find for each point in query_points the k closest points in target_points
 *  -put the found indices and distances into the out_indices and out_dists matrices
 *
 * TODO:
 *  -implement this exercise
 *  -make a typedef for nanoflann::KDTreeEigenMatrixAdaptor
 *  -build the kd_tree
 *  -loop over the query_points and find the k closest points in the kd_tree with the function knnSearch
 *  -put the indices and distances in the respective output matrices
 *
 * \param[in] target_points The target points organized as a n x d matrix
 * \param[in] query_points The query points organized as a m x d matrix
 * \param[in] k The number of nearest neighbors
 * \param[out] out_indices The found indices in the target points; |query_points| rows with k entries each
 * \param[out] out_dists The distances to the found points in the target points; |query_points| rows with k entries each-
 */
using eidx = Eigen::DenseIndex;
using knn_kdt = nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixXd, 3, nanoflann::metric_L2>;

// void findKNN(const Eigen::MatrixXd &target_points,
//              const Eigen::MatrixXd &query_points,
//              const int k,
//              Eigen::MatrixXi &out_indices,
//              Eigen::MatrixXd &out_dists)
// {
//     knn_kdt kdtree(3, std::cref(target_points));
//     std::vector<size_t> indices(query_points.rows() * k);
//     std::vector<double> distances2(query_points.rows() * k);
//     nanoflann::KNNResultSet<double> res(k * query_points.rows());
//     res.init(indices.data(),distances2.data());
//     kdtree.index->findNeighbors(res, query_points.data(), nanoflann::SearchParams(10));
//     out_indices.resize(query_points.rows(), k);
//     out_dists.resize(query_points.rows(), k);
//     for(eidx q = 0; q < query_points.rows(); ++q)
//     {
//         for(eidx j = 0; j < k; ++j)
//         {
//             out_indices(q, j) = static_cast<eidx>(indices[static_cast<size_t>(q * k + j)]);
//             out_dists(q, j) = sqrt(distances2[q * k + j]);
//         }
//     }
// }

void findKNN(const Eigen::MatrixXd &target_points,
             const Eigen::MatrixXd &query_points,
             const int k,
             Eigen::MatrixXi &out_indices,
             Eigen::MatrixXd &out_dists)
{
    knn_kdt kdtree(3, target_points);
    kdtree.index->buildIndex();
    std::vector<eidx> indices(static_cast<size_t>(k), -1);
    std::vector<double> distances2(static_cast<size_t>(k), std::numeric_limits<double>::max());
    out_indices.resize(query_points.rows(), k);
    out_dists.resize(query_points.rows(), k);
    double qp[3];    
    for(eidx q = 0; q < query_points.rows(); ++q)
    {
        qp[0] = query_points.row(q)(0);
        qp[1] = query_points.row(q)(1);
        qp[2] = query_points.row(q)(2);
        kdtree.index->knnSearch(&qp[0], static_cast<size_t>(k), &indices[0], &distances2[0]);
        for(eidx j = 0; j < k; ++j)
        {
            out_indices(q, j) = indices[static_cast<size_t>(j)];
            out_dists(q, j) = sqrt(distances2[static_cast<size_t>(j)]);
        }
    }
}