#ifndef ATCG1_LIB_KNN_QUERY_H
#define ATCG1_LIB_KNN_QUERY_H

#include <Eigen/Dense>

void findKNN(const Eigen::MatrixXd &target_points,
             const Eigen::MatrixXd &query_points,
             const int k,
             Eigen::MatrixXi &indices,
             Eigen::MatrixXd &distances);

#endif //ATCG1_LIB_KNN_QUERY_H
