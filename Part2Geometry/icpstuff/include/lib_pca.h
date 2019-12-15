#ifndef ATCG1_LIB_PCA_H
#define ATCG1_LIB_PCA_H

#include <Eigen/Core>

void center_data(Eigen::MatrixXd &X);

void pca(Eigen::MatrixXd &X, Eigen::VectorXd& stddevs, Eigen::MatrixXd& dirs);


#endif //ATCG1_LIB_PCA_H
