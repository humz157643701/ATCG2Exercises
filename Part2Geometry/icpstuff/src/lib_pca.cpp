#include <iostream>
#include <Eigen/Eigenvalues>
#include "lib_pca.h"

/**
 * \brief Centers the given data X
 *
 * what it does:
 *  -compute the mean of X where mean is (1 x d)
 *  -substract the from X
 *
 * TODO:
 *  -implement this function
 *
 * \param[in/out] X The data points organized as a n x d matrix
 */
void center_data(Eigen::MatrixXd &X)
{
    X.rowwise() -= X.colwise().mean();
}


/**
 * \brief Performs the principal component analysis on the data matrix X
 *
 * TODO:
 *  -implement this function
 *  -construct the covarianve matrix
 *  -compute the eigendecomposition of the covariance matrix
 *  -scale the eigenvectors to the length of the standard-deviations
 *
 * \param[in] X The data points organized as a n x d matrix
 * \param[out] stddevs The standard deviations of X as a 1 x d vector
 * \param[out] dirs The principal directions of the data X as a d x d matrix
 */
void pca(Eigen::MatrixXd &X, Eigen::VectorXd& stddevs, Eigen::MatrixXd& dirs)
{
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
    auto pcs = svd.matrixV();
    auto svals = svd.singularValues();
    for(int i = 0; i < svals.rows(); ++i)
    {
        svals(i) = sqrt((svals(i) * svals(i)) / (static_cast<double>(X.rows()) - 1.0));
    }
    stddevs = svals;
    dirs = pcs;
}