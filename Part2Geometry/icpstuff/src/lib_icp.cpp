#include "lib_icp.h"
#include "lib_knn_query.h"
#include <iostream>
#include <vector>
#include <igl/opengl/glfw/Viewer.h>
/**
 * \brief Compute a rotation matrix around the z-axis which rotates points by degree
 *
 * what it does:
 *  -compute a rotation matrix rot which is (3 x 3)
 *
 * \param[in] degree The angle specifying how much we rotate around z-axis in degrees
 * \param[out] rot The rotation matrix (3 x 3)
 */
void rotZ(double degree, Eigen::Matrix3d &rot)
{
    rot = Eigen::AngleAxisd(degree * 2.0 * M_PI / 360.0 , Eigen::Vector3d::UnitZ());
    std::cout << "rot " << degree << " degrees around Z axis" << std::endl;
}

/**
 * \brief Apply a rotation and a transformation in the icp sense R * p + t
 *
 * what it does:
 *  -apply the rotation and translation to the points  R * p + t
 *
 * TODO:
 *  -implement this function
 *  -first rotate the points and then translate
 *
 * \param[in/out] points The points we want to transform
 * \param[in] rot The rotation matrix (3 x 3)
 * \param[in] trans The translation vector (1 x 3)
 */
void apply_icp_transformation(Eigen::MatrixXd &points,
                              const Eigen::MatrixXd &rot,
                              const Eigen::RowVector3d &trans)
{
    points *= rot.transpose();
    points.rowwise() += trans;
}

/**
 * \brief Apply a rotation and a transformation
 *
 * what it does:
 *  -apply the rotation and translation to the points in a local coordinate system: R * (p-p_mean) + t + p_mean
 *
 * TODO:
 *  -implement this function
 *  -center the points
 *  -rotate the points
 *  -move the points back
 *  -apply the translation
 *
 * \param[in/out] points The points we want to transform
 * \param[in] rot The rotation matrix (3 x 3)
 * \param[in] trans The translation vector (1 x 3)
 */
void apply_transformation(Eigen::MatrixXd &points,
                          const Eigen::MatrixXd &rot,
                          const Eigen::RowVector3d &trans)
{
    Eigen::RowVector3d c = points.colwise().mean();
    points.rowwise() -= c;
    points *= rot.transpose();
    points.rowwise() += c;
    points.rowwise() += trans;
}

/**
 * \brief Compute the optimal rotation and translation which minimizes the icp energy
 *
 * what it does:
 *  -compute the optimal rotation and translation which minimizes the icp energy:
 *      E(R,t) = ||q - R * p + t||²
 *  -using the SVD of the matrix H we can compute the optimal rot and trans which minimizes the squared distances
 *
 * TODO:
 *  -implement this function
 *  -sort the source and target according to the correspondences
 *  -compute the matrix H
 *  -compute the SVD of H by using the JacobiSVD iteration of eigen
 *  -compute the rot matrix
 *  -compute the trans matrix
 *
 * \param[in] source_points The points we want to transform eventually (n x 3)
 * \param[in] target_points The points we want to transform the source_points to (m x 3)
 * \param[in] corrs_source_target The corresponding indices (n x 2) 1. column (source) correspond to 2. column (target)
 * \param[in/out] rot The rotation matrix (3 x 3)
 * \param[in/out] trans The translation vector (1 x 3)
 * \param[in] debug Specify wether to print messages
 */
void optimal_transformation(const Eigen::MatrixXd &source_points,
                            const Eigen::MatrixXd &target_points,
                            const Eigen::MatrixXi &corrs_source_target,
                            Eigen::Matrix3d &rot,
                            Eigen::RowVector3d &trans,
                            const bool debug)
{
    Eigen::MatrixXd S(corrs_source_target.rows(), source_points.cols());
    Eigen::MatrixXd T(corrs_source_target.rows(), target_points.cols());
    for(Eigen::DenseIndex i = 0; i < corrs_source_target.rows(); ++i)
    {
        S.row(i) = Eigen::RowVector3d(source_points.row(corrs_source_target(i, 0)));
        T.row(i) = Eigen::RowVector3d(target_points.row(corrs_source_target(i, 1)));
    }
    Eigen::RowVector3d S_mean = S.colwise().mean();
    Eigen::RowVector3d T_mean = T.colwise().mean();
    Eigen::Matrix3d H = (S.rowwise() - S_mean).transpose() * (T.rowwise() - T_mean);
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d F = Eigen::MatrixXd::Identity(3, 3);
    F(2, 2) = (svd.matrixV() * svd.matrixU().transpose()).determinant();
    Eigen::Matrix3d R = svd.matrixV() * F * svd.matrixU().transpose();
    Eigen::RowVector3d t = T_mean - (S_mean * R.transpose());
    rot = R;
    trans = t;
}

Eigen::MatrixXi makeIndexPairs(const Eigen::MatrixXi& index, const Eigen::MatrixXd& distances, double epsilon)
{
    std::vector<Eigen::RowVector2i> pairs;
    for(Eigen::DenseIndex i = 0; i < index.rows(); ++i)
    {
        if(distances(i) <= epsilon)
        {
            pairs.push_back(Eigen::RowVector2i(i, index(i)));
        }
    }
    Eigen::MatrixXi indexpairs(static_cast<Eigen::DenseIndex>(pairs.size()), 2);
    for(size_t i = 0; i < pairs.size(); ++i)
    {
        indexpairs.row(static_cast<Eigen::DenseIndex>(i)) = pairs[i];
    }
    return indexpairs;
}


/**
 * \brief The icp algorithm
 *
 * what it does:
 *  -compute the optimal rotation and translation which minimizes the icp energy:
 *      E(R,t) = ||q - R * p + t||²
 *  -using the SVD of the matrix H we can compute the optimal rot and trans which minimizes the squared distances
 *
 * TODO:
 *  -implement this function
 *  -estimate the correspondences by using findKNN
 *  -calculate the optimal rotation and transformation by using optimalTransform
 *  -apply the transformation to the source points
 *  -iterate until max_iters is reached or the energy change is less than epsilon
 *
 * \param[in] source_points The points we want to transform eventually (n x 3)
 * \param[in] target_points The points we want to transform the source_points to (m x 3)
 * \param[in] max_iters The maximum number of iterations
 * \param[in] epsilon The accuracy we want to guarantee unless max_iters is reached
 * \param[in/out] rot The final rotation matrix (3 x 3)
 * \param[in/out] trans The final translation vector (1 x 3)
 * \param[in] debug Specify wether to print messages
 */
void icp(const Eigen::MatrixXd &source_points,
         const Eigen::MatrixXd &target_points,
         const int max_iters,
         const double epsilon,
         Eigen::Matrix3d &rot,
         Eigen::RowVector3d &trans,
         const double mincd,
         const bool showIntermediateSteps)
{
    Eigen::MatrixXd ws_source_points(source_points);
    int it = 0;
    double err = std::numeric_limits<double>::infinity();
    Eigen::MatrixXi knnindex;
    Eigen::MatrixXd knndistances;
    Eigen::Matrix3d tmp_rot;
    Eigen::RowVector3d tmp_trans;
    tmp_rot.setIdentity();
    tmp_trans.setZero();
    rot.setIdentity();
    trans.setZero();
    while(err > epsilon && it < max_iters)
    {        
        findKNN(target_points, ws_source_points, 1, knnindex, knndistances);
        err = knndistances.mean();
        std::cout << "Iteration: " << it << " MAE: " << err << "\n";
        if(knnindex.rows() == 0)
            break;
        auto ip = makeIndexPairs(knnindex, knndistances, mincd);

        if(showIntermediateSteps)
        {
            Eigen::MatrixXd S(ip.rows(), source_points.cols());
            Eigen::MatrixXd T(ip.rows(), target_points.cols());
            for(Eigen::DenseIndex i = 0; i < knnindex.rows(); ++i)
            {
                S.row(i) = Eigen::RowVector3d(ws_source_points.row(i));
                T.row(i) = Eigen::RowVector3d(target_points.row(knnindex(i, 0)));
            }
            igl::opengl::glfw::Viewer viewer;
            viewer.data().clear();
            viewer.data().point_size = 2.0;
            viewer.data().add_points(S, Eigen::RowVector3d(1.0, 0.0, 0.0));
            viewer.data().add_points(T, Eigen::RowVector3d(1.0, 1.0, 1.0));
            viewer.data().add_points(target_points, Eigen::RowVector3d(0.0, 0.0, 0.6));
            viewer.launch();
        }

        optimal_transformation(ws_source_points, target_points, ip, tmp_rot, tmp_trans);
        apply_icp_transformation(ws_source_points, tmp_rot, tmp_trans);
        trans = trans * tmp_rot.transpose() + tmp_trans;
        rot = tmp_rot * rot;        
        ++it;
    }
}