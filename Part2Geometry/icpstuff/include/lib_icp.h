#ifndef ATCG1_LIB_ICP_H
#define ATCG1_LIB_ICP_H

#include <Eigen/Dense>

void rotZ(double degree, Eigen::Matrix3d &rot);

void apply_icp_transformation(Eigen::MatrixXd &points,
                              const Eigen::MatrixXd &rot,
                              const Eigen::RowVector3d &trans);

void apply_transformation(Eigen::MatrixXd &points,
                          const Eigen::MatrixXd &rot,
                          const Eigen::RowVector3d &trans);

void optimal_transformation(const Eigen::MatrixXd &source_points,
                            const Eigen::MatrixXd &target_points,
                            const Eigen::MatrixXi &corrs_source_target,
                            Eigen::Matrix3d &rot,
                            Eigen::RowVector3d &trans,
                            const bool debug = false);

void icp(const Eigen::MatrixXd &source_points, const Eigen::MatrixXd &target_points, const int max_iters,
         const double epsilon, Eigen::Matrix3d &rot, Eigen::RowVector3d &trans, const double, const bool);

#endif //ATCG1_LIB_ICP_H
