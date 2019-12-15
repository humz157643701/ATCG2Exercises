#include <cstdlib>
#include <iostream>

//set this to supress libigl viewer help
#define IGL_VIEWER_VIEWER_QUIET

#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <Eigen/Core>
#include <algorithm>

#include "readXYZ.h"
#include "readASCIIMATRIX.h"

#include "lib_pca.h"
#include "lib_knn_query.h"
#include "lib_icp.h"

/**
 * \brief Scales the data X into the unit hypercube
 *
 * sometimes data has a funny scale and it is outside of the libigl-viewer-frustrum
 * this is a convenience function to scale the data into the unit hypercube
 *
 * \params[in/out] X The data points organized as a n x d matrix
 */
void scale_to_unit_cube(Eigen::MatrixXd &X)
{
    double s = 1.0 / X.colwise().maxCoeff().maxCoeff();
    X = (X.array() * s).matrix();
}

void scale_to_unit_cube(Eigen::MatrixXd& A, Eigen::MatrixXd& B)
{
    double s = 1.0 / std::max(A.colwise().maxCoeff().maxCoeff(), B.colwise().maxCoeff().maxCoeff());
    A *= s;
    B *= s;   
}


/**
 * \brief First exercise
 *
 * what it does:
 *  -load filename.xyz point cloud into V
 *  -V is (n x 3)
 *
 * TODO:
 *  -implement this exercise
 *  -complete the center_data and pca function
 *  -calculate pca(V, stddevs, dirs)
 *  -scale dirs by stddevs
 *  -plot V and the scaled dirs using the libigl viewer
 *
 * \param[in] filename The path to the .xyz ascii point cloud file as string
 */
void exercise1(const std::string &filename)
{
    Eigen::MatrixXd V;
    bool load_success = igl::readXYZ(filename, V);
    if (!load_success)
    {
        std::cerr << "could not load file: " << filename << std::endl;
        return;
    }
    scale_to_unit_cube(V);
     // subtract mean
    center_data(V);
    // a) draw figurine    
    igl::opengl::glfw::Viewer viewer;
    viewer.data().point_size = 2.0;
    viewer.data().add_points(V, Eigen::RowVector3d(1.0, 1.0, 1.0));    

    // b)   
    Eigen::MatrixXd pcs;
    Eigen::VectorXd stddevs;
    pca(V, stddevs, pcs);
    Eigen::Matrix3d colors;
    colors <<   1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0;
    viewer.data().add_edges(Eigen::Matrix3d().Zero(), pcs.transpose(), colors);
    viewer.launch();
}

/**
 * \brief Second exercise
 *
 * what it does:
 *  -load filename.txt point cloud into V
 *  -V is (n x 100)
 *
 * TODO:
 *  -implement this exercise
 *  -complete the center_data and pca function
 *  -calculate pca(V, stddevs, dirs)
 *  -project V onto dirs
 *  -decide based on the eigenvalues which axes are most relevant to plot
 *  -plot most relevant axes using the libigl viewer
 *
 * \param[in] filename The path to the .txt ascii multi-dimensional point cloud file as string
 */
void exercise2(const std::string &filename)
{
    Eigen::MatrixXd V;

    bool load_success = igl::readASCIIMATRIX(filename, V);
    if (!load_success)
    {
        std::cerr << "could not load file: " << filename << std::endl;
        return;
    }

    // normalize data
    scale_to_unit_cube(V);
    // center data
    center_data(V);
    // pca
    Eigen::MatrixXd pcs;
    Eigen::VectorXd stddevs;
    pca(V, stddevs, pcs);
    // project data matrix onto principal components
    Eigen::MatrixXd projected = V * pcs;

    // slice projected data, take 3 components (those seem to be the most important.
    // all following components share roughtly the same, low-ish variance.)
    // Also, 3 dimensions can be displayed nicely using libigl

    igl::opengl::glfw::Viewer viewer;
    viewer.data().point_size = 2.0;
    viewer.data().add_points(projected.leftCols(3), Eigen::RowVector3d(1.0, 1.0, 1.0));
    viewer.launch();

    // no idea what that is
}


/**
 * \brief Third exercise
 *
 * what it does:
 *  -load filename_source.off point cloud into V_source
 *  -load filename_target.obj point cloud into V_target
 *  -V_source is (n x 3), V_target is (m x 3)
 *
 * TODO:
 *  -implement this exercise
 *  -calculate a rotation matrix and a translation vector
 *  -apply the transformation to V_source by using the function apply_transformation
 *  -complete the icp function and apply icp to V_source
 *  -plot the starting V_source as a e.g. red point cloud
 *  -plot the target V_target as e e.g. green point cloud
 *  -plot the final transformed V_source as a mesh
 *
 * \param[in] filename_source The path to the .off mesh file as string
 * \param[in] filename_target The path to the .obj mesh file as string
 */
void exercise3(const std::string &filename_source, const std::string &filename_target){
    Eigen::MatrixXd V_source, V_target, V_start;
    Eigen::MatrixXi F_source, F_target;

    bool load_success_source = igl::readOFF(filename_source, V_source, F_source);

    if (!load_success_source)
    {
        std::cerr << "could not load file: " << filename_source << std::endl;
        return;
    }

    bool load_success_target = igl::readOBJ(filename_target, V_target, F_target);
    if (!load_success_target)
    {
        std::cerr << "could not load file: " << filename_target << std::endl;
        return;
    }

    // Normalize point clouds. Global scale is irrelevant for this task.
    scale_to_unit_cube(V_source, V_target);

    igl::opengl::glfw::Viewer viewer1;
    Eigen::RowVector3d tt(0.3, 0.0, 0.0);
    Eigen::Matrix3d rt = Eigen::AngleAxisd(1.5, Eigen::Vector3d(1.0, 1.0, 1.0).normalized()).toRotationMatrix();
    apply_transformation(V_source, rt, tt);
    viewer1.data().point_size = 2.0;
    viewer1.data().add_points(V_source, Eigen::RowVector3d(1.0, 0.0, 0.0));
    viewer1.data().add_points(V_target, Eigen::RowVector3d(0.0, 0.0, 1.0));
    viewer1.launch();
    // apply icp
    Eigen::Matrix3d rot;
    Eigen::RowVector3d trans;
    icp(V_source, V_target, 500, 0.01, rot, trans, 1.0, false);
    apply_icp_transformation(V_source, rot, trans);
    // show result
     igl::opengl::glfw::Viewer viewer2;
    viewer2.data().clear();
    viewer2.data().point_size = 2.0;
    viewer2.data().add_points(V_source, Eigen::RowVector3d(1.0, 0.0, 0.0));
    viewer2.data().add_points(V_target, Eigen::RowVector3d(0.0, 0.0, 1.0));
    viewer2.launch();

    
}


/**
* \brief The main function called when running this program
*
* what it does:
*  -check provided filenames
*  -run both exercises in a row
*
*  \param[in] argc The number of arguments to the binary
*  \param[in] argv The array of arguments to the binary
*/
int main(int argc, char *argv[])
{
    std::string filename1, filename2, filename3, filename4;

    if (argc == 5)
    {
        filename1 = argv[1];
        filename2 = argv[2];
        filename3 = argv[3];
        filename4 = argv[4];
    }
    else
    {
        std::cerr << "please call assignmentsheet1 like this: " << std::endl;
        std::cerr << "./assignmentsheet1 data/figurine.xyz data/measurements.txt data/maxear.off data/maxsimple.obj" << std::endl;
        return EXIT_FAILURE;
    }

    try
    {
        exercise1(filename1);
        exercise2(filename2);
        exercise3(filename3, filename4);
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    }
    
   

	return EXIT_SUCCESS;
}
