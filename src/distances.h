#include <CGAL/Simple_cartesian.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "io.h"
#include <chrono>


double euclidean_distance(const Eigen::Vector3d& a, const Eigen::Vector3d& b);

std::vector<double> compute_nn_distances(const std::vector<Eigen::Vector3d>& points) ;

std::pair<Eigen::MatrixXd, std::vector<std::vector<double>>> computeSINGDistances(
    const std::vector<Eigen::Vector3d>& points,
    const std::vector<Eigen::Vector3d>& normals,
    const std::string& filename = "",
    bool write = false,
    double density_weight = 0.0,
    double normals_weight = 0.0
);



