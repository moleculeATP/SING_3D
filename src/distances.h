#ifndef DISTANCES_H
#define DISTANCES_H

#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "io.h"
#include <chrono>


using Edge_list = std::vector<std::pair<std::pair<int,int>,double>>;
double euclidean_distance(const Eigen::Vector3d& a, const Eigen::Vector3d& b);

std::vector<double> compute_nn_distances(const std::vector<Eigen::Vector3d>& points) ;

std::pair<Eigen::SparseMatrix<double>, Edge_list> computeSINGDistances(
    const std::vector<Eigen::Vector3d>& points,
    const std::vector<Eigen::Vector3d>& normals,
    const std::string& filename,
    bool write,
    double density_weight = 0.0,
    double normals_weight = 0.0,
    double treshold = 2.
);

#endif // DISTANCES_H



