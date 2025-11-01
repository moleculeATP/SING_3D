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
double ellipsoid_distance(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Matrix3d& S);

std::vector<double> compute_nn_distances(const std::vector<Eigen::Vector3d>& points);
std::vector<Eigen::Matrix3d> compute_S_matrices(const std::vector<Eigen::Vector3d>& points, const std::vector<Eigen::Vector3d>& normals, double direction_weight);
std::vector<Eigen::Matrix3d> compute_S_matrices(const std::vector<Eigen::Vector3d>& points, double direction_weight);

std::pair<Eigen::SparseMatrix<double>, Edge_list> computeSINGDistances(
    const std::vector<Eigen::Vector3d>& points,
    const std::vector<Eigen::Vector3d>& normals,
    const std::string& filename,
    bool write,
    double density_weight = 0.0,
    double normals_weight = 0.0,
    double treshold = 2.
);

std::pair<Eigen::SparseMatrix<double>, Edge_list> computeAnisotropeDistances(
    const std::vector<Eigen::Vector3d>& points,
    const std::vector<Eigen::Vector3d>& normals,
    double direction_weight = 0.0,
    double treshold = 2.0
);

std::pair<std::vector<std::pair<int, int>>, Adjacency_matrix> extractSINGEdges(Distance_matrix dist_mat, double epsilon = 1.0);

#endif // DISTANCES_H



