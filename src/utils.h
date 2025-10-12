#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <iomanip>
#include <limits>
#include <tuple>
#include <utility>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Rips_complex.h>
#include <gudhi/Sparse_rips_complex.h>
#include <gudhi/Persistent_cohomology.h>

using Edge_list = std::vector<std::pair<std::pair<int,int>,double>>;
using Distance_matrix = Eigen::SparseMatrix<double>;
using Adjacency_matrix = Eigen::SparseMatrix<bool>;
using Simplex_tree = Gudhi::Simplex_tree<>;
using Rips_complex = Gudhi::rips_complex::Rips_complex<double>;
using Filtration_value = Simplex_tree::Filtration_value;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;
constexpr double INF = std::numeric_limits<double>::infinity();

std::pair<std::vector<std::pair<int, int>>, Adjacency_matrix> extractSINGEdges(Distance_matrix dist_mat, double epsilon = 1.0);

std::vector<std::tuple<int, double, double>>
compute_persistence_diagram(Edge_list edge_list, int num_pt, double threshold = 2.);

void print_barcode(const std::vector<std::tuple<int, double, double>>& diagram, int dim = 0);

void plot_barcode_terminal(const std::vector<std::tuple<int, double, double>>& diagram, double epsilon_max, int width = 50);

void export_persistence_csv(const std::vector<std::tuple<int, double, double>>& diagram,
                            const std::string& filename);

#endif // UTILS_H