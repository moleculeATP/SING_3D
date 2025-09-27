#include <vector>
#include <iomanip>
#include <limits>
#include <tuple>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Rips_complex.h>
#include <gudhi/Persistent_cohomology.h>

using Distance_matrix = std::vector<std::vector<double>>;
using Simplex_tree = Gudhi::Simplex_tree<>;
using Rips_complex = Gudhi::rips_complex::Rips_complex<double>;
using Filtration_value = Simplex_tree::Filtration_value;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;

std::pair<std::vector<std::pair<int, int>>, Eigen::MatrixXd> extractSINGEdges(const Eigen::MatrixXd& dist_mat, double epsilon = 1.0);

std::vector<std::tuple<int, double, double>>
compute_persistence_diagram(const Distance_matrix& distance_matrix);

void print_barcode(const std::vector<std::tuple<int, double, double>>& diagram, int dim = 0);

void plot_barcode_terminal(const std::vector<std::tuple<int, double, double>>& diagram, double epsilon_max, int width = 50);

void export_persistence_csv(const std::vector<std::tuple<int, double, double>>& diagram,
                            const std::string& filename);