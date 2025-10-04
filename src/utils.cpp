#include "utils.h"
#include <omp.h>


std::pair<std::vector<std::pair<int, int>>, Adjacency_matrix> extractSINGEdges(Distance_matrix dist_mat, double epsilon) {
    int n = dist_mat.rows();
    std::vector<std::pair<int, int>> edges;
    std::vector<Eigen::Triplet<bool>> triplets;

    Adjacency_matrix adj_mat(n, n);   

    // Containers thread-local
    int num_threads = 1;
    #ifdef _OPENMP
        num_threads = omp_get_max_threads();
    #endif

    std::vector<std::vector<std::pair<int,int>>> edges_private(num_threads);
    std::vector<std::vector<Eigen::Triplet<bool>>> triplets_private(num_threads);

    #pragma omp parallel for if(num_threads > 1) schedule(dynamic)
    for (int i = 0; i < n; i++) {
        int tid = 0;
        #ifdef _OPENMP
                tid = omp_get_thread_num();
        #endif

        for (int j = 0; j < i; j++) {
            double val = dist_mat.coeff(i, j);
            if (val > 0 && val <= epsilon) {
                edges_private[tid].emplace_back(i, j);
                triplets_private[tid].emplace_back(i, j, true);
                triplets_private[tid].emplace_back(j, i, true);
            }
        }
    }

    // Fusion des résultats
    for (int t = 0; t < num_threads; t++) {
        edges.insert(edges.end(), edges_private[t].begin(), edges_private[t].end());
        triplets.insert(triplets.end(), triplets_private[t].begin(), triplets_private[t].end());
    }

    adj_mat.setFromTriplets(triplets.begin(), triplets.end());
    return {edges, adj_mat};
}

// std::vector<std::tuple<int, double, double>> compute_persistence_diagram(const Distance_matrix& distance_matrix) {
//     // https://gudhi.inria.fr/doc/3.2.0/_rips_complex_2rips_distance_matrix_persistence_8cpp-example.html
    
//     Rips_complex rips_complex(distance_matrix, 2);
//     Simplex_tree simplex_tree;

//     // threshold is 30 % of max distance in the matrix
//     double max_dist = 0.0;
//     for (const auto& row : distance_matrix) {
//         for (const auto& val : row) {
//             if (val > max_dist) {
//                 max_dist = val;
//             }
//         }
//     }

//     double threshold = 2;
//     rips_complex.create_complex(simplex_tree, 3); // au moins 1 pour dim 0

//     Persistent_cohomology pcoh(simplex_tree);
//     pcoh.init_coefficients(2);
//     pcoh.compute_persistent_cohomology();

//     std::vector<std::tuple<int, double, double>> result;
//     auto pairs = pcoh.get_persistent_pairs();
//     for (const auto& interval : pairs) {
//         int dim = simplex_tree.dimension(std::get<0>(interval));
//         double birth = simplex_tree.filtration(std::get<0>(interval));
//         double death = (std::get<1>(interval) == Simplex_tree::null_simplex()) 
//                        ? std::numeric_limits<double>::infinity() 
//                        : simplex_tree.filtration(std::get<1>(interval));
//         result.emplace_back(dim, birth, death);
//     }
//     return result;
// }

// void print_barcode(const std::vector<std::tuple<int, double, double>>& diagram, int dim) {
//     std::cout << "Persistence barcode:\n";
//     for (const auto& [dim_, birth, death] : diagram) {
//         if (dim_ != dim) continue; // Ignore negative dimensions
//         std::cout << "dim " << dim << " : ["
//                   << birth << ", "
//                   << (death == std::numeric_limits<double>::infinity() ? "∞" : std::to_string(death))
//                   << ")\n";
//     }
// }

// void plot_barcode_terminal(const std::vector<std::tuple<int, double, double>>& diagram, double epsilon_max, int width) {
//     std::cout << "Persistence barcode (dim 0) up to epsilon = " << epsilon_max << ":\n";

//     for (const auto& [dim, birth, death] : diagram) {
//         int start = static_cast<int>((birth / epsilon_max) * width);
//         int end;
//         if (death == std::numeric_limits<double>::infinity())
//             end = width; // barre infinie
//         else
//             end = static_cast<int>((death / epsilon_max) * width);

//         std::cout << "[";
//         for (int i = 0; i < width; ++i) {
//             if (i >= start && i < end)
//                 std::cout << "-";
//             else
//                 std::cout << " ";
//         }
//         std::cout << "]\n";
//     }
// }

// void export_persistence_csv(const std::vector<std::tuple<int, double, double>>& diagram,
//                             const std::string& filename) {
//     std::ofstream out(filename);
//     out << "dim,birth,death\n";

//     for (const auto& [dim, birth, death] : diagram) {
//         out << dim << "," << birth << ",";
//         if (death == std::numeric_limits<double>::infinity())
//             out << "inf\n";
//         else
//             out << death << "\n";
//     }
// }