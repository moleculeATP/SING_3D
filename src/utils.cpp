#include "utils.h"
#include <omp.h>
#include <numeric>
#include <unordered_map>


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

// std::vector<std::tuple<int, double, double>> compute_persistence_diagram(Edge_list edge_list, int num_pt, double threshold) {
//     // https://gudhi.inria.fr/doc/3.2.0/_rips_complex_2rips_distance_matrix_persistence_8cpp-example.html
    
//     int num_points = edge_list.size();
//     std::vector<int> indices(num_points);
//     std::iota(indices.begin(), indices.end(), 0); // Remplit de 0 à num_points-1

//     auto distance = [&](int i, int j) -> double {
//         return edge_list[i][j]; // ta distance arbitraire
//     };

//     // Crée le Rips complex
//     Gudhi::rips_complex::Rips_complex<double> rips_complex(indices, epsilon, distance);

//     // // threshold is 30 % of max distance in the matrix
//     // double max_dist = 0.0;
//     // for (const auto& row : edge_list) {
//     //     for (const auto& val : row) {
//     //         if (val > max_dist) {
//     //             max_dist = val;
//     //         }
//     //     }
//     // }
    
//     Simplex_tree simplex_tree;
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

std::vector<std::tuple<int, double, double>> compute_persistence_diagram(
    Edge_list edge_list, int num_pt, double threshold) {

    int num_points = num_pt;
    std::vector<int> indices(num_points);
    std::iota(indices.begin(), indices.end(), 0);

    auto pair_hash = [](const std::pair<int,int>& p) {
    return std::hash<long long>()(((long long)p.first << 32) ^ (long long)p.second);
    };

    std::cout << "Check here" << std::endl;

    std::unordered_map<std::pair<int,int>, double, decltype(pair_hash)> dist_map(0, pair_hash);

    for (const auto& e : edge_list) {
        int i = e.first.first;
        int j = e.first.second;
        double d = e.second;
        dist_map[{i, j}] = d;
        dist_map[{j, i}] = d; 
    }

    std::cout << "Check here" << std::endl;
    
    auto distance = [&](int i, int j) -> double {
    auto it = dist_map.find({i, j});
    if (it != dist_map.end())
        return it->second;
    return INF; 
    };

    std::cout << "Check here" << std::endl;


    Gudhi::rips_complex::Rips_complex<double> rips_complex(indices, threshold, distance);

    // std::vector<std::vector<double>> dist_matrix(num_points, std::vector<double>(num_points, INF));
    // for (const auto& e : edge_list) {
    //     dist_matrix[e.first.first][e.first.second] = e.second;
    //     dist_matrix[e.first.second][e.first.first] = e.second;
    // }

    // Rips_complex rips_complex(dist_matrix, threshold);

    
    std::cout << "Rips complex created." << std::endl;

    Simplex_tree simplex_tree;
    rips_complex.create_complex(simplex_tree, 1);
    std::cout << "Simplicial complex created with " << simplex_tree.num_simplices() << " simplices." << std::endl;

    Persistent_cohomology pcoh(simplex_tree);
    pcoh.init_coefficients(2);
    pcoh.compute_persistent_cohomology(0);
    std::cout << "Persistent cohomology computed." << std::endl;

    std::vector<std::tuple<int, double, double>> result;
    auto pairs = pcoh.get_persistent_pairs();
    for (const auto& interval : pairs) {
        int dim = simplex_tree.dimension(std::get<0>(interval));
        double birth = simplex_tree.filtration(std::get<0>(interval));
        double death = (std::get<1>(interval) == Simplex_tree::null_simplex()) 
                       ? std::numeric_limits<double>::infinity() 
                       : simplex_tree.filtration(std::get<1>(interval));
        result.emplace_back(dim, birth, death);
    }
    return result;
}



void print_barcode(const std::vector<std::tuple<int, double, double>>& diagram, int dim) {
    std::cout << "Persistence barcode:\n";
    for (const auto& [dim_, birth, death] : diagram) {
        if (dim_ != dim) continue; // Ignore negative dimensions
        std::cout << "dim " << dim << " : ["
                  << birth << ", "
                  << (death == std::numeric_limits<double>::infinity() ? "∞" : std::to_string(death))
                  << ")\n";
    }
}

void plot_barcode_terminal(const std::vector<std::tuple<int, double, double>>& diagram, double epsilon_max, int width) {
    std::cout << "Persistence barcode (dim 0) up to epsilon = " << epsilon_max << ":\n";

    for (const auto& [dim, birth, death] : diagram) {
        int start = static_cast<int>((birth / epsilon_max) * width);
        int end;
        if (death == std::numeric_limits<double>::infinity())
            end = width; // barre infinie
        else
            end = static_cast<int>((death / epsilon_max) * width);

        std::cout << "[";
        for (int i = 0; i < width; ++i) {
            if (i >= start && i < end)
                std::cout << "-";
            else
                std::cout << " ";
        }
        std::cout << "]\n";
    }
}

void export_persistence_csv(const std::vector<std::tuple<int, double, double>>& diagram,
                            const std::string& filename) {
    std::ofstream out(filename);
    out << "dim,birth,death\n";

    for (const auto& [dim, birth, death] : diagram) {
        out << dim << "," << birth << ",";
        if (death == std::numeric_limits<double>::infinity())
            out << "inf\n";
        else
            out << death << "\n";
    }
}