#include "distances.h"
#include "chrono"
#include <omp.h>

double euclidean_distance(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    return (a - b).norm();
}

std::vector<double> compute_nn_distances(const std::vector<Eigen::Vector3d>& points) {
    std::vector<CGAL::Point_3<CGAL::Simple_cartesian<double>>> cgal_points;
    cgal_points.reserve(points.size());
    for (const auto& p : points)
        cgal_points.emplace_back(p(0), p(1), p(2));

    CGAL::Kd_tree<CGAL::Search_traits_3<CGAL::Simple_cartesian<double>>> tree(cgal_points.begin(), cgal_points.end());

    std::vector<double> nn_distances;
    for (const auto& p : cgal_points) {
        CGAL::Orthogonal_k_neighbor_search<CGAL::Search_traits_3<CGAL::Simple_cartesian<double>>> search(tree, p, 2);
        auto it = search.begin();
        ++it;
        nn_distances.push_back(std::sqrt(it->second));
    }

    return nn_distances;
}

std::pair<Eigen::SparseMatrix<double>, std::vector<Edge>> computeSINGDistances(
    const std::vector<Eigen::Vector3d>& points,
    const std::vector<Eigen::Vector3d>& normals,
    const std::string& filename,
    bool write,
    double density_weight,
    double normals_weight,
    double treshold
){
    auto start = std::chrono::high_resolution_clock::now();
    int n = points.size();

    std::vector<Eigen::Vector3d> normals_normalized(n);
    for (int i = 0; i < n; i++) {
        normals_normalized[i] = normals[i].normalized();
    }

    Eigen::SparseMatrix<double> mat(n, n);

    std::vector<Eigen::Triplet<double>> triplets;
    int alloc_size = std::min(300 * n, 10000000);
    triplets.reserve(alloc_size);

    std::vector<Edge> edges;
    edges.reserve(alloc_size);

    std::vector<double> nn = compute_nn_distances(points);

    // Thread-safe containers
    std::vector<std::vector<Eigen::Triplet<double>>> triplets_private;
    std::vector<std::vector<Edge>> edges_private;

    int num_threads = 1;
    #ifdef _OPENMP
        num_threads = omp_get_max_threads();
    #endif
    triplets_private.resize(num_threads);
    edges_private.resize(num_threads);

    std::cout << "Using " << num_threads << " threads." << std::endl;

    #pragma omp parallel for if(num_threads > 1) schedule(dynamic)
    for (int i = 0; i < n; i++) {
        int tid = 0;
        #ifdef _OPENMP
                tid = omp_get_thread_num();
        #endif
        for (int j = 0; j < i; j++) {
            double distance = euclidean_distance(points[i], points[j]) /
                              (nn[i] + nn[j]) *
                              std::pow(std::max(nn[i], nn[j]) / std::min(nn[i], nn[j]), density_weight);

            if (normals_weight > 0.0) {
                distance *= (1.0 + normals_weight * std::sqrt(2 * (1.0 - std::abs(normals_normalized[i].dot(normals_normalized[j])))));
            }

            if (distance <= treshold){
                triplets_private[tid].emplace_back(i, j, distance);
                triplets_private[tid].emplace_back(j, i, distance);
                edges_private[tid].emplace_back(i, j, distance);
            }
        }

        if (write && i % 10000 == 0) {
            #pragma omp critical
            std::cout << "Saving not implemented yet" << std::endl;
        }

        if (i % 10000 == 0) {
            #pragma omp critical
            std::cout << "Processed " << i << " / " << n << " points." << std::endl;
        }
    }

    // Merge thread-local results
    for (int t = 0; t < num_threads; t++) {
        triplets.insert(triplets.end(), triplets_private[t].begin(), triplets_private[t].end());
        edges.insert(edges.end(), edges_private[t].begin(), edges_private[t].end());
    }

    mat.setFromTriplets(triplets.begin(), triplets.end());

    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end - start).count();
    std::cout << "Time " << elapsed << " s" << std::endl;

    return {mat, edges};
}