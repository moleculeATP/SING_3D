#include "distances.h"
#include "chrono"
#include <omp.h>

double euclidean_distance(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    return (a - b).norm();
}

double ellipsoid_distance(
    const Eigen::Vector3d& a,
    const Eigen::Vector3d& b,
    const Eigen::Matrix3d& S  // SPD matrix
) {
    Eigen::Vector3d diff = a - b;
    return std::sqrt(diff.transpose() * S * diff);
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

// std::vector<Eigen::Matrix3d> compute_S_matrices(
//     const std::vector<Eigen::Vector3d>& points,
//     double direction_weight
// ){
//      size_t n = points.size();
//     std::vector<Eigen::Matrix3d> S_list(n);

//     std::vector<CGAL::Point_3<CGAL::Simple_cartesian<double>>> cgal_points;
//     cgal_points.reserve(n);
//     for (const auto& p : points) cgal_points.emplace_back(p(0), p(1), p(2));

//     CGAL::Kd_tree<CGAL::Search_traits_3<CGAL::Simple_cartesian<double>>> tree(cgal_points.begin(), cgal_points.end());

//     const int K = 5; // A CHANGER PLUS TARD

//     for (size_t i = 0; i < n; i++) {
//         CGAL::Orthogonal_k_neighbor_search<CGAL::Search_traits_3<CGAL::Simple_cartesian<double>>> search(tree, cgal_points[i], K+1);
//         std::vector<Eigen::Vector3d> neighbors;
//         neighbors.reserve(K+1);

//         int count = 0;
//         for (auto it = search.begin(); it != search.end() && count < K+1; ++it) {
//             if (it->first != cgal_points[i]) { 
//                 neighbors.emplace_back(it->first.x(), it->first.y(), it->first.z());
//                 count++;
//             }
//         }
//         neighbors.push_back(points[i]); 

//         Eigen::Vector3d mean = Eigen::Vector3d::Zero();
//         for (const auto& p : neighbors) mean += p;
//         mean /= neighbors.size();

//         Eigen::Matrix3d C = Eigen::Matrix3d::Zero();
//         for (const auto& p : neighbors) {
//             Eigen::Vector3d diff = p - mean;
//             C += diff * diff.transpose();
//         }
//         C /= neighbors.size();

//         Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(C);
//         Eigen::Vector3d eigenvalues = es.eigenvalues();
//         Eigen::Matrix3d eigenvectors = es.eigenvectors();

//         Eigen::Vector3d inv_axes;
//         double lambda_max = eigenvalues.maxCoeff();
//         for (int d = 0; d < 3; d++) {
//             inv_axes(d) = 1.0 / (std::pow(eigenvalues(d) / lambda_max, direction_weight) + 1e-8);
//         }

//         // S = V * D * V^T
//         Eigen::Matrix3d D = inv_axes.asDiagonal();
//         S_list[i] = eigenvectors * D * eigenvectors.transpose();
//     }

//     return S_list;   
// }

std::vector<Eigen::Matrix3d> compute_S_matrices(
    const std::vector<Eigen::Vector3d>& points,
    double direction_weight
){
    size_t n = points.size();
    std::vector<Eigen::Matrix3d> S_list(n);

    std::vector<CGAL::Point_3<CGAL::Simple_cartesian<double>>> cgal_points;
    cgal_points.reserve(n);
    for (const auto& p : points)
        cgal_points.emplace_back(p(0), p(1), p(2));

    CGAL::Kd_tree<CGAL::Search_traits_3<CGAL::Simple_cartesian<double>>> tree(cgal_points.begin(), cgal_points.end());

    const int K = 20;
    const double epsilon = 10;

    for (size_t i = 0; i < n; i++) {
        CGAL::Orthogonal_k_neighbor_search<CGAL::Search_traits_3<CGAL::Simple_cartesian<double>>> search(
            tree, cgal_points[i], K + 1);

        std::vector<Eigen::Vector3d> neighbors;
        std::vector<double> weights;
        neighbors.reserve(K + 1);
        weights.reserve(K + 1);

        int count = 0;
        for (auto it = search.begin(); it != search.end() && count < K + 1; ++it) {
            if (it->first == cgal_points[i]) continue;
            Eigen::Vector3d q(it->first.x(), it->first.y(), it->first.z());
            double dist = (q - points[i]).norm();
            if (dist <= epsilon) {
                double w = std::exp(-dist * dist / (2 * epsilon * epsilon));
                neighbors.push_back(q);
                weights.push_back(w);
            }
            count++;
        }
        neighbors.push_back(points[i]);
        weights.push_back(1.0);

        double weight_sum = 0.0;
        Eigen::Vector3d mean = Eigen::Vector3d::Zero();
        for (size_t j = 0; j < neighbors.size(); j++) {
            mean += weights[j] * neighbors[j];
            weight_sum += weights[j];
        }
        mean /= weight_sum;

        Eigen::Matrix3d C = Eigen::Matrix3d::Zero();
        for (size_t j = 0; j < neighbors.size(); j++) {
            Eigen::Vector3d diff = neighbors[j] - mean;
            C += weights[j] * diff * diff.transpose();
        }
        C /= weight_sum;

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(C);
        Eigen::Vector3d eigenvalues = es.eigenvalues();
        Eigen::Matrix3d eigenvectors = es.eigenvectors();

        Eigen::Vector3d inv_axes;
        double lambda_max = eigenvalues.maxCoeff();
        for (int d = 0; d < 3; d++) {
            inv_axes(d) = 1.0 / (std::pow(eigenvalues(d) / lambda_max, direction_weight) + 1e-8);
        }

        Eigen::Matrix3d D = inv_axes.asDiagonal();
        S_list[i] = eigenvectors * D * eigenvectors.transpose();
    }

    return S_list;
}



// std::vector<Eigen::Matrix3d> compute_S_matrices(
//     const std::vector<Eigen::Vector3d>& points,
//     const std::vector<Eigen::Vector3d>& normals,
//     double direction_weight
// ){
//     size_t n = points.size();
//     std::vector<Eigen::Matrix3d> S_list(n);

//     std::vector<CGAL::Point_3<CGAL::Simple_cartesian<double>>> cgal_points;
//     cgal_points.reserve(n);
//     for (const auto& p : points) cgal_points.emplace_back(p(0), p(1), p(2));

//     CGAL::Kd_tree<CGAL::Search_traits_3<CGAL::Simple_cartesian<double>>> tree(cgal_points.begin(), cgal_points.end());

//     const int K = 5;

//     for (size_t i = 0; i < n; i++) {
//         CGAL::Orthogonal_k_neighbor_search<CGAL::Search_traits_3<CGAL::Simple_cartesian<double>>> search(tree, cgal_points[i], K+1);
//         std::vector<Eigen::Vector3d> neighbors;
//         neighbors.reserve(K+1);

//         int count = 0;
//         for (auto it = search.begin(); it != search.end() && count < K+1; ++it) {
//             if (it->first != cgal_points[i]) {
//                 neighbors.emplace_back(it->first.x(), it->first.y(), it->first.z());
//                 count++;
//             }
//         }
//         neighbors.push_back(points[i]);

//         Eigen::Vector3d mean = Eigen::Vector3d::Zero();
//         for (const auto& p : neighbors) mean += p;
//         mean /= neighbors.size();

//         Eigen::Vector3d normal = normals[i].normalized();
//         Eigen::Matrix3d P = Eigen::Matrix3d::Identity() - normal * normal.transpose();

//         std::vector<Eigen::Vector3d> projected;
//         projected.reserve(neighbors.size());
//         for (const auto& p : neighbors) projected.push_back(P * (p - mean));

//         Eigen::Matrix3d C = Eigen::Matrix3d::Zero();
//         for (const auto& v : projected) C += v * v.transpose();
//         C /= projected.size();

//         Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(C);
//         Eigen::Vector3d eigenvalues = es.eigenvalues();
//         Eigen::Matrix3d eigenvectors = es.eigenvectors();

//         Eigen::Matrix3d V;
//         V.col(2) = normal;

//         Eigen::Vector3d axis1 = eigenvectors.col(2) - normal.dot(eigenvectors.col(2)) * normal;
//         Eigen::Vector3d axis2 = eigenvectors.col(1) - normal.dot(eigenvectors.col(1)) * normal;
//         axis1.normalize();
//         axis2.normalize();
//         V.col(0) = axis1;
//         V.col(1) = axis2;

//         double lambda_max = eigenvalues.maxCoeff();
//         Eigen::Vector3d inv_axes;
//         inv_axes(0) = 1.0 / (std::pow(eigenvalues(2) / lambda_max, direction_weight) + 1e-8);
//         inv_axes(1) = 1.0 / (std::pow(eigenvalues(1) / lambda_max, direction_weight) + 1e-8);
//         inv_axes(2) = 1.0 / (std::pow(eigenvalues(0) / lambda_max, direction_weight) + 1e-8);

//         Eigen::Matrix3d D = inv_axes.asDiagonal();
//         S_list[i] = V * D * V.transpose();
//     }

//     return S_list;
// }
std::vector<Eigen::Matrix3d> compute_S_matrices(
    const std::vector<Eigen::Vector3d>& points,
    const std::vector<Eigen::Vector3d>& normals,
    double direction_weight
){
    size_t n = points.size();
    std::vector<Eigen::Matrix3d> S_list(n);

    int K = 20;
    double epsilon = 100;

    std::vector<CGAL::Point_3<CGAL::Simple_cartesian<double>>> cgal_points;
    cgal_points.reserve(n);
    for (const auto& p : points) cgal_points.emplace_back(p(0), p(1), p(2));

    CGAL::Kd_tree<CGAL::Search_traits_3<CGAL::Simple_cartesian<double>>> tree(cgal_points.begin(), cgal_points.end());

    for (size_t i = 0; i < n; i++) {
        CGAL::Orthogonal_k_neighbor_search<CGAL::Search_traits_3<CGAL::Simple_cartesian<double>>> search(tree, cgal_points[i], K+1);
        
        std::vector<Eigen::Vector3d> neighbors;
        neighbors.reserve(K+1);

        int count = 0;
        for (auto it = search.begin(); it != search.end() && count < K+1; ++it) {
            Eigen::Vector3d neighbor(it->first.x(), it->first.y(), it->first.z());
            if ((neighbor - points[i]).norm() <= epsilon) {
                neighbors.push_back(neighbor);
            }
            count++;
        }
        neighbors.push_back(points[i]);

        Eigen::Vector3d normal = normals[i].normalized();
        Eigen::Matrix3d P = Eigen::Matrix3d::Identity() - normal * normal.transpose();

        double weight_sum = 0.0;
        Eigen::Vector3d mean = Eigen::Vector3d::Zero();
        std::vector<double> weights;
        weights.reserve(neighbors.size());

        for (const auto& p : neighbors) {
            Eigen::Vector3d diff = p - points[i];
            double w = std::exp(-diff.squaredNorm() / (2 * epsilon * epsilon));
            weights.push_back(w);
            mean += w * p;
            weight_sum += w;
        }
        mean /= weight_sum;

        Eigen::Matrix3d C = Eigen::Matrix3d::Zero();
        for (size_t j = 0; j < neighbors.size(); j++) {
            Eigen::Vector3d diff = neighbors[j] - mean;
            diff = P * diff;
            C += weights[j] * diff * diff.transpose();
        }
        C /= weight_sum;

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(C);
        Eigen::Vector3d eigenvalues = es.eigenvalues();
        Eigen::Matrix3d eigenvectors = es.eigenvectors();

        Eigen::Matrix3d V;
        V.col(2) = normal;

        Eigen::Vector3d axis1 = eigenvectors.col(1) - normal.dot(eigenvectors.col(1)) * normal;
        axis1.normalize();
        V.col(0) = axis1;
        V.col(1) = normal;

        double lambda_max = eigenvalues.maxCoeff();
        Eigen::Vector3d inv_axes;
        inv_axes(0) = 1.0 / (std::pow(eigenvalues(2) / lambda_max, direction_weight) + 1e-8);
        inv_axes(1) = 1.0 / (std::pow(eigenvalues(1) / lambda_max, direction_weight) + 1e-8);
        inv_axes(2) = inv_axes(1);

        Eigen::Matrix3d D = inv_axes.asDiagonal();
        S_list[i] = V * D * V.transpose();
    }

    return S_list;
}





std::pair<Eigen::SparseMatrix<double>, Edge_list> computeSINGDistances(
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
    if(normals_weight > 0.){
        if (normals.size() != points.size()) {
            throw std::invalid_argument("Points and normals must have the same size.");
        }
        for (int i = 0; i < n; i++) {
            normals_normalized[i] = normals[i].normalized();
        }
    }

    Eigen::SparseMatrix<double> mat(n, n);

    std::vector<Eigen::Triplet<double>> triplets;
    int alloc_size = std::min(300 * n, 10000000);
    triplets.reserve(alloc_size);

    Edge_list edges;
    edges.reserve(alloc_size);

    std::vector<double> nn = compute_nn_distances(points);

    // Thread-safe containers
    std::vector<std::vector<Eigen::Triplet<double>>> triplets_private;
    std::vector<Edge_list> edges_private;

    std::cout << "dqzdqd" << std::endl;
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
                edges_private[tid].emplace_back(std::make_pair(i, j), distance);
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

std::pair<Eigen::SparseMatrix<double>, Edge_list> computeAnisotropeDistances(
    const std::vector<Eigen::Vector3d>& points,
    const std::vector<Eigen::Vector3d>& normals,
    double direction_weight,
    double treshold
){
    int n = points.size();
    Eigen::SparseMatrix<double> mat(n, n);
    Edge_list edges;

    std::vector<Eigen::Triplet<double>> triplets;
    int alloc_size = std::min(300 * n, 10000000);
    triplets.reserve(alloc_size);

    edges.reserve(alloc_size);

    std::vector<Eigen::Matrix3d> S_matrices = compute_S_matrices(points, direction_weight);
    // std::vector<Eigen::Matrix3d> S_matrices = compute_S_matrices(points, normals, direction_weight);
    write_Ellipse_field("Ellipse_field.obj", points, S_matrices);

    for(int i = 0; i < n; i++){
        for(int j = 0; j < i; j++){
            double distance = std::max(ellipsoid_distance(points[i], points[j], S_matrices[i]),
                                  ellipsoid_distance(points[j], points[i], S_matrices[j]));

            if(distance <= treshold){
                triplets.emplace_back(i, j, distance);
                triplets.emplace_back(j, i, distance);
                edges.emplace_back(std::make_pair(i, j), distance);
            }            
        }
    }
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return {mat, edges}; 
}

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

    for (int t = 0; t < num_threads; t++) {
        edges.insert(edges.end(), edges_private[t].begin(), edges_private[t].end());
        triplets.insert(triplets.end(), triplets_private[t].begin(), triplets_private[t].end());
    }

    adj_mat.setFromTriplets(triplets.begin(), triplets.end());
    return {edges, adj_mat};
}