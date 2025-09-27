#include "distances.h"

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

std::pair<Eigen::MatrixXd, std::vector<std::vector<double>>> computeSINGDistances(
    const std::vector<Eigen::Vector3d>& points,
    const std::vector<Eigen::Vector3d>& normals,
    const std::string& filename,
    bool write,
    double density_weight,
    double normals_weight
){
    int n = points.size();
    Eigen::MatrixXd dist_mat = Eigen::MatrixXd::Zero(n, n);
    std::vector<std::vector<double>> lower_tri;
    std::vector<double> nn(n, 0.0);


    // Linear scan
    // for (int i = 0; i < n; i++) {
    //     double min_dist = std::numeric_limits<double>::max();
    //     for (int j = 0; j < n; j++) {
    //         if (i != j) {
    //             double dist = euclidean_distance(points[i], points[j]);
    //             if (dist < min_dist) {
    //                 min_dist = dist;
    //             }
    //         }
    //     }
    //     nn[i] = min_dist;
        
    nn = compute_nn_distances(points);

    for (int i = 0; i < n; i++) {
        std::vector<double> dists;
        for (int j = 0; j < i; j++) {
            double distance = euclidean_distance(points[i], points[j]) / (nn[i] + nn[j]) * 
                                std::pow(std::max(nn[i], nn[j]) / std::min(nn[i], nn[j]), density_weight) *
                                (1.0 + normals_weight * sqrt(2 * (1.0 - abs(normals[i].dot(normals[j])))));
            dist_mat(i, j) = distance;
            dist_mat(j, i) = distance;
            dists.push_back(distance);
        }
        
        lower_tri.push_back(dists);

        if(write) {
            std::cout << "Saving not implemented yet" << std::endl;
        }
    }

    return std::pair<Eigen::MatrixXd, std::vector<std::vector<double>>>{dist_mat, lower_tri};
}