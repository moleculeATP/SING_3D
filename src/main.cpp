#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "io.h"


// computing SING distanceusing linear scan for the moment

double euclidean_distance(const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
    return (a - b).norm();
}

// to do : KD tree implementation
std::pair<Eigen::MatrixXd, std::vector<std::vector<double>>> computeSINGDistances(
    const std::vector<Eigen::VectorXd>& points,
    const std::string& filename = "",
    bool write = false,
    double density = 0.0
){
    int n = points.size();
    Eigen::MatrixXd dist_mat = Eigen::MatrixXd::Zero(n, n);
    std::vector<std::vector<double>> lower_tri;
    std::vector<double> nn(n, 0.0);

    // compute nearest neighbor for each point
    for (int i = 0; i < n; i++) {
        double min_dist = std::numeric_limits<double>::max();
        for (int j = 0; j < n; j++) {
            if (i != j) {
                double dist = euclidean_distance(points[i], points[j]);
                if (dist < min_dist) {
                    min_dist = dist;
                }
            }
        }
        nn[i] = min_dist;
    }

    for (int i = 0; i < n; i++) {
        std::vector<double> dists;
        for (int j = 0; j < i; j++) {
            double distance = euclidean_distance(points[i], points[j]) / (nn[i] + nn[j]) * 
                                std::pow(std::max(nn[i], nn[j]) / std::min(nn[i], nn[j]), density);
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


std::pair<std::vector<std::pair<int, int>>, Eigen::MatrixXd> extractSINGEdges(const Eigen::MatrixXd& dist_mat, double epsilon = 1.0) {
    int n = dist_mat.rows();
    std::vector<std::pair<int, int>> edges;
    Eigen::MatrixXd adj_mat = Eigen::MatrixXd::Zero(n, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            if (dist_mat(i, j) <= epsilon) {
                edges.push_back({i, j});
                adj_mat(i, j) = 1.0;
                adj_mat(j, i) = 1.0;
            }
        }
    }

    return {edges, adj_mat};
}

int main() {
    std::cout << "Sing 3D" << std::endl;

    std::string filename = "../datas/letter_A.obj";
    std::vector<Eigen::VectorXd> points = read_obj_cloud_points(filename);
    std::cout << "Number of points: " << points.size() << std::endl;

    auto [dist_mat, lower_tri] = computeSINGDistances(points, "", false, 5.0);
    std::cout << "Distance matrix computed." << std::endl;

    auto [edges, adj_mat] = extractSINGEdges(dist_mat, .8);
    std::cout << "Number of edges: " << edges.size() << std::endl;

    std::string out_filename = "../outputs/output_A.obj";
    write_neighboring_graph(out_filename, points, adj_mat);

    return 0;
}
