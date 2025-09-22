#include "utils.h"

std::pair<std::vector<std::pair<int, int>>, Eigen::MatrixXd> extractSINGEdges(const Eigen::MatrixXd& dist_mat, double epsilon) {
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