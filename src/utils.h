#include <vector>
#include <Eigen/Dense>

std::pair<std::vector<std::pair<int, int>>, Eigen::MatrixXd> extractSINGEdges(const Eigen::MatrixXd& dist_mat, double epsilon = 1.0);