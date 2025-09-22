# include  <Eigen/Dense>

std::vector<Eigen::VectorXd> read_obj_cloud_points(const std::string& file_name);
void write_neighboring_graph(const std::string& file_name, const std::vector<Eigen::VectorXd> points, const Eigen::MatrixXd& adj_mat);