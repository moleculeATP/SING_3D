# include  <Eigen/Dense>

std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> 
        read_obj_cloud_points(const std::string& file_name);
void write_neighboring_graph(const std::string& file_name, const std::vector<Eigen::Vector3d> points, const Eigen::MatrixXd& adj_mat);