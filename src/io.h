#ifndef IO_H
#define IO_H

# include  <Eigen/Dense>
# include  <Eigen/Sparse>
# include "utils.h"

std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> 
        read_obj_cloud_points(const std::string& file_name);
void write_neighboring_graph(const std::string& file_name, const std::vector<Eigen::Vector3d> points, Adjacency_matrix adj_mat);

#endif // IO_H