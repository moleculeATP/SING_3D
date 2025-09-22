// just some custom io tools. When this will be merge with RsR, we will use their io tools.

#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <ostream>

#include "io.h"

std::vector<Eigen::Vector3d> read_obj_cloud_points(const std::string& file_name){
    std::ifstream file(file_name);
    std::vector<Eigen::Vector3d> points = {};
    if (!file.is_open()) {
        std::cerr << "Error: could not open file " << file_name << std::endl;
        return points;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.substr(0, 2) == "v ") {
            std::istringstream iss(line.substr(2));
            std::vector<double> coords;
            double coord;
            while (iss >> coord) {
                coords.push_back(coord);
            }
            Eigen::Vector3d point;
            for (size_t i = 0; i < 3; ++i) {
                point(i) = coords[i];
            }
            points.push_back(point);
        }
    }

    file.close();
    return points;
}

void write_neighboring_graph(const std::string& file_name, const std::vector<Eigen::Vector3d> points, const Eigen::MatrixXd& adj_mat) {
    // write in a .obj format whit verices as v and edges as l. A new file is created or overwritten if it already exists.
    std::cout << "Writing neighboring graph to " << file_name << std::endl;
    std::ofstream file(file_name);
    
    for (const auto& point : points) {
        file << "v";
        for (int i = 0; i < point.size(); ++i) {
            file << " " << point(i);
        }
        file << "\n";
    }
    int n = adj_mat.rows();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            if (adj_mat(i, j) > 0) {
                file << "l " << i + 1 << " " << j + 1 << "\n"; // obj format is 1-indexed
            }
        }
    }
    file.close();
}