// just some custom io tools. When this will be merge with RsR, we will use their io tools.

#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <ostream>

#include "io.h"

std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> 
        read_obj_cloud_points(const std::string& file_name){

    std::ifstream file(file_name);
    std::vector<Eigen::Vector3d> points = {};
    std::vector<Eigen::Vector3d> normals = {};
    if (!file.is_open()) {
        std::cerr << "Error: could not open file " << file_name << std::endl;
        return std::make_pair(points, normals);
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
        } else if (line.substr(0, 3) == "vn ") {
            std::istringstream iss(line.substr(2));
            std::vector<double> coords;
            double coord;
            while (iss >> coord) {
                coords.push_back(coord);
            }
            Eigen::Vector3d normal;
            for (size_t i = 0; i < 3; ++i) {
                normal(i) = coords[i];
            }
            normals.push_back(normal);
        }
    }

    file.close();
    return std::make_pair(points, normals);
}

void write_neighboring_graph(const std::string& file_name, const std::vector<Eigen::Vector3d> points, Adjacency_matrix adj_mat){ 
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
            if (adj_mat.coeff(i, j) > 0) {
                file << "l " << i + 1 << " " << j + 1 << "\n"; // obj format is 1-indexed
            }
        }
    }
    file.close();
}

void write_Ellipse_field(const std::string& file_name, const std::vector<Eigen::Vector3d> points, const std::vector<Eigen::Matrix3d> S_matrices){
    // write in a .obj format whit verices as v and axes of the ellipses as l. A new file is created or overwritten if it already exists.
    std::cout << "Writing ellipse field to " << file_name << std::endl;
    std::ofstream file(file_name);
    int n = points.size();

    for(int i = 0; i < n; i++){
        file << "v " << points[i](0) << " " << points[i](1) << " " << points[i](2) << "\n";

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(S_matrices[i]);
        Eigen::Vector3d eigenvalues = es.eigenvalues();
        Eigen::Matrix3d eigenvectors = es.eigenvectors();

        for(int d = 0; d < 3; d++){
            Eigen::Vector3d axis = eigenvectors.col(d) * std::sqrt(1.0 / eigenvalues(d));
            axis *= 1; // SCALE FOR VISUALIZATION
            file << "v " << points[i](0) + axis(0) << " " << points[i](1) + axis(1) << " " << points[i](2) + axis(2) << "\n";
            file << "l " << i * 4 + 1 << " " << i * 4 + 2 + d << "\n";
        }
    }
    file.close();
}