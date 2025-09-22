#include <CGAL/Simple_cartesian.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <chrono>
#include "io.h"
#include "distances.h"
#include "utils.h"



int main() {
    std::cout << "Sing 3D" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    float epsilon = 1.2;
    float p = 3.0;

    std::string filename = "../datas/letter_A.obj";
    std::vector<Eigen::Vector3d> points = read_obj_cloud_points(filename);
    std::cout << "Number of points: " << points.size() << std::endl;

    auto [dist_mat, lower_tri] = computeSINGDistances(points, "", false, p);
    std::cout << "Distance matrix computed." << std::endl;

    auto [edges, adj_mat] = extractSINGEdges(dist_mat, epsilon);
    std::cout << "Number of edges: " << edges.size() << std::endl;

    std::string out_filename = "../outputs/output_A.obj";
    write_neighboring_graph(out_filename, points, adj_mat);


    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Time spend : " << elapsed.count() << " s\n";

    return 0;
}
