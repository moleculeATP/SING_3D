#include <CGAL/Simple_cartesian.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include "io.h"
#include "distances.h"
#include "utils.h"



int main(int argc, char* argv[]) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <epsilon> <p> <input_filename> <out_graph> <out_csv>\n";
        return 1;
    }

    double epsilon = std::atof(argv[1]);
    float p = std::atof(argv[2]);
    float q = std::atof(argv[3]);
    std::string filename = argv[4];
    std::string out_graph = argv[5];
    std::string out_csv = argv[6];
    double threshold = 5.;

    auto start = std::chrono::high_resolution_clock::now();

    std::cout << "Sing 3D" << std::endl;
    std::cout << "epsilon = " << epsilon << ", p = " << p << ", filename = " << filename << std::endl;

    auto [points, normals] = read_obj_cloud_points(filename);
    std::cout << "Number of points: " << points.size() << std::endl;

    auto [dist_mat, full_edges_list] = computeSINGDistances(points, normals, "", false, p, q, threshold);
    std::cout << "Distance matrix computed." << std::endl;

    int num_pt = points.size();
    auto [edges, adj_mat] = extractSINGEdges(dist_mat, epsilon);
    std::cout << "Number of edges: " << edges.size() << std::endl;

    // write_neighboring_graph(out_graph, points, adj_mat);

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Total time spend : " << elapsed.count() << " s\n";

    auto diag = compute_persistence_diagram(full_edges_list, num_pt, threshold);
    // print_barcode(diag, 1);
    // plot_barcode_terminal(diag, 1.0, 100);

    std::cout << "Diagram size: " << diag.size() << std::endl;
    export_persistence_csv(diag, out_csv);

    return 0;
}
