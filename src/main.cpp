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
    if (argc < 9) {
        std::cerr << "Usage: " << argv[0] << " <epsilon> <mode> <p> <input_filename> <out_graph> <out_csv>\n";
        return 1;
    }

    double epsilon;
    float p, q, w;
    std::string filename, out_graph, out_csv;
    std::string mode = argv[2];
    double threshold;

 
    epsilon = std::atof(argv[1]);
    p = std::atof(argv[3]);
    q = std::atof(argv[4]);
    w = std::atof(argv[5]);
    filename = argv[6];
    out_graph = argv[7];
    out_csv = argv[8];
    threshold = std::atof(argv[9]);
    


    auto start = std::chrono::high_resolution_clock::now();

    std::cout << "Sing 3D" << std::endl;

    auto [points, normals] = read_obj_cloud_points(filename);
    std::cout << "Number of points: " << points.size() << std::endl;
     
    Eigen::SparseMatrix<double> dist_mat;
    Edge_list full_edges_list;
    if (mode == "sing") {
        std::cout << "Using SING distances with p=" << p << std::endl;
        auto res = computeSINGDistances(points, normals, "", false, p, q, threshold);
        dist_mat = res.first;
        full_edges_list =  res.second;
    } else if(mode == "anisotropic") {
        std::cout << "Using Anisotropic distances with w=" << w << std::endl;
        auto res = computeAnisotropeDistances(points, normals, w, threshold);
        dist_mat = res.first;
        full_edges_list = res.second;
    }
    
    
    std::cout << "Distance matrix computed." << std::endl;

    int num_pt = points.size();
    auto [edges, adj_mat] = extractSINGEdges(dist_mat, epsilon);
    std::cout << "Number of edges: " << edges.size() << std::endl;

    write_neighboring_graph(out_graph, points, adj_mat);

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
