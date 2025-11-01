import subprocess

epsilon = 20
mode = "anisotropic" # anisotropic or sing
p = 0
q = 0
w = 1.
thresold = 200
dim_to_display = 0
FILE_NAME = "light_tree_cleaned"
filename = "../datas/" + FILE_NAME + ".obj"
out_graph = "../outputs/graphs/" + FILE_NAME + ".obj"
out_csv = "../outputs/persistence/" + FILE_NAME + ".csv"

dim = str(dim_to_display)
subprocess.run(["make"], cwd="build")
subprocess.run(["./sing_3D", str(epsilon), mode, str(p), str(q), str(w), filename, out_graph, out_csv, str(thresold)], cwd="build")
# subprocess.run(["python", "visu.py", out_csv, dim], cwd="python")
