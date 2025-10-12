import subprocess

epsilon = 10
p = 3
q = 0
dim_to_display = 0
FILE_NAME = "teaser"
filename = "../datas/" + FILE_NAME + ".obj"
out_graph = "../outputs/graphs/" + FILE_NAME + ".obj"
out_csv = "../outputs/persistence/" + FILE_NAME + ".csv"

dim = str(dim_to_display)
subprocess.run(["make"], cwd="build")
subprocess.run(["./sing_3D", str(epsilon), str(p), str(q), filename, out_graph, out_csv], cwd="build")
subprocess.run(["python", "visu.py", out_csv, dim], cwd="python")
