# SING_3D  
An implementation of the **SING algorithm in 3D**.


<summary><strong style="font-size: 1.5em;">&nbsp;Linux / macOS (conda environment)</strong></summary>

---

### 0. Prerequisites

Install the required system tools:

```bash
sudo apt install cmake g++
```

Then create and activate a conda environment with the required libraries:

```bash
conda create -n sing3d_env -c conda-forge gudhi eigen boost cgal
conda activate sing3d_env
```

Make sure `$CONDA_PREFIX` is correctly set once the environment is activated (it will be used by CMake to find the packages).

---

### 1. Installing

- Clone the repository:

```bash
git clone https://github.com/your_username/SING_3D.git
cd SING_3D
```

- Compile:

```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j 12
```

---

### 2. Test the installation

Once compiled, the executable will be available in the `build/` directory:

```bash
./sing_3D
```

If it runs without errors, the installation was successful.

---

### 3. Example

After building, you can test a minimal example (replace with your own data):

```bash
./sing_3D path/to/input_file.txt <epsilon> <p> <input_filename> <out_graph> <out_csv>
```

It will compute de SING graph and persistence diagram.
You can use the python code provided to generate plot of persistence. 
You can also use the script 
```bash
python run_program_linux.py
```
to automaticly compile, run the c++ code on the specified arguments and display the diagrams.

</details>