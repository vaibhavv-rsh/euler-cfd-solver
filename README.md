# euler-cfd-solver

A small C++ solver for the compressible Euler equations (shock-tube problems).

This repository contains simple 1D and 2D Euler solvers with a parallel 2D variant that uses OpenMP. It uses `yaml-cpp` (bundled under `external/`) for parameter loading and `CMake` as the build system.

**Features**

- 1D shock-tube solver (`1d_shock_tube`).
- 2D shock-tube solver (`2d_shock_tube`).
- OpenMP-parallel 2D solver (`parallel_2d_shock_tube`).

**Prerequisites**

- A C++17 compiler (GCC/Clang) with OpenMP support for the parallel build.
- `cmake` (minimum required in the project is `2.8.12`).
- No external `yaml-cpp` install is required — the project builds the copy in `external/yaml-cpp`.

**Build**
From the repository root:

```bash
mkdir -p build && cd build
cmake ..
make -j$(nproc)
```

To build a Release binary:

```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

**Available executables**

- `1d_shock_tube` — 1D Euler solver (built from `src/euler_1d.cpp`).
- `2d_shock_tube` — 2D Euler solver (built from `src/euler_2d.cpp`).
- `parallel_2d_shock_tube` — 2D solver with OpenMP parallelization (built from `src/euler_2d_parallel.cpp`).

After building, the binaries are in the `build/` directory. Example runs below assume you are inside `build/`.

**Run examples**
Use the YAML parameter files in `config/`.

```bash
# 1D example
./1d_shock_tube

# 2D example
./2d_shock_tube

# Parallel 2D (uses OpenMP)
# run with 4 threads (set `OMP_NUM_THREADS` as needed)
OMP_NUM_THREADS=4 ./parallel_2d_shock_tube
```

**Output**
Simulation output files are written into `solutions_2d/` (CSV and VTK files are included in the repository as samples). Filenames follow the pattern `solution_<nx>x<ny>_<step>.{csv,vtk}`.

**Configuration**

- `config/params.yaml` — parameter file used by the 1D solver.
- `config/params_2d.yaml` — parameter file used by the 2D solvers.

**Development notes**

- The project uses `external/yaml-cpp` and adds it with `add_subdirectory`. If you update or replace that folder, ensure `CMakeLists.txt` is updated accordingly.
- The parallel target links to `OpenMP::OpenMP_CXX`; make sure your compiler provides OpenMP support when building `parallel_2d_shock_tube`.

**License**
This project is provided under the terms in the `LICENSE` file.

**Contributing / Contact**
Feel free to open issues or pull requests. For quick questions, refer to the source files in `src/` and headers in `include/`.

# euler-cfd-solver
