# TISC

**Tectonics, Isostasy, Surface Processes, and Climate Planform Modeling**

TISC is a C program for Linux and macOS that simulates the evolution of a landscape under the influence of various geological and climatic processes. It couples tectonic deformation, isostatic adjustment of the lithosphere, erosion and sedimentation by surface processes (fluvial, hillslope, and glacial), and the influence of climate. 

> *Note: The pseudo-2D (vertical cross-section) version of TISC is called **tAo** and is distributed separately.*

* **Main Author:** Daniel Garcia-Castellanos (CSIC)
* **Details:** See authorship, disclaimers, and citation info in `doc/TISC_Documentation.md`.
* **License:** See details in the `doc/` directory.

---

## Installation Prerequisites

### Native Compilation (macOS / Linux)
To compile and run TISC, you will need:
*   A C compiler (e.g., `gcc` or `clang`).
*   A Fortran compiler (e.g., `gfortran`), only if the thin-sheet tectonic model is used.
*   `make` utility.
*   (Optional) GMT 4 or Python 3 for graphical output.

---

## Building TISC

1. Uncompress the downloaded file into a directory named `tisc/`.
2. Compile by typing `make` in the main directory. This will create the executable `tisc` in the `bin/` directory.
3. Add `tisc/bin/` and `tisc/script/` to your system `PATH` and define the `tisc_dir` variable in your environment:
   
   ```sh
   export tisc_dir=/path/to/your/tisc
   export PATH=$PATH:$tisc_dir/bin:$tisc_dir/script
   ```

## Getting Started

* Try the examples in the `tisc/demo/` directory to ensure that the code is working properly and that the graphic output (`*.ps` or `*.jpg`/`*.svg`) is generated correctly.
* Read the `doc/TISC_Documentation.md` file to learn how to configure and use TISC.