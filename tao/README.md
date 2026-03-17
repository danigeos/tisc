# tAo

**1D (Pseudo-2D, Vertical Cross-Section) Modelling of Flexural Isostasy, Erosion, Orographic Precipitation, and Foreland Basin Formation Using Finite Differences**

tAo is a C program for Linux and macOS designed to calculate the flexure of the lithosphere and, more generally, to numerically simulate the formation of foreland basins in pseudo-2D (vertical cross section). 

> *Note: The pseudo-3D version of tAo is called **TISC** and is distributed separately.*

* **Main Author:** Daniel Garcia-Castellanos (CSIC)
* **Details:** See authorship, disclaimers, and citation info in `doc/tAo_Documentation.md`.
* **License:** See details in the `doc/` directory.

---

## Installation Prerequisites

### Option A: Native Compilation (macOS / Linux)
tAo has been developed in C and is ready for compilation in a macOS or Linux environment. You will also need **GMT 4** or **Python** for graphic output.
```sh
sudo apt-get update
sudo apt-get install gmt ksh csh tofrodos
```

### Option B: Virtual Machine (Ubuntu + GMT4)
1. Install VirtualBox.
2. Download the pre-configured VM.
3. Import the VM with Ubuntu 20.04 + GMT 4.5.18 (`VirtualBox > File > Import...`).

---

## Building tAo

1. Uncompress the downloaded file into a directory named `tao/`.
2. Check the compiling options in `tao/config.mk`.
3. Compile by typing `make` in the main directory. This will create the executable `tao` in the `bin/` directory.
4. Add `tao/bin/` and `tao/script/` to your system `PATH` and define the `tao_dir` variable in your environment:
   ```sh
   export tao_dir=/path/to/tao
   export PATH=$PATH:$tao_dir/bin:$tao_dir/script
   ```

## Getting Started

* Try the examples in the `tao/demo/` directory to ensure that the code is working properly and that the graphic output (`*.ps` or `*.png`) is generated correctly.
* Read the `doc/tAo_Documentation.md` file to learn how to configure and use tAo.
