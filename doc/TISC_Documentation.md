# TISC Documentation

This document provides a comprehensive user manual and technical documentation for the TISC (Tectonics, Isostasy, Surface Processes, and Climate) model.

# Part 1: User Manual

## 1. Introduction

TISC is a numerical model that simulates the evolution of a landscape under the influence of various geological and climatic processes. It couples tectonic deformation, isostatic adjustment of the lithosphere, erosion and sedimentation by surface processes (fluvial, hillslope, and glacial), and the influence of climate. TISC is designed to be a flexible tool for researchers and students in geology, geophysics, and geomorphology.

## 2. Installation

### 2.1. Dependencies

To compile and run TISC, you will need:
*   A C compiler (e.g., `gcc` or `clang`).
*   A Fortran compiler (e.g., `gfortran`), only if the thin-sheet tectonic model is used.
*   `make` utility.

These tools are standard on most Linux and macOS systems.

### 2.2. Compilation

To compile the code, simply run the `make` command in the root directory of the project:

```sh
make
```

This will create the executables in the `bin/` directory.

### 2.3. Environment Setup

For the model to run correctly, you need to add the `bin/` and `script/` directories to your system's `PATH`. You also need to set the `tisc_dir` environment variable to the root directory of the TISC project.

The `make` command provides instructions on how to do this for `bash` and `csh` shells.

For **bash** (and compatible shells like `zsh`), add the following lines to your `~/.bashrc` or `~/.zshrc` file:

```sh
export tisc_dir=/path/to/your/tisc/folder
export PATH=$PATH:$tisc_dir/bin:$tisc_dir/script
```

Replace `/path/to/your/tisc/folder` with the actual absolute path to the TISC project directory.

For **csh** or **tcsh**, add these lines to your `~/.cshrc` file:

```csh
setenv tisc_dir /path/to/your/tisc/folder
set path = ($path $tisc_dir/bin $tisc_dir/script)
```

## 3. Running TISC

### 3.1. Basic Syntax

The basic syntax to run a TISC simulation is:

```sh
tisc [options] project_name
```

Where `project_name` is the base name for your project's input files (without extension). For example, if your parameters file is `my_project.PRM`, you would run `tisc my_project`.

### 3.2. Command-Line Options

TISC accepts several command-line options to control its behavior:

*   `-h`: Display help and information.
    *   `-hp`: Display the content of the default `template.PRM` file.
    *   `-hc`: Display the cleaned (without comments) `template.PRM` file.
    *   `-hu`: Display the `template.UNIT` file.
*   `-V[level]`: Set the verbose level (0-5). Higher levels print more information.
*   `-f[reformat]`: Reformat the `.PRM` file.
*   `-F[resume_filename]`: Resume a previous run from a `.all` file.
*   `-Q load_file_name`: Direct mode. Solves flexure for a given load file.
*   `-q parameter=value`: Override a parameter from the `.PRM` file.

### 3.3. Direct Mode (`-Q`)

This mode is for solving the flexural response to a single load without running a full time-dependent simulation. It's useful for quick tests of isostatic models.

### 3.4. Resuming a Run (`-F`)

You can resume a run that was previously stopped. TISC periodically saves its state in a `project_name.all` file. To resume, use the `-F` flag:

```sh
tisc -F project_name.all
```

## 4. Input Files

TISC uses a set of input files to define a simulation.

### 4.1. `*.PRM` (Parameters file)

This is the main configuration file for a project. It defines the model domain, grid, densities, time stepping, and parameters for the physical models. The file `doc/template.PRM` contains all available parameters with explanations.

### 4.2. `*.UNIT` (Unit files)

These files define tectonic or sedimentary events that occur at specific times. They are named `project_name1.UNIT`, `project_name2.UNIT`, etc. Each file represents a load (e.g., a thrust sheet, a volcanic extrusion, an erosional event) and contains parameters like the time of the event, density of the material, and a spatial distribution of the thickness.

### 4.3. Initial Condition Files

*   `*.ZINI`: Initial topography grid.
*   `*.WINI`: Initial deflection grid.
*   `*.RIV`: Defines initial river paths that are carved into the initial topography.

### 4.4. Time-Dependent Boundary Conditions

*   `*.SLV`: Sea-level variations over time.
*   `*.RAIN`: Spatially variable precipitation grid.
*   `*.EET`: Spatially variable effective elastic thickness (Te).

## 5. Output Files

TISC generates a variety of output files in the project directory.

*   `*.hrz`: Horizon data in a grid format. Each column represents a geological horizon.
*   `*.pfl`: Profile data along a cross-section defined in a `*.PRFL` file.
*   `*.xyw`: Drainage network data, including water discharge, sediment load, and topography for each cell.
*   `*.lak`: Information about lakes, including their extent, volume, and outlets.
*   `*.st`: Erosion and sedimentation rates for each cell.
*   `*.all`: Binary file containing the complete state of the model, used for resuming runs.
*   `*.ps`, `*.jpg`: Graphical outputs of the model state at different time steps.

---

# Part 2: Technical Documentation

## 1. Code Architecture

The TISC codebase is organized into several C source files, each responsible for a specific part of the model.

*   `tisc.c`: This is the main file containing the `main()` function. It orchestrates the overall simulation. It parses command-line arguments, calls the input routines, and manages the main time loop. Within the time loop, it calls functions for tectonic loading, isostatic flexure, and surface processes in sequence.

*   `tiscio.c`: This file handles most of the input and output operations. It contains functions for reading the main parameters file (`read_file_parameters`), unit files (`read_file_unit`), initial conditions, and for writing the various output files (`write_file_Blocks`, `write_file_drainage`, etc.).

*   `tisclib.c`: This library file contains many of the core numerical and physical subroutines. This includes memory allocation (`Allocate_Memory`), matrix definitions for the flexure solver (`defineLESalmostdiagonalmatrix`), and helper functions for managing the geological `Blocks`.

*   `surf_proc.c`: This is a large and central file containing the complete model for surface processes. It calculates water flow, erosion, transport, and sedimentation. It includes the logic for drainage network definition, lake formation and evolution, fluvial erosion/sedimentation, hillslope diffusion, and glacial processes.

## 2. Core Physical Models

### 2.1. Lithospheric Flexure

TISC models the isostatic response of the lithosphere to loads using a thin elastic or viscoelastic plate model.

#### 2.1.1. Governing Equation

The vertical deflection `w` of a thin plate under a vertical load `q` is described by the biharmonic equation:

∇²(D∇²w) + (ρ_m - ρ_i)gw = q

Where:
*   `D` is the flexural rigidity of the plate, given by `D = E * Te³ / (12 * (1 - ν²))`.
*   `E` is Young's modulus.
*   `Te` is the effective elastic thickness of the lithosphere (a key input parameter).
*   `ν` is Poisson's ratio.
*   `ρ_m` is the density of the underlying asthenosphere (`densasthen`).
*   `ρ_i` is the density of the infill material (`densinfill` or `densenv`).
*   `g` is the acceleration due to gravity.
*   `q` is the vertical load per unit area.

#### 2.1.2. Implementation

*   The `Elastic_Deflection` function in `tisc.c` is responsible for solving the flexure equation.
*   The equation is discretized using a finite-difference scheme. The coefficients of the resulting system of linear equations are assembled in the `defineLESalmostdiagonalmatrix` function (in `tisclib.c`).
*   The system is then solved using a custom solver for "almost diagonal" matrices (`solveLESalmostdiagonal`).
*   The primary variables in the code related to flexure are:
    *   `w`: A 2D array holding the deflection.
    *   `Dq`: A 2D array for the incremental load applied during a time step.
    *   `q`: A 2D array for the total accumulated load.
    *   `Te_default`: The default effective elastic thickness. `EET` is the 2D array for spatially variable Te.

#### 2.1.3. Viscoelastic Relaxation

If `isost_model = 2`, TISC simulates a viscoelastic plate. This is implemented in the `Viscous_Relaxation` function, which calculates the relaxation of the flexural stresses over time, causing the deflection to evolve from an elastic response towards a state of local isostasy. The relaxation time is controlled by the `tau` parameter.

---

### 2.2. Surface Processes: Water and Sediment Routing

The core of TISC's landscape evolution capabilities lies in `surf_proc.c`. This module simulates the movement of water and sediment across the surface, leading to the development of drainage networks and the carving and filling of basins. The entire process is managed by the `Surface_Transport` function.

#### 2.2.1. Drainage Network Definition (`Define_Drainage_Net`)

Before water and sediment can be routed, a drainage network must be established. TISC uses a priority-first search algorithm based on topography.

1.  **Sorting**: All grid cells are sorted by elevation in descending order. This is stored in the `sortcell` array.
2.  **Bottom-up Traversal**: The code iterates through the cells from the lowest to the highest. For each cell, it determines where it drains to.
3.  **Flow Direction**: For a given cell, the algorithm searches its 8 neighbors to find the path of steepest descent. This determines the flow direction.
4.  **Pits and Flats**: If a cell is lower than all its neighbors, it is a "pit" and becomes the seed of a new lake. If a cell has neighbors at the same elevation (a flat area), it is also treated as a potential lake area.
5.  **Node Types**: Each cell is classified into a `type` within the `drainage` structure:
    *   `'R'` (River): A standard cell that drains to a lower neighbor.
    *   `'L'` (Lake): A cell that is part of a lake body.
    *   `'E'` (Exit): A lake cell that is an outlet (spillway).

This process results in a complete drainage network, where every cell has a defined downstream neighbor.

#### 2.2.2. Lake Dynamics

Lakes are a critical feature of the TISC model. They form in topographic depressions and can be either exorheic (open, with an outlet) or endorheic (closed, internally drained).

##### Lake Identification

Lakes are identified during the `Define_Drainage_Net` process.
*   A pit (a cell lower than all its neighbors) initiates a new lake.
*   The lake expands by flooding adjacent cells until it finds one or more outlets (saddles in the topography).
*   All cells belonging to a lake are tagged with a lake ID number in `drainage[i][j].lake`. The properties of each lake are stored in the `Lake` structure array.

##### Endorheic vs. Exorheic Lakes

The distinction between an open and closed lake is determined by water balance, which is calculated in the `Calculate_Discharge` function.

*   An **exorheic (open) lake** has enough water input to fill up to its lowest outlet and spill over. Its water level (`Lake[il].alt`) is fixed by the elevation of its lowest outlet.
*   An **endorheic (closed) lake** does not have enough water to spill. Its water level is determined by the point where evaporation from the lake surface balances the total water inflow.

##### Determining Endorheic Lake Extent

This is a key process in TISC, handled in `Calculate_Discharge`.

1.  **Initial Assumption**: All lakes are initially treated as full up to their spillways.
2.  **Water Balance Calculation**: For each lake `il`, the code calculates the total water input `Lake_Input_Discharge(il)` (from river inflow and direct precipitation) and the total potential evaporation `lake_evap` (sum of `evaporation` over all cells in the lake).
3.  **Shrinking the Lake**: If `Lake_Input_Discharge(il) < lake_evap`, the lake is potentially endorheic and cannot maintain its level. The model then enters a loop to shrink the lake.
4.  **`Attempt_Delete_Node_From_Lake`**: The code iteratively removes the highest-elevation node from the lake's edge.
5.  **Splitting Lakes (`Divide_Lake`)**: If removing a node would disconnect the lake into two or more separate water bodies, the `Divide_Lake` function is called to handle the split.
6.  **Re-calculating Balance**: After a node is removed (and the lake area shrinks), the water balance is re-evaluated. This process continues until `Lake_Input_Discharge(il) >= lake_evap`, at which point the lake has reached a stable, smaller extent where inflow equals evaporation. The final water level is the elevation of the highest remaining lake cell.

#### 2.2.3. Water Balance (`Calculate_Discharge`)

The water balance for each grid cell is calculated by routing water downstream, accounting for various sources and sinks. The governing principle is conservation of mass:

`Q_out = Q_in + P - E - Q_underground`

*   `Q_in`: Discharge flowing into the cell from upstream neighbors.
*   `P` (Precipitation): Calculated in `Calculate_Precipitation_Evaporation`. It can be spatially uniform, dependent on elevation (`Krain`), or based on a sophisticated orographic precipitation model (`hydro_model = 2` or `3`).
*   `E` (Evaporation): A constant rate `evaporation_ct` is applied over lakes.
*   `Q_underground` (Seepage): A simple Darcy-like law is used to simulate water loss to the subsurface, controlled by the `permeability` parameter.
*   `Q_out`: The final discharge leaving the cell (`drainage[i][j].discharge`).

#### 2.2.4. Sediment Balance and Transport (`Fluvial_Transport`)

The transport of sediment is coupled with the water flow. The change in sediment mass in a cell is:

`dM_sed/dt = Eroded_mass + Sed_in - Sed_out`

##### Erosion Models (`erosed_model`)

TISC supports several models for fluvial erosion, selected by the `erosed_model` parameter. A key concept is the **transport capacity** (`transp_capacity_eq`), which is the maximum sediment load a river can carry, typically modeled as a function of discharge `Q` and slope `S`.

*   **Undercapacity Model (e.g., `erosed_model = 2`)**: Erosion or deposition is proportional to the difference between the current sediment load (`drainage[i][j].masstr`) and the transport capacity.
    *   If `transp_capacity_eq > masstr`, erosion occurs.
    *   If `transp_capacity_eq < masstr`, deposition occurs.
    *   Key variables: `K_river_cap` (for transport capacity), `erodibility` (length scale for erosion), `l_fluv_sedim` (length scale for deposition).

*   **Stream Power Law (e.g., `erosed_model = 6`)**: Erosion rate is a power-law function of discharge and slope, often representing shear stress. This model also accounts for the inhibiting effect of sediment load.
    *   Erosion `~ Q^m * S^n * (1 - masstr/transp_capacity_eq)`
    *   Key variables: `erodibility` (coefficient for the stream power law).

The actual erosion is performed by the `Erode()` function, which removes material from the uppermost geological `Blocks`.

##### Sedimentation

Deposition occurs when the sediment load exceeds the transport capacity.
*   **In Rivers**: The excess sediment is deposited, raising the bed elevation.
*   **In Lakes (`Lake_Fill`)**: When a river enters a lake, its velocity drops, and it deposits its sediment load. The `Lake_Fill` function distributes the incoming sediment (`drainage[row][col].masstr`) within the lake basin, starting from the river mouth and propagating inwards. For endorheic lakes, if any sediment remains after the distribution process, it is deposited uniformly over the entire lake area.

The deposition is handled by the `Sediment()` function, which adds thickness to the topmost sediment `Block`.
