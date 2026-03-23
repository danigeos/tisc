# TISC Documentation

This document provides a comprehensive user manual and technical documentation for the TISC (Tectonics, Isostasy, Surface Processes, and Climate) model.

## Table of Contents
* Part 1: User Manual
  * 1. Introduction
  * 2. Installation
  * 3. Running TISC
  * 4. Input Files
  * 5. Output Files
* Part 2: Technical Documentation
  * 1. Code Architecture
  * 2. Core Physical Models
  * 3. Computational Performance and Scaling

---

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

TISC and tAo share a unified build system and a common codebase (`tao+tisc_commons`). 
To compile the code, simply run the `make` command in the root directory of the project. 
You can configure compilers and optimization flags in the `config.mk` file:

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

```text
tisc  project  -B<bound_type> -D[xmin/xmax/ymin/ymax] -d<dx>[/<dy>]
      -e<Kriv>[/<Kdif>] -f[2] -F[<file>] -h[i|u|p|c] -l -M<lith_type>
      -N<Nx>[/<Ny>] -o -P[0|1|2|3][geom] -p<rain>[/<Krain>][/<evap>] -Q<file>
      -q<param>=<value> -R<random_topo>[/<seed>] -r<e|c|i|m|a><dens>
      -s<solver> -T<eet> -t<i|f|d|e|v|r><time> -V[<level>] 
      -v<num>/<vx>/<vy>
```

**Options:**
*   `project` is the root name for all files (i.e., if modelname='test' then the parameters file will be 'test.PRM' and the first load/fault 'test1.UNIT'). 
*   `-B` Sets boundary conditions (see doc/template.PRM file).
*   `-D` Set the model region (domain) [m].
*   `-d` To set dx and dy grid cell dimensions.
*   `-e` To set the constants of fluvial transport capacity [kg/m3] and diffusive transport [m2/a]. `-e0` suppresses surface transport calculations.
*   `-F` To continue (resume) a previous run by reading `<file>`. Default file_name is 'project.all'. If specified, `<file>` must be a binary file originally written by TISC in 'project.all' after completing a run. Put this option always *after* the project name in the command.
*   `-f` Reformats the parameters file 'project.PRM' according to the present version (using 'template.PRM'). Writes to stnd. output. Add '2' to shorten. Use '-f' after the projectname. 
*   `-h` To show this information file. Append 'p' to print the default parameters file in stdout; 'c' to print a clean version of it; 'u' to print an example of load file; 'i' (default) to print this help file.
*   `-l` To artificially fill with sediment all lakes except sea.
*   `-M` Sets the plate model (isost_model): 0 means no isostasy; 1 elastic plate (with `-T0` produces local isostasy); 2 for viscoelastic plate (following Nadai , 1963).
*   `-N` To set Nx and Ny of gridding. Ny<Nx optimises computation time.
*   `-o` Redirects standard output to file projectname.out.
*   `-P` Sets the plotting mode:
    *   `-P0` (default) no plotting.
    *   `-P1` produce a postscript graphic display using GMT (tisc.gmt.job) at the end of the run.
    *   `-P2` make it every time step and convert the PS files into JPG using ImageMagick. Append geometry (default is: `-P2"-trim -density 120"`, which crops the rectangle used and then saves it into JPG format with 120dpi).
    *   `-P3` use the python graphic script tisc.plot.py instead at every time step.
*   `-p` Sets the amount of runoff [l/m2/a]=[mm/a], its proportionality with height [mm/a/km], and lake evaporation. `-p0` suppresses the hydrological calculations.
*   `-Q` With this option TISC will just calculate the elastic deflection of `<file>` (3 columns: x[m],y[m],pressure[Pa]). Will write deflection to the standard output and exit the program.
*   `-q` Sets any other parameter value listed in `./tisc/doc/template.PRM` file. Note: `<value>` must be given in internal program units (usually IS).
*   `-R` Add random noise to the initial topography between -random_topo/2 and +random_topo/2. seed is an integer.
*   `-r` Sets density of environment (e), crust (c), infill (i), mantle (m) or asthenosphere (a).
*   `-S` Move block number `<b>` by `<n>` positions (use with `-F`).
*   `-s` Add 'm' to use 'mathlib' to solve the flexural equations.
*   `-T` Sets elastic thickness to `<eet_value>`. 0 means local isostasy (faster calculation).
*   `-t` Sets initial (i), final (f), increment (d), surface processes (e), viscous relaxation (v) or sediment-Block-record (r) times.
*   `-V` Means verbose mode with additional run-time prints. Add `<level>` for additional i/o information. See doc/template.PRM file for level specification.
*   `-v` Changes velocity to `<vx>,<vy>` [km/Myr] for block number `<num>` or, if num<0, to all blocks with density=-<num> (use with `-F`).

### 3.3. Direct Mode (`-Q`)

With this option TISC will just calculate the elastic deflection of `<file>` (3 columns: x[m],y[m],pressure[Pa]). Will write deflection to the standard output and exit the program. This mode is for solving the flexural response to a single load without running a full time-dependent simulation. It's useful for quick tests of isostatic models.

### 3.4. Resuming a Run (`-F`)

You can resume a run that was previously stopped. TISC periodically saves its state in a `project.all` binary file. To resume, use the `-F` flag:

```sh
tisc -F project.all
```

## 4. Input Files

Input parameter values are taken from: (a) parameters file `project.PRM`, (b) command line options, (c) other input files named `project.*` (with an uppercase extension), and (d) from the `project.all` file if `-F` option is used. (b) values override (a); (c) overrides (a) and (b); (d) overrides (a) and (c). Default units are I.S. 

All the following input files accept comment lines to be inserted everywhere. These lines will simply be ignored if a parameter name and a numerical value cannot be read at the beginning.
All input files are read in the beginning of the execution, at the initial time, except for `*.UNIT` (which are read at their respective specified times) and for `*.NDEF` (read at every time step dt_eros). 

*   **`project.PRM`**: Parameters file. This file is compulsory unless using the `-Q` option. The parameter values defined in it can be also changed with the command line options listed above. `tisc/doc/template.PRM` file is an explained example containing the default parameter values. The internal use of each parameter is also somewhat explained in the include file `tisc.h`. xmin, xmax, ymin, ymax give the coordinates of the center of the boundary cells or pixels (i.e., the location of their associated nodes). Equivalent to node-registration mode in GMT or NETCDF. The conceptual limit of the model is at xmin-dx/2, ymax+dx/2, ... For example, if the resolution of a model is needed to be reduced to dx/2,dy/2 keeping the same domain parameters, then the new Nx,Ny should not be Nx*2,Ny*2, but (Nx-1)*2+1,(Ny-1)*2+1. If (topoest != 0) then the load will rest at zero level and plate subsidence will be filled with 'densinfill' material. Do not combine it with erosion.
*   **`projectN.UNIT`**: File modifying the geometry/architecture of the model (i.e, introducing new bodies, moving blocks, or faults). `N` is the number of this load; first load file should be: `project1.UNIT`; `doc/template.UNIT` file is an explained example. Or type `tisc -hu` to get that file in the screen. `N` must be sorted consistent with the UNIT parameter 'time', so that the time in UNIT number N must be smaller than the time in UNIT N+1. TISC deformation consists of moving predefined blocks of prescribed velocity and geometry by an entire number of cells at each time step. The minimum shift is either 0 or 1 cells. If your time step dt is 0.5 Myr, velocity is 5 km/Myr and cell size is 10 km, the block will move only once every 4 time steps. So, depending on the grid cell size, the motion will become more or less abrupt. 
*   **`project.EET`**: Equivalent Elastic Thickness file (optional, if not present EET will be taken from the `*.PRM` file). Note that abrupt EET lateral changes can lead to instabilities in the calculation of flexural deflection. The file format is:
    ```text
    [z_default <value>]      #Default value of the elastic thickness.
    [mode_interp  <value>]   #See explanation in doc/template.UNIT.
    x1  y1  thickness1
    x2  y2  thickness2
    x3  y3  thickness3
    ...
    Or:
    [mode_interp  4]
    thickness1
    x1  y1
    x2  y2
    x3  y3
    thickness2
    x2  y2
    x1  y1
    ...
    ```
*   **`project.PRFL`**: This file (optional) defines a x-y polygon corresponding to the cross section that will be written in a `*.pfl` file. Default is no cross section.
*   **`project.ZINI`**: Initial topography (optional). Similar format to `project.EET`.
*   **`project.WINI`**: Initial deflection file (optional). Similar format to `project.EET`. Positive deflection downwards.
*   **`project.RAIN`**: Runoff distribution (optional). Similar format to `project.EET`. Units: mm/a. Negative values and undefined cells are interpreted as no-data and substituted by the automatic calculation from the parameters in `*.PRM`. The values in `*.RAIN` are kept constant through the entire model evolution.
*   **`project.NDEF`**: Allows to change properties of specific nodes. The closest node to x,y will be changed: 'param' determines what parameter is changed at that node: 1: to add water discharge at x,y (m3/s, max. is Amazon at 120e3). 2: to add sediment discharge (kg/s, max. is Ganges at 52e3). This file, in contrast to all others, is read not once at startup but once for every surface process iteration.
*   **`project.SLV`**: Evolution of the sea level and the level separating erosion and sedimentation levels (optional, read only if water_load is activated in the PRM file):
    ```text
    time1   s_level1    [eroslevel1]
    time2   s_level2
    time3   s_level3    [eroslevel2]
    ```
*   **`project.INS`**: Evolution of the insolation (optional, use Laskar's, f.e.). Only for hydro_model=1. Adds a factor insolation/mean_insolation multiplying `<rain+Krain*z>` (see rain_amp in `*.PRM`) and also multiplying the discharges in `*.NDEFS`. Thus, absolute insolation values or their units are irrelevant. The deviations of insolation from its average value will be amplified in rain by a factor `<rain_amp>` (if not 0).
*   **`project.REC`**: Horizon recording times file (optional, defaults to dt_record parameter, see `*.PRM` file).
*   **`project.RIV`**: Initial river paths that must be respected by the initial topography. When introducing a Digital Elevation Model through the `*.ZINI` file, the grid may not respect the path of the rivers because of its limited resolution, which generates unrealistic barriers. This file must be introduced in order to allow the program to modify the topographic grid and fit it to the drainage network. Introduce river in x-y format starting from the highest point. Different segments can be separated by lines starting with '>'.
*   **`projectN.VISC`**: Only required if a thin sheet unit has been defined. Contains the viscosity x,y distribution for unit `N`, but only its thermal component (~ .1e7). The total viscosity will be calculated as viscosity[Pa*s] = viscTer/str.rate. Format is equivalent to file `project.EET`. 
*   **`project.TSBC`**: Only required if defining a thin_sheet unit. Boundary conditions of the thin sheet velocity or stress field for the units is defined to undergo a thin-sheet-like deformation.
*   **`project.CMP`**: This file defines x-y polygons (in km) separated with a '>', used only to be plotted in the postscript image, called from graphic script like script/tao.gmt.job, just for comparison with model results. No effect on the model results.
*   **`project.INT`**: This file defines a x-y polygon (in km) of interest where certain statistics should be calculated during the execution of graphic scripts, such as the volume of sediment. `project.INT2` `project.INT3` etc are also read. Used only during the image script (`script/tisc.gmt.job`) after tisc, just post-processing. No effect on the model results.

### 4.1. Interpolation Modes
Most TISC input files define 3D surfaces or 2D planform distributions of a variable. To facilitate this, TISC accepts different ways of interpolation between given points:

*   **`0`**: 3-column file containing x,y,z in 'Nx x Ny' skyline-sorted rows. x,y are ignored (no interpolation performed). Useful to pass a surface generated externally with SURFER.
*   **`1`**: 3-column file containing x,y,z. Inverse distance interpolation with the points supplied.
*   **`2`**: 3-column file containing x,y,z. Inverse square distance interpolation with the points supplied.
*   **`3`**: 3-column file containing x,y,z. The nearest given point is assigned each cell.
*   **`4`**: Polygon interpolation mode (see format below).
*   **`5`**: Binary (short int) skyline array of z values (no interpolation performed).
*   **`6`**: 1-column file containing z-values in 'Nx x Ny' skyline-sorted rows (no interpolation performed).
*   **`7`**: Same as 4, but nodes falling out of all polygons are interpolated with the distance to each polygon (no default value assigned).
*   **`8`**: 3-column file containing x,y,z. Each cell gets a value if some point is given inside, otherwise it gets the default value. If more than one value is given then they will be averaged.

If `mode_interp` is set to `4` or `7`, then x-y polygons are read with an associated z-value. In this case the format is:
```text
mode_interp  4
z_value_pol_1	[n|c]
x_pol1_1      y_pol1_1
x_pol1_2      y_pol1_2
x_pol1_3      y_pol1_3
...
z_value_pol_2	[n|c]
x_pol2_1      y_pol2_1
x_pol2_2      y_pol2_2
x_pol2_3      y_pol2_3
...
```
For `mode_interp=4` or `7`, if `c` is specified for a polygon, the z value inside that polygon will be constant. If `n` is specified (default), the z value interpolated for all nodes falling inside that polygon but outside the next polygon will vary linearly with the distance to both polygons. So, tricky as it seems, the order of the polygons matters!

## 5. Output Files

The screen output while TISC runs shows a number of statistics including the balance of water, sediment and vertical force over the entire model domain. Total rain in the model must equal total evaporation plus water exiting the model through the boundaries. Except for the effect of the boundary conditions, the flexural restitutive force must equal the total load. This can be checked in the screen output at each time step, for example:

```text
T= 9.5000 My
  2D prof. :  max =     799.8 m     min =    -531.8 m
  topogr.  :  max= 2179.1 m	min= 0.0 m	vol=1.18e+12 m3
  load: -2.28e+17 N   restit_force: -9.68e+15 N
  deflect. :  max= 114.8 m	min= -18800.3 m	vol=-7.25e+12 m3
  precipit :  max= 0.3 m/yr	min= 0.1 m/yr	integr=5.40e+10 m/yr*m2
  evaporat.:  max= 0.0 m/yr	min= 0.0 m/yr	integr=0.00e+00 m/yr*m2
  rain_now : +2.75e+02 m3/s  evap_wat: +0.00e+00 m3/s outp_water: +2.75e+02 m3/s
  lake 1/3 : 1.68e-04 km3 3.24e+00 km2   10 m  -22,20  1 out @ -23,20 40.5 m3/s
  river max:    54.72 m3/s     9.46 kg/s @  -26.0,-20.8 km, 0.1 m
  topo_diff_eros_max=    -1.10 mm/yr    sedim_max:     0.53 mm/yr
  eros_nosd: +9.44e+15 N     sedim_inc: -2.91e+13 N   outp_seds:   +9.47e+15 N  
  eros-sed.:  max= 1.62 mm/yr	min= -0.53 mm/yr	difference=1.67e+06 m3/yr
  sediment vol: 4.09e+03 km3
```

At the final time step, a list of the generated Blocks is shown including their properties:  
```text
2 Blocks:
No. Density Age Volume Surf.  Vel_x Vel_y Shft_x Shft_y erosL
     kg/m3  My  1e3km3 1e3km2 km/My km/My   km     km    [m] 
 1:  2200   5.0    0.0   0.6   0.0   0.0    0.0    0.0 6.0e+04  S
 0:  2200   0.0    0.0   0.8   0.0   0.0    0.0    0.0 6.0e+04  S  1st Block
 -:  2850   0.0    -     -     0     0      0      0   1.2e+05  -  basement
```

Output files contain the information required for graphic output, and have a lowercase extension:

*   **`project.hrz`**: Final geometry of all Blocks (elevat./depth of all the horizons).
*   **`project.xyzt`**: Isostatic (flexural) subsidence (deflection w).
*   **`project.all`**: Binary file with all the information to resume the model (see `-F` option). Valid for checkpointing purposes.
*   **`project.xyw`**: Drainage network. Mass transport rate. Type of node ('R'=river; 'L' for lake; 'E' for lake outlet; 'S' for sea).
    *   Column #1 x(km) 
    *   Column #2 y(km) 
    *   Column #3 water(m3/s) 
    *   Column #4 sed[kg/s] 
    *   Column #5 type 
    *   Column #6 topo[m] 
    *   Column #7 x-to 
    *   Column #8 y-to 
    *   Column #9 topo-to 
    *   Column #10 preciptation (runooff) [mm/y] 
    *   Column #11 lake/sea evaporation[mm/y] 
    *   Column #12 ice thickness[m] 
    *   Column #13 ice sediment load[m]
    *   Column #14 maximum swimming distance between neigbor cells [km]
*   **`project.bas`**: Similar to the previous, but sorting the nodes hierarchically according to hydrographic (sub)basins. Within each basin, nodes are sorted downstream (every entry drains to another entry listed later in the file.
    *   Column #1 x[km] 
    *   Column #2 y[km]
    *   Column #3 dischg[m3/s]
    *   Column #4 masstr[kg/s] 
    *   Column #5 type
    *   Column #6 topo[m]
    *   Column #7 length[km] Length along channel from river mouth.
    *   Column #8 x-to
    *   Column #9 y-to
    *   Column #10 topo-to[m]
    *   Column #11 length-to[km]
    *   Column #12 eros[m/My]
    *   Column #13 level
    *   Column #14 nodes_above
*   **`project.lak`**: Lake information.
*   **`project.st`**: Total erosion/sedimentation and present rate (in equivalent meters of denscrust).
*   **`project.eeth`**: Interpolated values for elastic thickness entered via `*.EET`.
*   **`project.vel`**: Velocity and viscosity distribution from the thin-sheet Block.
*   **`project.pfl`**: Cross section along path defined in `*.PRFL`.
*   **`project.ps`**, **`project.jpg`**, **`project.svg`**: Graphical plots. A postscript produced with GMT 4.0 software (using unix scripts: `./tisc/script/tisc.*.gmt.job`), with graphics of results of TISC. The GMT and shell commands generating this image are always searched by TISC in a file named `tisc.gmt.job` in the user's path. The SVG/JPG outputs are made by the python script `script/tisc.plot.py`, using matplotlib.

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

## 3. Computational Performance and Scaling

The computational time of TISC is heavily dependent on the grid dimensions ($N_x$ and $N_y$). If we assume a square grid where $N = N_x = N_y$, the total number of grid cells is $V = N^2$. The two most expensive operations scale as follows:

### 3.1. Isostasy (Flexure Solver)
The elastic and viscoelastic flexure equations are solved using a finite-difference band solver.
*   **Time Complexity:** The solver processes a matrix of size $V \times V$. Because of the 2D finite difference stencil, the matrix is banded with a bandwidth $B \approx \mathcal{O}(\min(N_x, N_y))$. The time consumption for a band solver scales as $\mathcal{O}(V \cdot B^2)$. Substituting the dimensions, this yields **$\mathcal{O}(N_x N_y \cdot \min(N_x, N_y)^2)$**.
*   For a square grid ($N_x = N_y = N$), the flexure solver time scales as **$\mathcal{O}(N^4)$**.
*   **Optimization Tip:** Because the solver time is proportional to the square of the smallest dimension, keeping $N_y < N_x$ (or vice versa) dramatically reduces computation time while maintaining the same total number of cells (e.g., a $1000 \times 100$ grid solves exponentially faster than a $316 \times 316$ grid).

### 3.2. Surface Processes
The surface transport algorithms (`surf_proc.c`) normally scale linearly or log-linearly, but gracefully degrade into quadratic limits in complex topographies.
*   **Base Complexity:** Sorting the grid and routing water on standard terrain takes $\mathcal{O}(V \log V)$, which translates to **$\mathcal{O}(N^2 \log N)$**.
*   **Worst-Case Complexity:** Perfectly flat areas (such as lake beds) require topological sorting to ensure proper sediment routing, and evaporating endorheic lakes require recursive node-deletion and component checking (`Divide_Lake`). Both of these sub-algorithms feature nested loops that scale quadratically with the size of the flat area or lake ($A$). In the absolute worst-case scenario (e.g., a massive shrinking lake covering the entire grid, $A \approx V$), the time consumption degrades to $\mathcal{O}(V^2)$, which is **$\mathcal{O}(N^4)$**.

### 3.3. Memory (RAM) Consumption
The memory footprint of TISC is generally dominated by the arrays required for the flexure matrix solver.
*   A standard band matrix requires storing $V \times B$ elements.
*   This results in a memory complexity of **$\mathcal{O}(N_x N_y \cdot \min(N_x, N_y))$**.
*   For a square grid ($N_x = N_y = N$), RAM usage scales as **$\mathcal{O}(N^3)$**.
