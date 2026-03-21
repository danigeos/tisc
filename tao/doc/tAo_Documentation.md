# tAo Documentation

This document provides a comprehensive user manual and technical documentation for the tAo (1D/Cross-section Tectonics, Isostasy, and Surface Processes) model.

## Table of Contents
* Part 1: User Manual
  * 1. Introduction
  * 2. Installation
  * 3. Running tAo
  * 4. Input Files
  * 5. Output Files
* Part 2: Technical Documentation
  * 1. Code Architecture
  * 2. Core Physical Models
* Part 3: Additional Information
  * 1. Graphic Output
  * 2. Algorithms and Numerical Methods
  * 3. Glossary
  * 4. Authorship & Credits
  * 5. Acknowledgements
  * 6. Disclaimer

---

# Part 1: User Manual

## 1. Introduction

The program `tAo` is a pseudo-2D (vertical cross-section / 1D profile) numerical model designed to calculate the interaction between lithospheric flexure, tectonic loading (e.g., by thrusting), and surface transport (erosion and deposition). The code puts a particular emphasis on the formation of foreland basin systems.

It allows the application of different rheological behaviors for the lithospheric plate:
* Pure Elastic.
* Viscoelastic.
* Depth & temperature-dependent elastic-plastic rheology.

**This program is particularly designed to:**
* Calculate lithospheric flexure directly from the lithospheric stratification and load distribution (option -Q).
* Model foreland basin formation and the large-scale sedimentary infill geometries.
* Solve simple 1D flexural problems in other geological contexts such as trenches, seamounts, etc.
* Also to be used externally (from a unix script) as a flexure calculator, as part of a broader job.
* Calculate surface processes (erosion/sedimentation) in 1D.

## 2. Installation

### 2.1. Dependencies
To compile and run tAo, you will need:
*   A C compiler (e.g., `gcc` or `clang`).
*   `make` utility.
*   (Optional) GMT 4 or Python 3 for graphical output.

### 2.2. Compilation
tAo and TISC share a unified build system. Run the `make` command in the root directory (where `config.mk` is located):
```sh
make
```
This will compile the source code and place the `tao` executable in the `bin/` directory.

### 2.3. Environment Setup
You must add the `bin/` and `script/` directories to your system `PATH` and set the `tao_dir` environment variable.

For **bash/zsh**:
```sh
export tao_dir=/path/to/your/tao/folder
export PATH=$PATH:$tao_dir/bin:$tao_dir/script
```

## 3. Running tAo

### 3.1. Basic Syntax
The standard execution of tAo requires a project name corresponding to the `*.PRM` file:
```sh
tao [options] project_name
```

### 3.2. Command-Line Options
Command-line options override variables initialized inside the `*.PRM` file.

```text
tao  project  -A[1|2] -B<bound_type> -D[x0/xf] -d<dx> -F[file] -f[2] 
      -h[i|u|p]
      -M<lih_type>[t] -m<app_mom> -N<Nx> -o -P[0|1|2|3][geom] -p<tec_force>
      -q<param=value> -Q<file>  -r<a|c|i|m|a><density> -S<b>/<n> 
      -s<app_force> -T<eet> -t<i|f|d|v|r><time> -V[<level>] 
      -v[<num>/<vel>]
```

**Options:**
*   `project` is the root name for the project files (i.e., if project=test then the parameters file should be test.PRM and the first load test1.UNIT.
*   `-A` Switches gravity calculations on (1 for Bouguer; 2 for Free Air).
*   `-B` Sets boundary conditions type (0 to 5, see template.PRM file).
*   `-D` Set the model region (domain: left and right coordinates) [m]. Modifies parameters x0,xf,xmin,xmax.
*   `-d` Sets the x increment used for finite differences method.
*   `-F` To continue (resume) a previous model reading `<file>`. `<file>` is a binary file written by tAo as 'project.all' (default file is 'project.all'). Put this option always *after* the project name in the command.
*   `-f` Reformats the parameters file 'project.PRM' according to the present version (using 'template.PRM'). Writes to stnd. output. Add '2' to shorten.
*   `-h` To show this information file. Append 'p' to print the default parameters file in stdout; 'c' to print a clean version of it; 'u' to print an example of load file; 'i' (default) to print this help file.
*   `-M` Sets the plate model. (0 to 6, see isost_model parameter template.PRM). Add 't' to supress accounting for stress history.
*   `-m` Sets the boundary moment applied at left boundary. >0 means clockwise.
*   `-N` Sets the number of lateral gridding points (Nx).
*   `-o` Redirects standard output to file projectname.out.
*   `-P` Sets the plotting mode:
    *   `-P0` (default) no plotting.
    *   `-P1` produce a postscript graphic display using GMT (tao.gmt.job) at the end of the run.
    *   `-P2` make it every time step and convert the PS files into JPG using ImageMagick. Append geometry (default is: `-P2"-trim -density 120"`, which crops the rectangle used and then saves it into JPG format with 120dpi).
    *   `-P3` use the python graphic script tao.plot.py instead at every time step.
*   `-p` Sets the tectonic horizontal force in the plate. >0 means compressive.
*   `-q` Sets a value for any parameter listed in the parameter file `*.PRM`. Internal tAo units are expected (usually IS), rather than those units used in the PRM.
*   `-Q` Quick flexural calculation. With this option tAo will just calculate the elastic deflection in response to `<file>` (2 columns: x[m],pressure[Pa]). Will write deflection to the standard output and exit the program. With this option, tAo doesn't need a projectname, since no other files are read nor written.
*   `-r` Sets density of environment (e), crust (c), infill (i), mantle (m), or asthenosphere (a).
*   `-S` Move block number `<b>` by `<n>` positions (use with `-F`).
*   `-s` Sets the applied vertical shear force at left boundary. >0 means upwards.
*   `-T` Sets elastic thickness to `<eet_value>`. For the case of multilayered rheology plate, eet_value=-1 indicates that the initial Te must not be constant or read from the `*.EET` file, but from the `*.eeth` file.
*   `-t` Sets initial (i), final (f), increment (d), viscous relaxation (v) or sediment-block-record (r) times. Time is in My and goes from negative to positive values (Timeini<Timefinal).
*   `-V` Means verbose mode with additional run-time prints. Add `<level>` for additional i/o information. See doc/template.PRM file for level specification.
*   `-v` Changes velocity to `<vel>` for block number `<num>` or, if `<num><0`, to all blocks with density=-<num> (use with `-F`).

### 3.3. Runtime Signals
When running, the program prints single-character signals indicating the current calculation loop:
* `l`: calculating new load.
* `b`: defining a new internal Block.
* `e/s`: calculating erosion/sedimentation.
* `ft`: calculating fluvial transport.
* `e`: calculating elastic deflection.
* `v`: calculating viscoelastic deflection.
* `Rh`: Moment-curvature iteration for multi-layered plates.

## 4. Input Files

Input files are named with an uppercase extension. Inputs are read sequentially from: (a) parameters file `project.PRM`, (b) command line options, (c) other input files named `project.*` (files starting with the project name with an uppercase extension), and (d) from the `project.all` file when `-F` option is used. (b) values override (a); (c) overrides (a) and (b); (d) overrides all. Default units are I.S. 

*   **`project.PRM`**: Parameters file. This file is compulsory unless using the `-Q` option. The parameter values defined in it can be also changed with the command line options listed above. `tao/doc/template.PRM` file is an explained example containing the default parameter values. The internal use of each parameter is also somewhat explained in the include file `tao.h`.
*   **`projectN.UNIT`**: File modifying the geometry/architecture of the model setup (adds tectonic Blocks, vertical loads, moving blocks, faults, etc). `N` is the number of this load: The first load file should be named `project1.UNIT`. The time of each UNIT (`doc/template.UNIT`) must increase with `N`. tAo deformation consists of moving predefined blocks of prescribed velocity and geometry by an entire number of cells at each time step. The minimum shift is either 0 or 1 cells. If your time step dt is 0.5 Myr, velocity is 5 km/Myr and cell size is 10 km, the block will move only once every 4 time steps. So, depending on the grid cell size, the motion will become more or less abrupt. 
*   **`project.ZINI`**: Initial height (optional).

Lithosphere rigidity or elastic thickness (for the isostatic calculations) can be entered in 4 ways: 1) Constant EET provided in `*.PRM`; 2) Laterally varying EET defined in `*.EET`; 3) EET calculated from an input yield stress envelope in `*.YSE`; and 4) EET calculated from the temperature distribution in the lithospheric plate and the crustal geometry (`*.TMP`, `*.CRUST`, `*.UCRUST`). 1) and 2) apply only when isost_model=1,2. 3) and 4) apply for isost_model>=3. 

*   **`project.EET`**: Elastic thickness file (optional, defaults to default_Te parameter, see `*.PRM` file). Also used as initial Te for the multilayered plate case.
    ```text
    x1  thickness1
    x2  thickness2
    x3  thickness3
    ...
    ```
    With: x3>x2>x1. 
*   **`project.YSE`**: YSE file (only in case isost_model>=2; overrides YSE calculations using temperature). Format (note different units relative to the output file `*.ysen`):
    ```text
    z1[m]  yieldcompres1[Pa]  yieldextens1[Pa]
    z2     yieldcompres2      yieldextens2
    ...
    ```
*   **`project.TMP`**: Temperature file (only in case isost_model>=2). If a YSE input file is given, then this TMP file is ignored.
    ```text
    x1[m]   #this x is optional, if not present, only one geotherm will be read 
    z1[m]   temperature1[C] #(z is depth, positive downwards)
    z2      temperature2
    ...
    x2[m]   #x location for the second geotherm
    z1[m]   temperature1[C]
    z2      temperature2
    ```
    If a single number is read in a line, it will be taken as the x position for the z-T geotherm listed below, and then additional geotherms at other x locations can be supplied. Temperature will be linearly interpolated (laterally and in depth) in between them.
*   **`project.CRUST`**: Crust thickness. Used only if `*.TMP` is provided. If this file is not present, then the parameter 'crust' in `*.PRM` will be adopted (constant along x).
*   **`project.UCRUST`**: Upper crust thickness Used only if `*.TMP` is provided. If this file is not present, then the parameter 'ucrust' in `*.PRM` will be adopted (constant along x).
*   **`project.SLV`**: Sea and erosion/sedimentation levels along time file (optional):
    ```text
    time1   level1  [eroslevel1]
    time2   level2
    ...
    ```
*   **`project.WINI`**: Initial deflection file (optional). Positive deflection downwards.
*   **`project.REC`**: Horizon recording time file (optional, defaults to an interval of dt_record, defined in the `*.PRM` file):
*   **`project.CMP`**: An x-z file to be compared with results (optional). The results with real profiles digitized in this file. This file is not used by tAo but by the GMT job `tao.gmt.job` which produces a postscript image of it. Similar format (two columns) but x is expected in km (z in m).

All these files are interpolated (except the `.CMP` one) and accept comment lines to be inserted everywhere. These lines will be simply ignored when two floats cannot be read in them. Attention: x3>x2>x1. tAo will linearly interpolate between the given points.

## 5. Output Files

Outputs are written to files that start with the project name with an lowercase extension. Please be careful because all these files are removed when running again the same project with tAo. Note that x units in the output files are in km.

*   **`project.xzt`**: An output file with the deflection and the topography at each loading event is written.
*   **`project.pfl`**: A final profile of the model containing the final model geometry (elevation of the top and bottom of each Block, from basament to surface).
*   **`project.eros`**: Erosion (>0) and sedimentation (<0) at each position. Note that each value is the result of adding the thickness of eros/sedim at each time step at that particular location. A final value of 1000 m of erosion may be the result of 1000 m deposition + 2000 m erosion. Besides, this is recorded where it happened, so the values of eros./sedim. are not shifted with the movement of the sediment or mother rock units.
*   **`project.xg`**: A file with the gravity anomaly.
*   **`project.ps`**, **`project.jpg`**, **`project.png`**: Graphical plots. A multipage postscript produced with GMT 4.0 software (using a unix script: `tao.gmt.job`), with graphics of all the results of tAo: Distribution of loads/units/sediment and the deflected basament compared with the `project.CMP` file if it exists. Bouguer gravity anomaly and geoid anomaly if requested. Stress distribution when rheological considerations have been requested. Yield stress envelope and stress distribution at the maximum moment point. The GMT and shell commands generating this image are always searched by tAo in a file named `tao.gmt.job` in the user's path. Python graphic outputs utilize `tao.plot.py`.
*   **`project.ysen`**: Calculated yield stress envelope and stress distribution at the maximum moment point. Only if a multilayered plate was selected (isost_model>2).
*   **`project.strs`**: (only if isost_model>2) x-z cross-sectional grid of stress, temperature, and yield constraints.
*   **`project.temp`**: Temperature distribution within lithosphere (only isost_model>2 and with swith_verbose=1).
*   **`project.eeth`**: Elastic thickness distribution resulting calculated when isost_model>2. If this file is present and the default Te value is a signal `-1`, then this `*.eeth` file is used as the initial value of the moment-curvature iteration.
*   **`project.all`**: Binary file to be read when the model is going to be resumed (continued) with option `-F`.

---

# Part 2: Technical Documentation

## 1. Code Architecture

The tAo codebase shares much of its numerical heritage with its 3D counterpart (TISC) but is specifically optimized for 1D arrays and vertical cross-sectional deformation.

*   `tao.c`: The core driver. It manages memory allocation, command-line parsing, the primary time loop, and the sequence of physics calls.
*   `taoio.c`: Handles all specialized file input/output, parsing strings, and formatting grid outputs to disk.
*   `taolib.c`: Contains the numerical subroutines, flexural matrix construction, yield stress calculations, and geometric block management. 
*   `taosp.c`: Manages the 1D Surface Processes (erosion, sedimentation, fluvial routing).

## 2. Core Physical Models

The conceptual model assumes the classical thin plate approach, computing subsidence through elastic, viscoelastic, or depth-dependent-elastic-plastic rheologies.

### 2.1. Lithospheric Flexure

The 1D flexural equation is constructed as a finite difference matrix and solved using a band solver. tAo supports highly advanced rheological modeling via the `isost_model` parameter:

*   **0:** No isostasy.
*   **1:** Pure elastic plate model (constant or spatially variable `Te`).
*   **2:** Visco-elastic plate model (relies on the `tau` parameter for relaxation time).
*   **3:** Elastic-plastic oceanic lithosphere.
*   **4-6:** Elastic-plastic continental lithosphere (decoupled, never decoupled, or auto-coupling).

#### Elasto-Plastic Rheology (`isost_model` >= 3)
For advanced modes, the effective elastic thickness is dynamically calculated based on the lithosphere's thermal gradient (`project.TMP`) and crustal thicknesses (`project.CRUST`). 

The code calculates a Yield Stress Envelope (YSE) that restricts the maximum stress that can be sustained at a given depth. A moment-curvature iterative algorithm (`Rheo_Flex_Iter`) organically decreases the local rigidity (`D`) of the plate wherever the bending stresses exceed the calculated yield stress, mimicking brittle failure (faulting) in the upper crust and ductile flow in the lower crust/mantle.

### 2.2. Surface Processes

tAo simulates erosion and deposition across the 1D profile using parameter `erosed_model`:
*   **Mode 1 (Diffusion & Constant Rate):** Background weathering/sea-sedimentation and basic hillslope mass transport proportional to the topographical slope. 
*   **Mode 2 (Beaumont Undercapacity):** Simulates a river where sediment transport capacity (`qeq`) is proportional to discharge and slope. The erosion/deposition rate reacts linearly to the difference between actual load and capacity.
*   **Mode 3 (Tucker & Slingerland):** Uses hybrid stream power laws.
*   **Mode 6 (Basal Shear Stress):** Solves detachment-limited and transport-limited regimes simultaneously, dynamically scaling with an explicit basin width (`riverbasinwidth`) to account for upstream 3D accumulation in a 2D profile.

Sediments are mechanically tracked, deposited as distinct geological `Blocks`, and automatically dragged and folded along with any underlying moving tectonic thrust faults defined in the `.UNIT` files.

---

# Part 3: Additional Information

## 1. Graphic Output
If the `-P` flag is utilized, tAo produces a multi-panel visual summary of the simulation at each time step.
*   If GMT is requested (`-Pc`), tAo calls the `tao.gmt.job` shell script to generate highly customizable PostScript (`.ps`) files showing layer geometries, gravity anomalies, and stress envelopes.
*   If Python is requested (`-Pp`), it will utilize `tao.plot.py` to natively render `.png` image sequences.

## 2. Algorithms and Numerical Methods
An explicit finite difference scheme is used to solve the spatial partial differential equations as well as the time-dependent processes. Matrices arising from the differential equations are solved using a highly optimized classical band solver.

Gravity anomalies are calculated using the classic Talwani et al. (1959) algorithm. Geoid anomalies are calculated using Chapman's (1979) methodology.

## 3. Glossary
*   **Hidden load:** Subduction, lithospheric mantle thickening, plate interaction forces, etc., leading to the need of additional loads to fit the observed deflection.
*   **Project:** Group of input files corresponding to a work project.
*   **switch_topoestable:** Strongly simplifies the modeling by having the read load stay at zero level and the space filled up with a constant density 'densinfill'.
*   **Topographic load:** Remains on the top of the profile and neither it nor the whole profile subsides.
*   **Unit:** Input files (`*.UNIT`) to provide geometry of faults or loads.
*   **Blocks:** Bodies written in `*.hrz`. A Block corresponds to a 2D body with a density, velocity, age and thickness profile associated.

## 4. Authorship & Credits
tAo is a freeware program developed since Sept. 1993 by Daniel Garcia-Castellanos. Details on copyright are in `tao/doc/LICENSE.txt`.

When showing results from this code please cite at least one of the following articles:
*   Garcia-Castellanos, D., M. Fernandez & M. Torne, 1997. Numerical modeling of foreland basin formation: a program relating thrusting, flexure, sediment geometry and lithosphere rheology. Computers & Geosciences, 23 (9), 993-1003.
*   Garcia-Castellanos, D., 2007. The role of climate in high plateau formation. Insights from numerical experiments. EPSL, 257, 372-390.

## 5. Acknowledgements
The present code has been developed mainly in the Instituto de Ciencias de la Tierra Jaume Almera (Barcelona) and the Vrije Universiteit (Amsterdam). Feedback with users have been key and special mention must be made to: J.M. Gaspar-Escribano, Brian Shiro, Erica Emri.

## 6. Disclaimer
tAo is a tool developed for my own research with no professional programmers involved, the code is thus not so user-friendly, and I cannot always guarantee my assistance in case of trouble.