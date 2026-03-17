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
Run the `make` command in the `tao` root directory:
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

*   `-A[1|2]`: Switches gravity calculations on (1 for Bouguer; 2 for Free Air).
*   `-B<bound_type>`: Sets boundary conditions type (0 to 5).
*   `-D[x0/xf]`: Sets the model region domain bounds in meters.
*   `-d<dx>`: Sets the `x` spatial increment used for the finite differences grid.
*   `-F[file]`: Resumes a previous model reading from a `.all` file.
*   `-f[2]`: Reformats the parameters file `project.PRM` to match the current version syntax. Add `2` to shorten.
*   `-h[i|u|p|c]`: Show help menu. `p` prints the default PRM, `c` prints clean PRM, `u` prints a sample UNIT file.
*   `-M<lith_type>[t]`: Sets the plate model (0 to 6). Add `t` to suppress accounting for stress history.
*   `-m<app_mom>`: Boundary moment applied at the left boundary (>0 means clockwise).
*   `-N<Nx>`: Sets the number of lateral gridding points in the 1D array.
*   `-o`: Redirects standard output to `projectname.out`.
*   `-P[c[geom]|p]`: Graphic output. `p` uses Python script `tao.plot.py`. `c` uses GMT.
*   `-p<tec_force>`: Horizontal tectonic force (>0 means compressive).
*   `-q<param=value>`: Overrides any specific parameter value listed in the PRM file.
*   `-Q<file>`: Quick flexural calculation mode. Calculates the elastic deflection of a supplied load file and exits immediately.
*   `-r<a|c|i|m|e><density>`: Sets density of asthenosphere (a), crust (c), infill (i), mantle (m), or environment (e).
*   `-S<b>/<n>`: Moves block number `b` by `n` positions (used with `-F`).
*   `-s<app_force>`: Applied vertical shear force at the left boundary (>0 means upwards).
*   `-T<eet>`: Sets elastic thickness. `-1` means initial Te will be read from the `*.eeth` file instead.
*   `-t<i|f|d|v|r><time>`: Sets times: initial (i), final (f), increment (d), viscous relaxation (v), or record (r).
*   `-V[<level>]`: Sets the verbose level (0-5) for standard output prints.
*   `-v[<num>/<vel>]`: Changes velocity for block number `<num>` (if `<num>` < 0, applies to all blocks with density = `-num`).

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

All inputs are read sequentially from: (a) `project.PRM`, (b) command line options, (c) other files named `project.*` with an uppercase extension, and (d) from the `project.all` file when the `-F` option is used.

*   `project.PRM`: The master parameters file.
*   `projectN.UNIT`: Defines tectonic architecture, moving blocks, thrusts, and vertical loads. Numbered sequentially.
*   `project.ZINI`: Defines the initial bathymetry/topography.
*   `project.EET`: Variable Elastic Thickness spatial distribution file.
*   `project.YSE`: Direct Yield Stress Envelope file (overrides temperature-based rheology).
*   `project.TMP`: Temperature file containing vertical geotherms.
*   `project.CRUST`: Variable crust thickness.
*   `project.UCRUST`: Variable upper crust thickness.
*   `project.SLV`: Time-series of Sea Level and erosion/sedimentation baselines.
*   `project.WINI`: Initial baseline isostatic deflection.
*   `project.REC`: Specific timeline horizon recording triggers.
*   `project.CMP`: An external x-z polyline used solely to compare tAo outputs graphically against real-world observations.

## 5. Output Files

Output files are written with lowercase extensions.
*   `project.pfl`: Final unified profile of the model containing the horizons of all blocks from basement to surface.
*   `project.xzt`: Deflection and topography states captured at each horizon recording event.
*   `project.eros`: Cumulative erosion (>0) and sedimentation (<0) recorded at each point. 
*   `project.xg`: Computed gravity anomalies.
*   `project.ps`, `project.jpg`: Graphical plots.
*   `project.ysen`: Calculated yield stress envelope and stress distribution at the maximum moment point.
*   `project.strs`: x-z cross-sectional grid of stress, temperature, and yield constraints.
*   `project.temp`: Output temperature distribution.
*   `project.eeth`: Effective elastic thickness distribution calculated dynamically when `isost_model` > 2.
*   `project.all`: Binary snapshot of the run, used for resuming the model.

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