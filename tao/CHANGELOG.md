# tAo Changelog

## tAo_2026-03-14
* Full refactorisation of the code using Gemini. Optimisation. 

## tAo_2025-11-17
* Python graphic output available through `script/tao.plot.py`.

## tAo_2019-12-06
* Gradual sediment compaction now implemented (see new parameter `compact_depth`).
* For `isost_model>=3`: yield stress is now written in `*.strs` output file. 
* Main graphic script `tao.gmt.job` is now accepting a number as 2nd argument, to predefine the x coordinate for the stress distribution and envelope.

## tAo_2017-05-19
* Average erodibility is now calculated for the shallowest meters of bedrock.
* Erosion boundary conditions 0 and 1 are swapped.
* `-v` command line option to change velocity of Blocks based on density.
* Changed `UNIT1D` and Units structure to `BLOCK1D` and Blocks, to avoid confusion with input files (`*.UNIT`).

## tAo_2014-09-21
* Erosion models 6 and 7 implemented.
* Graphic scripts now distinguish old sediment in cross section.
* Improved sediment deformation. `deform_sed` parameter added. Sediment above moving Blocks is transported and piled in front.
* Change in units of precipitation and sed. load. Incorporation of `riverbasinwidth`.

## tAo_2012-08-29
* Change in output extensions that differ from inputs only in their case, such as `*.eet` -> `*.eeth` and `*.yse` -> `*.ysen`. This is to avoid conflict in Mac computers. Successfully compiled for Mac IOS 10.7.5.
* Solved bug during `-F` (project resume) option.

## tAo_2009-12-07
* Additional erosion models implemented.

## tAo_2006-11-02
* Removed bugs in `-Q` mode, now pressure is entered instead of thickness.

## v4.0 (2002)
* Fluvial sediment transport. 
* Loads and faults now defined in `*n.UNIT` input files (instead of `*n.CRG` file), and the file accepts now a list of parameters for the new Block. 
* Files `*.H0` and `*.W0` are now called `*.ZINI` and `*.WINI`, respectively.
* (2004) Orographic precipitation implemented following a conservative scheme (`runoff=precipitation-evaporation`).
* (2005) A memory bug solved in surface_transport (this was producing NaN results in a quite unpredictable manner).

## v3.4 to v1.0 (1993 - 1999)
* **v3.4 (1999):** Change in param. file keeping compatibility. Possible to resume a model that was stopped (option `-F`). Option 'c' in loads to cut Blocks with faults.
* **v3.3 (1998):** Stress history is accounted for.
* **v3.21 (1997):** Allows `sea_level` variations along time. Allows erosion/sedimentation variations along time.
* **v3.2 (1996):** Inclusion of thrust loading and sediment deformation. Command-line options.
* **v3.0 (1995):** Multilayered elastic-plastic rheology. Geoid anomaly calculations.
* **v2.0 (1994):** Inclusion of moving loads and surface processes. Gravity anomaly calculations.
* **v1.0 (1993):** Proto-model or numerical algorithm to calculate thin plate elastic and viscoelastic flexure.
