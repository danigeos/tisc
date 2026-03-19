# TISC Changelog

## TISC_2026-03-14_factorized
* Full refactorisation of the code using Gemini. Optimisation and 
slight parallelisation. Complex endorheic-lake systems now solved about 20 times faster.
* Grain size tracking implementation (one average value). 
* Evaporite (Gypsum and Halite) salt precipitation implementation.

## TISC_2025-11-06
* Precipitation proportional to insolation file in `*.INS`.
## tAo_2025-11-17
* Python graphic output available through `script/tisc.plot.py`.

## TISC_2023-08-31
* Debugged in `tisclib.c` line 863: Initial water column was loaded as initial isostatic load, by Hanneke Heida. 
* Added input parameters for harmonic changes in precipitation: `rain_amp`, `rain_per`. 
* Repared bug in the screen output on lake properties.
## tAo_2019-12-06
* Gradual sediment compaction now implemented (see new parameter `compact_depth`).
* For `isost_model>=3`: yield stress is now written in `*.strs` output file. 
* Main graphic script `tao.gmt.job` is now accepting a number as 2nd argument, to predefine the x coordinate for the stress distribution and envelope.

## TISC_2019-05-03
* Underground water flow now implemented (see new parameter `permeability`).
* Sklar & Dietrich's (1994) erosion dependence on sediment load added by Michael Berry (UNM, Socorro), new `erosed_model` #8. 
* Critical stress for sediment mobilisation.

## TISC_2018-07-13
* Gradual sediment compaction now implemented (see new parameter `compact_depth`). 
* Bug solved for `erosed_model` 6 (`surf_proc.c` around line 1120), erosion was overestimated (or erodability underestimated) by about a factor 2.

## TISC_2017-05-19
## tAo_2017-05-19
* Average erodibility is now calculated for the shallowest meters of bedrock.
* Erosion boundary conditions 0 and 1 are swapped.
* `-v` command line option to change velocity of Blocks based on density.
* Changed `UNIT1D` and Units structure to `BLOCK1D` and Blocks, to avoid confusion with input files (`*.UNIT`).

## TISC_2016-08-30
* An average erodibility is now calculated for the shallowest meters of bedrock.
* River Chi (scaled length) is written in column 8 of `*.bas`. A script is available to plot chi.
* GitHub used for version upload. GitHub wiki. The symbolic links under `tisc/` are now hard links to allow proper git synchronization.
* Water divides: TISC now writes the swimming distance between neighboring cells in column #14 of drainage file `*.xyw`.

## TISC_2014-06-13
* `cuthrz` accepts now hrz files with up to 120 horizons.

## TISC_2014-05-24
* Corrected bug in calculation of lake-water weight.
* Corrected bug in lake splitting during evaporation.
* Corrected bug in smoothing orographic precipitation.
## tAo_2014-09-21
* Erosion models 6 and 7 implemented.
* Graphic scripts now distinguish old sediment in cross section.
* Constant-rate scheme for erosion/sedimentation (see `Ksedim` & `Keroseol` in `*.PRM`).
* Improved sediment deformation. `deform_sed` parameter added.
* Improved sediment deformation. `deform_sed` parameter added. Sediment above moving Blocks is transported and piled in front.
* Change in units of precipitation and sed. load. Incorporation of `riverbasinwidth`.

## TISC_2014-03-25
* Erosion boundary conditions 0 and 1 are swapped.
* Erosion model #6 now smoothly transits through `eq_transport_capacity=0`.
## tAo_2012-08-29
* Change in output extensions that differ from inputs only in their case, such as `*.eet` -> `*.eeth` and `*.yse` -> `*.ysen`. This is to avoid conflict in Mac computers. Successfully compiled for Mac IOS 10.7.5.
* Solved bug during `-F` (project resume) option.

## TISC_2012-08-29
* Corrected bug in file `tisc.c`, line 108, it now reads `dy=(ymax-ymin)/(Ny-1)`.
* More river incision models added. River channel width explicitly accounted.
* Change in output extensions that differ from inputs only in case size, such as `*.eet` -> `*.eeth` and `*.PFL` -> `*.PRFL`. This is to avoid conflict in Mac computers. Successfully compiled for Mac IOS 10.7.5.
* (2013) Orographic precipitation implemented following a conservative scheme (`runoff=precipitation-evaporation`).
## tAo_2009-12-07
* Additional erosion models implemented.

## TISC_2007-07-24
* Change in version nomenclature. Debug of drainage network calculation, now very robust even with evaporation. 
* Parameters `switch_debug` and `switch_verbose` now joint into `verbose_level`.
## tAo_2006-11-02
* Removed bugs in `-Q` mode, now pressure is entered instead of thickness.

## v4.2 (Nov 2004)
* Change of name to TISC (the executable is now `tisc`).

## v4.1 (Sept 2002)
* Node edition through file `*.NDEF`. 
* Isostatic effect of the weight of lake water. Landsliding. Willgoose's fluvial incision model.
## v4.0 (2002)
* Fluvial sediment transport. 
* Loads and faults now defined in `*n.UNIT` input files (instead of `*n.CRG` file), and the file accepts now a list of parameters for the new Block. 
* Files `*.H0` and `*.W0` are now called `*.ZINI` and `*.WINI`, respectively.
* Orographic precipitation following Roe et al. 2003. 
* (Feb 2004) Glaciar ice flow and glaciar erosion and transport (Knap et al.; Braun et al., 2001).
* (2004) Orographic precipitation implemented following a conservative scheme (`runoff=precipitation-evaporation`).
* (2005) A memory bug solved in surface_transport (this was producing NaN results in a quite unpredictable manner).

## v4.0 (Jul 2000)
* Time dependence of sea level (`*.SLV`) and definition of the horizons to be marked (`*.REC`). 
* Change in load input files (`*.CRG` -> `*.LOAD` -> `*.UNIT`): Faults cut sediments; Each unit can have different erodibility. 
* Drainage including lake evaporation WORKS!!! 
* Thin sheet units with dynamic deformation. Different BCs available for erosion. 
* (Jan 2002) A constant porosity P is assumed for the sediments, so that denssedim is the total mean density: `densedim = P*1020+(1-P)*bulkdensity` (water contributes as load, but not as detritic mass). Tracking of new users.

## v3.3 to v1.0 (1994 - 2000)
* **v3.3 (Feb 2000):** Inclusion of evaporation in lakes and occurrence of underburden lakes. New output files `*.lak` and `*.bas`.
* **v3.2 (Sept 1999):** Trustful river transport including lake deposition. `*.RIV` input file read. Possible to resume stopped models (`-F`).
* **v3.1 (1996):** Improvement on surface transport accounting for lakes and river power law. Arbitrary cross sections.
* **v3.0 (1995):** Thrust loading, surface processes and sediment deformation implemented.
* **v2.0 (1995):** Moving tectonic loads implemented.
* **v1.0 (1994):** Numerical algorithm to calculate thin plate elastic and viscoelastic flexure in 3D.
## v3.4 to v1.0 (1993 - 1999)
* **v3.4 (1999):** Change in param. file keeping compatibility. Possible to resume a model that was stopped (option `-F`). Option 'c' in loads to cut Blocks with faults.
* **v3.3 (1998):** Stress history is accounted for.
* **v3.21 (1997):** Allows `sea_level` variations along time. Allows erosion/sedimentation variations along time.
* **v3.2 (1996):** Inclusion of thrust loading and sediment deformation. Command-line options.
* **v3.0 (1995):** Multilayered elastic-plastic rheology. Geoid anomaly calculations.
* **v2.0 (1994):** Inclusion of moving loads and surface processes. Gravity anomaly calculations.
* **v1.0 (1993):** Proto-model or numerical algorithm to calculate thin plate elastic and viscoelastic flexure.