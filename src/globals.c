#include "types_n_defs.h"
#include "universal.h"
#include "geomodel.h"
#include "tao+tisc.h"

/* From universal.h */
float g = 9.81;
float pi = 3.14159265;
float sqrt2 = 1.4142136;
float CGU = 6.6732E-11;
float number_e = 2.7182818;
float densice = 920;
float denswater = 1020;
float E = 7E+10;
float nu = 0.25;
float Rearth = 6378e3;
float viscwater = 8.9e-4;

int Nx, Ny, Nz;
int switch_ps, verbose_level;
float dx, dy, dz, dxy;
char version[LENGTHVERS], version_input[LENGTHVERS];
bool switch_geograph_coor, switch_write_file;

/* From geomodel.h */
int grav_anom_type, isost_model, water_load;
float Te_default, crust_thick_default, upper_crust_thick_default;
float densasthen, densmantle, denscrust, densinfill, denssedim, densenv;
float sea_level = 0, temp_sea_level;
float Time = 0, Timeini = 0, Timefinal = 0;
float dt = 0, dt_eros, tau;
char projectname[MAXLONFICH], title[MAXLONLINE];

/* From tao+tisc.h */
int nloads, n_sea_level_input_points, n_eros_level_input_points, n_record_times;
int i_first_Block_load, i_Block_insert, numBlocks, nwrotenfiles, run_type;
float Kerosdif, Keroseol, Ksedim, critical_slope, K_river_cap, erodibility, erodibility_sed, critical_stress;
float initial_grain_size = 1.0;
float distance_half_grainsize = 5000.0;

/* Salt & Ion parameters */
float C_Ca_SEA = 0.048;
float C_SO4_SEA = 2.04;
float C_Na_SEA = 10.70;
float C_Cl_SEA = 19.30;
float C_Ca_RIV = 0.06;
float C_SO4_RIV = 0.01;
float C_Na_RIV = 0.05;
float C_Cl_RIV = 0.05;
float GYPSUM_PRECIP_CN = 5.25;
float HALITE_PRECIP_CN = 272.1;

float spl_m, spl_n, l_fluv_sedim, lost_rate, permeability, rain, Krain, relative_humidity, CXrain;
float rain_per, rain_amp, total_bedrock_eros_mass, total_sed_mass;
float total_precip_gypsum_rate = 0, total_precip_halite_rate = 0;
float total_accum_gypsum = 0, total_accum_halite = 0;
float zini, dt_record, sed_porosity, compact_depth, last_time_file_time, random_topo = 0;
float **var_sea_level, **var_eros_level, **var_insolation, *horiz_record_time;
bool switch_file_out = false, switch_gradual, switch_topoest, switch_write_file_Blocks, deform_sed;