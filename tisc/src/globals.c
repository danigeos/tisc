#include "tisc.h"
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

/* Definitions for tisc-specific global variables */
struct GRIDNODE *sortcell;
struct DRAINAGE **drainage;
struct LAKE_INFO *Lake;
struct BLOCK *Blocks;
int erosed_model, hydro_model, mode_interp, nbasins, nlakes, n_ice_flow, n_image, n_insolation_input_points;
float evaporation_ct, K_ice_eros, A_ice_rheo, A_ice_slide, dt_ice, total_rain, insolation_mean;
float Px, Py, Pxy, CYrain, windazimut, xmin, xmax, ymin, ymax;
float **D, **Dq, **Dw, **eros_now, **EET, **h_water, **h_last_unit, **ice_thickness, **ice_sedm_load;
float **ice_velx_sl, **ice_vely_sl, **ice_velx_df, **ice_vely_df, **q, **evaporation, **precipitation;
float **precipitation_snow, **precipitation_file, **topo, **accumul_erosion, **Blocks_base, **w;
char boundary_conds[5], eros_bound_cond[5], solver_type, gif_geom[MAXLENLINE];
int lake_instant_fill = 0, **lake_former_step;

bool switch_debug = false;

ParameterEntry param_table[] = {
    {"Nx", PARAM_TYPE_INT, &Nx, 0, false},
    {"Ny", PARAM_TYPE_INT, &Ny, 0, false},
    {"xmin", PARAM_TYPE_FLOAT, &xmin, 0, false},
    {"xmax", PARAM_TYPE_FLOAT, &xmax, 0, false},
    {"ymin", PARAM_TYPE_FLOAT, &ymin, 0, false},
    {"ymax", PARAM_TYPE_FLOAT, &ymax, 0, false},
    {"Te", PARAM_TYPE_FLOAT, &Te_default, 0, false},
    {"zini", PARAM_TYPE_FLOAT, &zini, 0, false},
    {"random_topo", PARAM_TYPE_FLOAT, &random_topo, 0, false},
    {"mode_interp", PARAM_TYPE_INT, &mode_interp, 0, false},
    {"densasthen", PARAM_TYPE_FLOAT, &densasthen, 0, false},
    {"densmantle", PARAM_TYPE_FLOAT, &densmantle, 0, false},
    {"denscrust", PARAM_TYPE_FLOAT, &denscrust, 0, false},
    {"densinfill", PARAM_TYPE_FLOAT, &densinfill, 0, false},
    {"denssedim", PARAM_TYPE_FLOAT, &denssedim, 0, false},
    {"densenv", PARAM_TYPE_FLOAT, &densenv, 0, false},
    {"sed_porosity", PARAM_TYPE_FLOAT, &sed_porosity, 0, false},
    {"compact_depth", PARAM_TYPE_FLOAT, &compact_depth, 0, false},
    {"boundary_conds", PARAM_TYPE_STRING, &boundary_conds, sizeof(boundary_conds), false},
    {"Px", PARAM_TYPE_FLOAT, &Px, 0, false},
    {"Py", PARAM_TYPE_FLOAT, &Py, 0, false},
    {"Pxy", PARAM_TYPE_FLOAT, &Pxy, 0, false},
    {"hydro_model", PARAM_TYPE_INT, &hydro_model, 0, false},
    {"rain", PARAM_TYPE_FLOAT, &rain, 0, false},
    {"Krain", PARAM_TYPE_FLOAT, &Krain, 0, false},
    {"relhumid", PARAM_TYPE_FLOAT, &relative_humidity, 0, false},
    {"windazimut", PARAM_TYPE_FLOAT, &windazimut, 0, false},
    {"CXrain", PARAM_TYPE_FLOAT, &CXrain, 0, false},
    {"CYrain", PARAM_TYPE_FLOAT, &CYrain, 0, false},
    {"rain_per", PARAM_TYPE_FLOAT, &rain_per, 0, false},
    {"rain_amp", PARAM_TYPE_FLOAT, &rain_amp, 0, false},
    {"evaporation", PARAM_TYPE_FLOAT, &evaporation_ct, 0, false},
    {"lost_rate", PARAM_TYPE_FLOAT, &lost_rate, 0, false},
    {"permeability", PARAM_TYPE_FLOAT, &permeability, 0, false},
    {"erosed_model", PARAM_TYPE_INT, &erosed_model, 0, false},
    {"Kerosdif", PARAM_TYPE_FLOAT, &Kerosdif, 0, false},
    {"Keroseol", PARAM_TYPE_FLOAT, &Keroseol, 0, false},
    {"Ksedim", PARAM_TYPE_FLOAT, &Ksedim, 0, false},
    {"critical_slope", PARAM_TYPE_FLOAT, &critical_slope, 0, false},
    {"K_river_cap", PARAM_TYPE_FLOAT, &K_river_cap, 0, false},
    {"erodibility", PARAM_TYPE_FLOAT, &erodibility, 0, false},
    {"erodibility_sed", PARAM_TYPE_FLOAT, &erodibility_sed, 0, false},
    {"critical_stress", PARAM_TYPE_FLOAT, &critical_stress, 0, false},
    {"grain_size", PARAM_TYPE_FLOAT, &initial_grain_size, 0, false},
    {"grain_size_decay", PARAM_TYPE_FLOAT, &distance_half_grainsize, 0, false},
    {"l_fluv_sedim", PARAM_TYPE_FLOAT, &l_fluv_sedim, 0, false},
    {"temp_sea_level", PARAM_TYPE_FLOAT, &temp_sea_level, 0, false},
    {"deform_sed", PARAM_TYPE_BOOL, &deform_sed, 0, false},
    {"K_ice_eros", PARAM_TYPE_FLOAT, &K_ice_eros, 0, false},
    {"dt_ice", PARAM_TYPE_FLOAT, &dt_ice, 0, false},
    {"n_ice_flow", PARAM_TYPE_INT, &n_ice_flow, 0, false},
    {"A_ice_rheo", PARAM_TYPE_FLOAT, &A_ice_rheo, 0, false},
    {"A_ice_slide", PARAM_TYPE_FLOAT, &A_ice_slide, 0, false},
    {"eros_bound_cond", PARAM_TYPE_STRING, &eros_bound_cond, sizeof(eros_bound_cond), false},
    {"Timeini", PARAM_TYPE_FLOAT, &Timeini, 0, false},
    {"Timefinal", PARAM_TYPE_FLOAT, &Timefinal, 0, false},
    {"tau", PARAM_TYPE_FLOAT, &tau, 0, false},
    {"dt", PARAM_TYPE_FLOAT, &dt, 0, false},
    {"dt_eros", PARAM_TYPE_FLOAT, &dt_eros, 0, false},
    {"dt_record", PARAM_TYPE_FLOAT, &dt_record, 0, false},
    {"isost_model", PARAM_TYPE_INT, &isost_model, 0, false},
    {"water_load", PARAM_TYPE_INT, &water_load, 0, false},
    {"switch_topoest", PARAM_TYPE_BOOL, &switch_topoest, 0, false},
    {"switch_files", PARAM_TYPE_BOOL, &switch_write_file, 0, false},
    {"switch_ps", PARAM_TYPE_INT, &switch_ps, 0, false},
    {"verbose_level", PARAM_TYPE_INT, &verbose_level, 0, false},

    // Salt and Ions
    {"C_Ca_SEA", PARAM_TYPE_FLOAT, &C_Ca_SEA, 0, false},
    {"C_SO4_SEA", PARAM_TYPE_FLOAT, &C_SO4_SEA, 0, false},
    {"C_Na_SEA", PARAM_TYPE_FLOAT, &C_Na_SEA, 0, false},
    {"C_Cl_SEA", PARAM_TYPE_FLOAT, &C_Cl_SEA, 0, false},
    {"C_Ca_RIV", PARAM_TYPE_FLOAT, &C_Ca_RIV, 0, false},
    {"C_SO4_RIV", PARAM_TYPE_FLOAT, &C_SO4_RIV, 0, false},
    {"C_Na_RIV", PARAM_TYPE_FLOAT, &C_Na_RIV, 0, false},
    {"C_Cl_RIV", PARAM_TYPE_FLOAT, &C_Cl_RIV, 0, false},
    {"GYPSUM_PRECIP_CN", PARAM_TYPE_FLOAT, &GYPSUM_PRECIP_CN, 0, false},
    {"HALITE_PRECIP_CN", PARAM_TYPE_FLOAT, &HALITE_PRECIP_CN, 0, false},

    // Old versions
    {"erodability", PARAM_TYPE_FLOAT, &erodibility, 0, true},
    {"erodability_sed", PARAM_TYPE_FLOAT, &erodibility_sed, 0, true},
    {"switch_verbose", PARAM_TYPE_INT, &verbose_level, 0, true},
    {"switch_debug", PARAM_TYPE_BOOL, &switch_debug, 0, true},
    {"alt0", PARAM_TYPE_FLOAT, &zini, 0, true},
    {"lith_type", PARAM_TYPE_INT, &isost_model, 0, true},
    {"erosed_type", PARAM_TYPE_INT, &erosed_model, 0, true},
    {"switch_hydro", PARAM_TYPE_INT, &hydro_model, 0, true},
    {"leng_fluv_eros", PARAM_TYPE_FLOAT, &erodibility, 0, true},
    {"leng_fluv_sedim", PARAM_TYPE_FLOAT, &l_fluv_sedim, 0, true},
    {"switch_erosed", PARAM_TYPE_INT, &erosed_model, 0, true},
    {"switch_sea", PARAM_TYPE_INT, &water_load, 0, true},
    {"l_fluv_eros", PARAM_TYPE_FLOAT, &erodibility, 0, true},
    {"l_fluv_eros_sed", PARAM_TYPE_FLOAT, &erodibility_sed, 0, true},
    {"dtmemounit", PARAM_TYPE_FLOAT, &dt_record, 0, true}
};

const int num_params = sizeof(param_table) / sizeof(ParameterEntry);