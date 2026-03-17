#include "tao.h"
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
float spl_m, spl_n, l_fluv_sedim, lost_rate, permeability, rain, Krain, relative_humidity, CXrain;
float rain_per, rain_amp, total_bedrock_eros_mass, total_sed_mass;
float zini, dt_record, sed_porosity, compact_depth, last_time_file_time, random_topo = 0;
float **var_sea_level, **var_eros_level, **var_insolation, *horiz_record_time;
bool switch_file_out = false, switch_gradual, switch_topoest, switch_write_file_Blocks, deform_sed;

/* Global variable definitions for tao */
int *sortcell;
struct DRAINAGE_1D  *drainage;
struct LAKE_INFO_1D *Lake;
struct BLOCK_1D	*Blocks ;

int 	imomentmax, nmax_input_points, nx_temp_input;
char	eros_bound_cond[2];
int	erosed_model, hydro_model, mode_interp, nbasins, nlakes, boundary_conds, n_image;

float 	evaporation_ct, riverbasinwidth, total_evap_water, total_lost_water, total_rain;
float   x0, xf, xmin, xmax, zmin, zmax, horz_force, vert_force, appmoment;

float	*crust_thick, *D, *Dq, *Dw, *h_water, *h_last_unit, *q, *Te;
float   *eros_now, *precipitation, *evaporation, *total_erosion, *topo;
float   *Blocks_base, *upper_crust_thick, *w;
float   **Temperature, **stress, **yieldcompres, **yieldextens;

char 	gif_geom[MAXLENLINE];
bool	switch_insert_load=false, switch_strs_history=false, switch_YSE_file=false;
bool switch_debug = false;

// Global parameter table for match_parameter function
ParameterEntry param_table[] = {
    {"Nx", PARAM_TYPE_INT, &Nx, 0, false},
    {"Nz", PARAM_TYPE_INT, &Nz, 0, false},
    {"x0", PARAM_TYPE_FLOAT, &x0, 0, false},
    {"xf", PARAM_TYPE_FLOAT, &xf, 0, false},
    {"xmin", PARAM_TYPE_FLOAT, &xmin, 0, false},
    {"xmax", PARAM_TYPE_FLOAT, &xmax, 0, false},
    {"zmin", PARAM_TYPE_FLOAT, &zmin, 0, false},
    {"zmax", PARAM_TYPE_FLOAT, &zmax, 0, false},
    {"Te", PARAM_TYPE_FLOAT, &Te_default, 0, false},
    {"crust", PARAM_TYPE_FLOAT, &crust_thick_default, 0, false},
    {"ucrust", PARAM_TYPE_FLOAT, &upper_crust_thick_default, 0, false},
    {"zini", PARAM_TYPE_FLOAT, &zini, 0, false},
    {"random_topo", PARAM_TYPE_FLOAT, &random_topo, 0, false},
    {"densasthen", PARAM_TYPE_FLOAT, &densasthen, 0, false},
    {"densmantle", PARAM_TYPE_FLOAT, &densmantle, 0, false},
    {"denscrust", PARAM_TYPE_FLOAT, &denscrust, 0, false},
    {"densinfill", PARAM_TYPE_FLOAT, &densinfill, 0, false},
    {"denssedim", PARAM_TYPE_FLOAT, &denssedim, 0, false},
    {"densenv", PARAM_TYPE_FLOAT, &densenv, 0, false},
    {"sed_porosity", PARAM_TYPE_FLOAT, &sed_porosity, 0, false},
    {"compact_depth", PARAM_TYPE_FLOAT, &compact_depth, 0, false},
    {"erosed_model", PARAM_TYPE_INT, &erosed_model, 0, false},
    {"Kerosdif", PARAM_TYPE_FLOAT, &Kerosdif, 0, false},
    {"Keroseol", PARAM_TYPE_FLOAT, &Keroseol, 0, false},
    {"Ksedim", PARAM_TYPE_FLOAT, &Ksedim, 0, false},
    {"hydro_model", PARAM_TYPE_INT, &hydro_model, 0, false},
    {"rain", PARAM_TYPE_FLOAT, &rain, 0, false},
    {"Krain", PARAM_TYPE_FLOAT, &Krain, 0, false},
    {"relhumid", PARAM_TYPE_FLOAT, &relative_humidity, 0, false},
    {"CXrain", PARAM_TYPE_FLOAT, &CXrain, 0, false},
    {"evaporation", PARAM_TYPE_FLOAT, &evaporation_ct, 0, false},
    {"riverbasinwidth", PARAM_TYPE_FLOAT, &riverbasinwidth, 0, false},
    {"lost_rate", PARAM_TYPE_FLOAT, &lost_rate, 0, false},
    {"K_river_cap", PARAM_TYPE_FLOAT, &K_river_cap, 0, false},
    {"erodibility", PARAM_TYPE_FLOAT, &erodibility, 0, false},
    {"erodibility_sed", PARAM_TYPE_FLOAT, &erodibility_sed, 0, false},
    {"l_fluv_sedim", PARAM_TYPE_FLOAT, &l_fluv_sedim, 0, false},
    {"temp_sea_level", PARAM_TYPE_FLOAT, &temp_sea_level, 0, false},
    {"deform_sed", PARAM_TYPE_BOOL, &deform_sed, 0, false},
    {"eros_bound_cond", PARAM_TYPE_STRING, &eros_bound_cond, sizeof(eros_bound_cond), false},
    {"Timeini", PARAM_TYPE_FLOAT, &Timeini, 0, false},
    {"Timefinal", PARAM_TYPE_FLOAT, &Timefinal, 0, false},
    {"tau", PARAM_TYPE_FLOAT, &tau, 0, false},
    {"dt", PARAM_TYPE_FLOAT, &dt, 0, false},
    {"dt_eros", PARAM_TYPE_FLOAT, &dt_eros, 0, false},
    {"dt_record", PARAM_TYPE_FLOAT, &dt_record, 0, false},
    {"boundary_conds", PARAM_TYPE_INT, &boundary_conds, 0, false},
    {"horz_force", PARAM_TYPE_FLOAT, &horz_force, 0, false},
    {"vert_force", PARAM_TYPE_FLOAT, &vert_force, 0, false},
    {"moment", PARAM_TYPE_FLOAT, &appmoment, 0, false},
    {"isost_model", PARAM_TYPE_INT, &isost_model, 0, false},
    {"water_load", PARAM_TYPE_BOOL, &water_load, 0, false},
    {"switch_topoest", PARAM_TYPE_BOOL, &switch_topoest, 0, false},
    {"grav_anom", PARAM_TYPE_INT, &grav_anom_type, 0, false},
    {"switch_files", PARAM_TYPE_BOOL, &switch_write_file, 0, false},
    {"switch_ps", PARAM_TYPE_INT, &switch_ps, 0, false},
    {"verbose_level", PARAM_TYPE_INT, &verbose_level, 0, false},

    // Old versions (for backward compatibility)
    {"erodability", PARAM_TYPE_FLOAT, &erodibility, 0, true},
    {"erodability_sed", PARAM_TYPE_FLOAT, &erodibility_sed, 0, true},
    {"switch_verbose", PARAM_TYPE_INT, &verbose_level, 0, true}, // Note: verbose_level is int, not bool
    {"switch_debug", PARAM_TYPE_BOOL, &switch_debug, 0, true},
    {"alt0", PARAM_TYPE_FLOAT, &zini, 0, true},
    {"lith_type", PARAM_TYPE_INT, &isost_model, 0, true},
    {"erosed_type", PARAM_TYPE_INT, &erosed_model, 0, true},
    {"switch_erosed", PARAM_TYPE_INT, &erosed_model, 0, true},
    {"switch_sea", PARAM_TYPE_BOOL, &water_load, 0, true},
    {"ymin", PARAM_TYPE_FLOAT, &zmin, 0, true}, // Assuming zmin/ymax are equivalent to ymin/ymax in 1D
    {"ymax", PARAM_TYPE_FLOAT, &zmax, 0, true}, // Assuming zmin/ymax are equivalent to ymin/ymax in 1D
    {"dtmemounit", PARAM_TYPE_FLOAT, &dt_record, 0, true}
};

const int num_params = sizeof(param_table) / sizeof(ParameterEntry);