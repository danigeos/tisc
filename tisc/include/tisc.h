/*
INCLUDE FILE FOR tisc.c
*/

//#include <stdio.h>
#include "geomodel.h"		/*General definitions and types for geophysical models*/

#define AUTHORSHIP		{ fprintf(stderr, "\n\t\t\t1995-2018, Daniel Garcia-Castellanos\n");}
#define	NmaxBlocks	250	/*Max number of Blocks*/
#undef	MATHLIB_SOLVER		/*Indicates whether linkage with MathLib will be done (define) or not (undef)*/

/* ===================================================================== */
/* NEW CONTEXT STRUCTS FOR REFACTORING (2D VERSION)                      */
/* ===================================================================== */

/* Holds static configuration that rarely changes during a time step */
typedef struct {
	int Nx, Ny;
	float dx, dy, dxy;
	float xmin, xmax, ymin, ymax;
	int hydro_model;
	int erosed_model;
	int verbose_level;
	float densenv, denssedim, denscrust, denswater, densice;
	float Px, Py, Pxy;
} ModelConfig;

/* Holds the active state of the simulation */
typedef struct {
	float Time, Timeini, Timefinal;
	float dt, dt_eros;
	float sea_level;
	float total_rain, total_sed_mass, total_bedrock_eros_mass;
	float total_precip_gypsum_rate, total_precip_halite_rate;
	float total_accum_gypsum, total_accum_halite;
	int numBlocks;
	int nlakes;
	float **topo;
} ModelContext;

#define GET_KREST(Krest, q_array, i, j) { /*Define the restoring force value.*/\
	    if (switch_topoest) {\
		/*If the current i,j knot is below the load then the compensation density is densinfill.*/\
		if (q_array[i][j] != 0 || Dq[i][j] != 0)   Krest = (densasthen-densinfill)*g;\
		/*Otherwise use the sediment density.*/\
		else	    Krest = (densasthen-denssedim)*g;\
	    } else    Krest = (densasthen-densenv)*g;}


extern struct GRIDNODE	*sortcell;
extern struct DRAINAGE	**drainage;
extern struct LAKE_INFO *Lake;
extern struct BLOCK 	*Blocks;

extern int erosed_model, hydro_model, mode_interp, nbasins, nlakes, n_ice_flow, n_image, n_insolation_input_points;

extern float evaporation_ct, K_ice_eros, A_ice_rheo, A_ice_slide, dt_ice, total_rain, insolation_mean;
extern float Px, Py, Pxy, CYrain, windazimut, xmin, xmax, ymin, ymax;

extern float **D, **Dq, **Dw, **eros_now, **EET, **h_water, **h_last_unit, **ice_thickness, **ice_sedm_load;
extern float **ice_velx_sl, **ice_vely_sl, **ice_velx_df, **ice_vely_df, **q, **evaporation, **precipitation;
extern float **precipitation_snow, **precipitation_file, **topo, **accumul_erosion, **Blocks_base, **w;

extern char boundary_conds[5], eros_bound_cond[5], solver_type, gif_geom[MAXLENLINE];

extern int lake_instant_fill, **lake_former_step;



extern float total_precip_gypsum_rate, total_precip_halite_rate;
extern float total_accum_gypsum, total_accum_halite;

/*functions at libreria.c*/
#include "libreria.h"

/*functions at surf_proc.c*/
int Add_Node_To_Lake (int row, int col, int i_lake);
int Add_Saddle_To_Lake (int row_sd, int col_sd, int row_tr, int col_tr, int i_lake);
int Attempt_Delete_Node_From_Lake (ModelConfig *cfg, ModelContext *ctx, int row, int col);
int Deallocate_Lake (int i_lake);
int Define_Lake (int i_lake);
int Delete_Node_From_Lake (ModelConfig *cfg, ModelContext *ctx, int row, int col);
int Diffusive_Eros (ModelConfig *cfg, ModelContext *ctx, float Kerosdif);
int Direct_mode(char *load_file_name);
int interpr_command_line_opts(int argc, char **argv); 
int Landslide_Transport (ModelConfig *cfg, ModelContext *ctx, float critical_slope);
int Divide_Lake (ModelConfig *cfg, ModelContext *ctx, int row, int col);
int Elastic_Deflection(ModelConfig *cfg, ModelContext *ctx);
int constant_rate_eros (ModelConfig *cfg, ModelContext *ctx, float Keroseol, float Ksedim, int water_load);
int Erode (ModelConfig *cfg, ModelContext *ctx, double d_mass, int row, int col);
int Surface_Transport 	(ModelConfig *cfg, ModelContext *ctx, float **topo_ant, int lake_instant_fill);
int inputs (int argc, char **argv);
int Lake_Fill (struct LAKE_INFO *Lake, ModelConfig *cfg, ModelContext *ctx, int row, int col, float hl, float dt_fv, int lake_instant_fill);
float Lake_Input_Discharge (ModelConfig *cfg, int ilake);
int Lake_Node_Number(int row, int col);
int Lake_Saddle_Number (int row, int col);
int match_parameter (char *str1, char *str2, int show, int replace, char *line);
float Minimum_Neg_Slope (ModelConfig *cfg, ModelContext *ctx, int i, int j, int *dr_row, int *dr_col);
int New_Lake ();
float Precipitation (int row, int col, int type);
int Damn_River_Node (int ia, int ja, int i,  int j);
int Rise_Damn_Node (int iia, int jja, int i, int j);
int Sediment (ModelConfig *cfg, ModelContext *ctx, double dh_sed, int row, int col, float grainsize);
int surface_processes (float **topo_ant, ModelConfig *cfg, ModelContext *ctx);
int tectload(ModelConfig *cfg, ModelContext *ctx);
int move_Blocks(ModelConfig *cfg, ModelContext *ctx);
int read_file_unit(ModelConfig *cfg, ModelContext *ctx);
int read_file_node_defs(ModelConfig *cfg, ModelContext *ctx, float dt_st);
int syntax();
int The_End(ModelConfig *cfg, ModelContext *ctx);
int Write_Ouput(ModelConfig *cfg, ModelContext *ctx);
int Unify_Lakes (ModelConfig *cfg, ModelContext *ctx, int i_lake, int i_lake_to_delete);
int Viscous_Relaxation(ModelConfig *cfg, ModelContext *ctx);

/*
extern void velocity_field_(
	double *elapsed_time, double *dt_aux, double *Lx, double *Ly, 
	int *n, int *m, int *nn, double *vissup, double *visinf, double *viscTer, double *viscosity, 
	int *nincogn, double *tallmax, double *alfa, int *nitermax, double *average_pressure, 
	double *vel_x_array, double *vel_y_array,
	char *tmpTSBCfilename, int *nbanda
);

extern void vertical_strain_rate_(double *Lx, double *Ly, int *n, int *m, int *nn,
	double *vel_x_array, double *vel_y_array, double *vert_strain_rate
);
extern void thicken_(double *dt_aux, int *n, int *m, int *nn, double *Lx, double *Ly,
	double *vel_x_array, double *vel_y_array, double *vert_strain_rate, double *layer_thickness,
	int *thicken_BC
);
*/

#include "tao+tisc.h"
