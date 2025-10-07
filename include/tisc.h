/*
INCLUDE FILE FOR tisc.c
*/

//#include <stdio.h>
#include "geomodel.h"		/*General definitions and types for geophysical models*/

#define AUTHORSHIP		{ fprintf(stderr, "\n\t\t\t1995-2018, Daniel Garcia-Castellanos\n");}
#define	NmaxBlocks	250	/*Max number of Blocks*/
#undef	MATHLIB_SOLVER		/*Indicates whether linkage with MathLib will be done (define) or not (undef)*/

#define GET_KREST(Krest, q, i, j) { /*Define the restoring force value.*/\
	    if (switch_topoest) {\
		/*If the current i,j knot is below the load then the compensation density is densinfill.*/\
		if (q[i][j])   Krest = (densasthen-densinfill)*g;\
		/*Otherwise use the sediment density.*/\
		else	    Krest = (densasthen-denssedim)*g;\
	    } else    Krest = (densasthen-densenv)*g;}


struct GRIDNODE	*sortcell;
struct DRAINAGE	**drainage;
struct LAKE_INFO *Lake;		/*Lake[0] does not exist; Lake[1] is the sea or the first normal lake.*/
struct BLOCK 	*Blocks;

int	erosed_model, 
	hydro_model, 
	mode_interp, 		
	nbasins, 
	nlakes, 		/*number of lakes >= 0 */
	n_ice_flow, 
	n_image,
	n_insolation_input_points;

float 	evaporation_ct, 		/*[m3/s/m2].*/
	K_ice_eros, 
	A_ice_rheo, 
	A_ice_slide, 
	dt_ice, 
	total_rain, 
	insolation_mean, 
	Px, Py, Pxy, 		/*Horizontal external load */
	CYrain, windazimut, 	/*[m], [degrees]*/
	xmin, xmax, ymin, ymax;	/*Model domain*/

float 	**D,			/*Equivalent rigidity grid*/
	**Dq, 
	**Dw, 
	**eros_now,  		/*Erosion during present dt at each node*/
	**EET,			/*Equivalent Elastic Thickness grid*/
	**h_water, 		/*Sea/lake water column thickness computed if water_load==1*/
	**h_last_unit,
	**ice_thickness, 
	**ice_sedm_load, 
	**ice_velx_sl, **ice_vely_sl, 
	**ice_velx_df, **ice_vely_df, 
	**q, 
	**evaporation,  	/*Evaporation distribution*/
	**precipitation,  	/*Runoff distribution (only liquid precipitation)*/
	**precipitation_snow,  	/*Snowfall distribution*/
	**precipitation_file,  	/*Runoff distribution defined in file*/
	**topo,  		/*Topography relative to z=0*/
	**accumul_erosion,  	/*Total erosion at each node*/
	**Blocks_base,		/*Base of Blocks[0] (the lowest) measured from the current position of the original datum surface (now deflected by isostasy*/
	**w;

char	boundary_conds[5],	/*Boundary conditions at N, S, E, W.*/
	eros_bound_cond[5],
	solver_type,		/*'l' for local solver; 'g' for genspd; 'm' for mathlib*/
	gif_geom[MAXLENLINE];

int	lake_instant_fill=0,
	**lake_former_step;	/**/



/*functions at libreria.c*/
#include "../lib/libreria.h"

/*functions at surf_proc.c*/
int Add_Node_To_Lake (int row, int col, int i_lake);
int Add_Saddle_To_Lake (int row_sd, int col_sd, int row_tr, int col_tr, int i_lake);
int Attempt_Delete_Node_From_Lake (int row, int col);
float calculate_topo(float **topo);
int Deallocate_Lake (int i_lake);
int Define_Lake (int i_lake);
int Delete_Node_From_Lake (int row, int col);
int Diffusive_Eros (float Kerosdif, float dt, float dt_eros);
int Direct_mode(char *load_file_name);
int interpr_command_line_opts(int argc, char **argv); 
int Landslide_Transport (float critical_slope, float dt, float dt_eros);
int Divide_Lake (int row, int col);
int Elastic_Deflection();
int constant_rate_eros (float **topo, float Keroseol, float Ksedim, float sea_level, BOOL switch_sea, float dt, float Time);
int Erode (double d_mass, int row, int col);
int Surface_Transport 	(float **topo, float **topo_ant, float dt, float dt_eros, BOOL switch_erosed, int lake_instant_fill);
int inputs (int argc, char **argv);
int Lake_Fill (struct LAKE_INFO *Lake, int row, int col, float hl, float dt_fv, int lake_instant_fill);
float Lake_Input_Discharge (int ilake);
int Lake_Node_Number(int row, int col);
int Lake_Saddle_Number (int row, int col);
int match_parameter (char *str1, char *str2, int show, int replace, char *line);
float Minimum_Neg_Slope (int i, int j, int *dr_row, int *dr_col);
int New_Lake ();
float Precipitation (int row, int col, int type);
int Damn_River_Node (int ia, int ja, int i,  int j);
int Rise_Damn_Node (int iia, int jja, int i, int j);
int Sediment (double dh_sed, int row, int col);
int surface_processes (float **topo_ant);
int tectload();
int move_Blocks();
int read_file_unit();
int syntax();
int The_End();
int Write_Ouput();
int Unify_Lakes (int i_lake, int i_lake_to_delete);
int Viscous_Relaxation();

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
