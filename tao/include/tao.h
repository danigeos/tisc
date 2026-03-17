/*
INCLUDE FILE FOR tao.c
*/

#include <stdbool.h>

#include "geomodel.h"			/*General definitions and types for geophysical models*/

/* ===================================================================== */
/* NEW CONTEXT STRUCTS FOR REFACTORING (1D VERSION)                      */
/* ===================================================================== */

typedef struct {
	int Nx, Nz;
	float dx, dz;
	float xmin, xmax, x0, xf;
	int hydro_model;
	int erosed_model;
	int verbose_level;
	float densenv, denssedim, denscrust, denswater, densasthen, densmantle, densinfill;
	float riverbasinwidth, sed_porosity;
} ModelConfig;

typedef struct {
	float Time, Timeini, Timefinal;
	float dt, dt_eros;
	float sea_level;
	float total_rain, total_sed_mass, total_bedrock_eros_mass;
	int numBlocks;
	int nlakes;
	float *topo;
} ModelContext;

#define	MAXETERR	10e6		/*Maximum error area allowed as Te convergence criteria in the stress distrib. algorithm (m2)*/
#define	MAX_Te_LOC_VAR	200		/*Maximum error in Te at any point. Also a convergence criteria in the stress distrib. algorithm (m2)*/
#define	NMAXRHEOITERS 	100		/*Maximum number of iterations to find Te convergence in the stress distrib. algorithm*/
#define	NmaxBlocks	400			/*Maximum number of Blocks to be recorded*/
#define Read_Param(x, y)  {char line[MAXLENLINE+200]; while ((sscanf(lineptr=fgets(line, MAXLENLINE+200, file), x, y )) < 1) \
	if (lineptr==NULL) {fprintf(stderr, "\nERROR reading parameters file: parameter not found.\n");	break;}}

#define PRINT_ARRAY_INFO_1D(array, name, units, unitsintegr) {\
	    int i,imin=SIGNAL,imax=SIGNAL; float max=-1e19, min=1e19, vol=0;\
	    for (i=(xmin-x0)/dx; i<Nx-(xf-xmax)/dx; i++) {\
	    	    vol += array[i]*dx;\
	    	    if (max<array[i])  {max=array[i]; imax=i;}\
	    	    if (min>array[i])  {min=array[i]; imin=i;}}\
	    PRINT_SUMLINE("%s:\tmax= %.2e %s @x= %.1f km\tmin= %.2e %s @x= %.1f km\tintegr= %.2e %s", name, max, units, (x0+imax*dx)/1e3, min, units, (x0+imin*dx)/1e3, vol, unitsintegr);}

/*convert between sediment thickness and sediment grain mass*/
#define MASS2SEDTHICK_1D(cfg, mass)	((mass) /(cfg->denssedim-cfg->sed_porosity*cfg->denswater)/cfg->dx/cfg->riverbasinwidth)	/*converts sediment mass into sediment thickness*/
#define THICK2SEDMASS_1D(cfg, thick)	((thick)*(cfg->denssedim-cfg->sed_porosity*cfg->denswater)*cfg->dx*cfg->riverbasinwidth)	/*converts sediment thickness into sediment mass*/
#define GET_KREST(Krest, q, i) { /*Define the restoring force value.*/\
	    if (switch_topoest) {\
		/*If the current i knot is below the load then the compensation density is densinfill.*/\
		if (q[i])   Krest = (densasthen-densinfill)*g;\
		/*Otherwise use the sediment density.*/\
		else	    Krest = (densasthen-denssedim)*g;\
	    } else    Krest = (densasthen-densenv)*g;}

extern int *sortcell;
extern struct DRAINAGE_1D  *drainage;
extern struct LAKE_INFO_1D *Lake;	/*Lake[0] does not exist; Lake[1] is the sea or the first normal lake.*/
extern struct BLOCK_1D	*Blocks ; 	/*Blocks array: Blocks[numBlocks] first Block is 0, last is numBlocks-1;*/

extern int 	imomentmax, 
	nmax_input_points,	/*Maximum number of input points in moving loads*/
	nx_temp_input;		/*Number of x points with defined geotherm in .TMP file*/

extern char	eros_bound_cond[2];
extern int	erosed_model, 
	hydro_model, 
	mode_interp, 		
	nbasins, 
	nlakes, 		/*number of lakes >= 0 */
	boundary_conds,
	n_image;

extern float 	evaporation_ct, 		/*[m3/s/m2].*/
	riverbasinwidth, 		/*[m]*/
	total_evap_water,
	total_lost_water,
	total_rain, 
	x0, xf, 			/*X horizontal direction logical limits of the model [m]*/
	xmin, xmax, zmin, zmax, 	/*Screen limits for the graphic output [m]*/
	horz_force, vert_force, appmoment; 	/*Horizontal tectonic force and Vertical boundary applied force (in N/m); Boundary applied moment (in N m/m)*/


/*Arrays:*/
extern float	*crust_thick,			/*Crust thickness array*/
	*D, 				/*Rigidity		[N m] */
	*Dq, 				/*Load increment	[N/m/m] */
	*Dw, 				/*Deflection increment	[m] */
	*h_water, 				/*Water load thickness	[m] */
	*h_last_unit,		 	/*Height at last unit file [m] */
	*q, 				/*load 			[N/m/m] */
	*Te,				/*E.E.T. array*/
	*eros_now,
	*precipitation,
	*evaporation,
	*total_erosion, 		/*Cumulated erosion at each point*/
	*topo,  			/*Topography over sea level*/
	*Blocks_base,			/*Base of Blocks[0] (the lowest) measured from the current position of the original datum surface (now deflected by isostasy).*/
	*upper_crust_thick,		/*Upper crust thickness array*/
	*w, 				/*Deflection		[m] */
	**Temperature, 			/*[�C]*/
	**stress, 			/*[Pa]*/
	**yieldcompres,			/*yield stress envelope -- compresive  [Pa]*/
	**yieldextens;			/*yield stress envelope -- extensive  [Pa]*/


extern char 	gif_geom[MAXLENLINE];


/*Boolean type variables:*/
extern bool	switch_insert_load,		/*Si to insert the load height beneath the first load Block*/
	switch_strs_history, 
	switch_YSE_file;			/*SI if yield stress is given directly from a .YSE file*/
extern bool switch_debug;




/*Function prototypes in libreria.c:*/
int yield_stress_envelope (
    ModelConfig *cfg,
    float *Temperatura,     /*Temperature array*/
    float z0, 
    float Uppercrustbase,
    float Lowercrustbase,
    float *yieldcompres,
    float *yieldextens,
    float *espmecanptr 
);

/*Function prototypes in libreria.c:*/
float **alloc_matrix  	(int num_fil, int num_col);
double **alloc_matrix_dbl 	(int num_fil, int num_col);
float geoidanompolig 		(float *, float *, int, float, float, float);
float gravanompolig		(float *, float *, int, float, float, float);
int	readinterplin		(FILE *, float *, int, float, float);
float interpol_in_xy_data 	(float *x_array, float *y_array, int n_x, float x);

/*Declaration of functions in taosp.c*/
int Surface_Transport (ModelConfig *cfg, ModelContext *ctx, float *topo);
int Fluvial_Transport (ModelConfig *cfg, ModelContext *ctx, float *topo, float dt_fv);
int Calculate_Discharge (int *sortcell, ModelConfig *cfg, ModelContext *ctx);
int Define_Drainage_Net (int *sortcell, ModelConfig *cfg, ModelContext *ctx);
int Lake_Fill (struct LAKE_INFO_1D *Lake, ModelConfig *cfg, ModelContext *ctx, int ind, float hl, float dt_fv);
int Calculate_Precipitation_Evaporation (ModelConfig *cfg, ModelContext *ctx);
int Add_Node_To_Lake (int row, int i_lake);
int Add_Saddle_To_Lake (int row_sd, int row_tr, int i_lake);
int Attempt_Delete_Node_From_Lake (ModelConfig *cfg, ModelContext *ctx, int row);
int constant_rate_eros (ModelConfig *cfg, ModelContext *ctx, float Keroseol, float Ksedim, int water_load, int n_eros_level_input_points, float **var_eros_level, float *eros_level);
int Deallocate_Lake (ModelContext *ctx, int i_lake);
int Define_Lake (int i_lake);
int Delete_Node_From_Lake (ModelConfig *cfg, ModelContext *ctx, int ln);
int Diffusive_Eros_1D (ModelConfig *cfg, ModelContext *ctx, float Kerosdif); 
int Divide_Lake (ModelConfig *cfg, ModelContext *ctx, int ind);
int Erode (ModelConfig *cfg, ModelContext *ctx, double d_mass, int ind);
float Lake_Input_Discharge (ModelConfig *cfg, int ilake);
int	Lake_Node_Number(int row);
int	Lake_Saddle_Number (int row);
int move_Blocks(ModelConfig *cfg, ModelContext *ctx);
int	New_Lake (ModelConfig *cfg, ModelContext *ctx);
int	Repare_Blocks(ModelConfig *cfg, ModelContext *ctx);
int	Sediment (ModelConfig *cfg, ModelContext *ctx, double d_mass, int ind);
int	tectload(ModelConfig *cfg, ModelContext *ctx);
int	Unify_Lakes (ModelConfig *cfg, ModelContext *ctx, int i_lake, int i_lake_to_delete);
int	Landslide_Transport (ModelConfig *cfg, ModelContext *ctx, float critical_slope);
int flexural_stats (ModelConfig *cfg, ModelContext *ctx, float *moment);
int Rheo_Flex_Iter (ModelConfig *cfg, ModelContext *ctx);
int solveLES (double **A, double *b, int Nx, int nds, int ndi, float *x);
float moment_calculator (ModelConfig *cfg, float 	d2wdx2, 
			float 	*yieldcompres, 
			float 	*yieldextens, 
			float 	*stress, 
			float 	decoupl_depth, 		/*In m. Only used when isost_model==4*/
			float 	*refstressdir, 
			int 	*ncapas); 		/*Number of decoupled layers*/
float moment_calculator_hist (ModelConfig *cfg,
			float 	d2wdx2, 
			float 	*yieldcompres, 
			float 	*yieldextens, 
			float 	*stress, 
			float 	decoupl_depth, 		/*In m. Only used when isost_model==4*/
			float 	*totalmoment,  		/*Total (cumulative) moment at this point*/
			int 	*nlayers); 		/*Number of decoupled layers*/

int write_file_Temperature_initial (ModelConfig *cfg, ModelContext *ctx);
int inputs(ModelConfig *cfg, ModelContext *ctx, int argc, char **argv);
int Elastoplastic_Deflection(ModelConfig *cfg, ModelContext *ctx);
int Viscous_Relaxation(ModelConfig *cfg, ModelContext *ctx);
int surface_processes(ModelConfig *cfg, ModelContext *ctx);
int Write_Ouput(ModelConfig *cfg, ModelContext *ctx);
int syntax();
int The_End(ModelConfig *cfg, ModelContext *ctx);
int interpr_command_line_opts(int argc, char **argv);
int Direct_mode(ModelConfig *cfg, ModelContext *ctx, char *load_file_name);
int read_file_unit(ModelConfig *cfg, ModelContext *ctx);
int read_file_parameters(int show, int reformat);
int read_file_resume(char *filename);
int read_file_sea_level();
int read_file_horiz_record_time();
int read_file_Te();
int read_file_Crust_Thick(float crust_thick_default);
int read_file_Upper_Crust_Thick(float crust_thick_supdefault);
int read_file_initial_topo(float *h);
int read_file_initial_deflection(float *w);
int read_file_Temperature(ModelConfig *cfg, ModelContext *ctx);
int match_parameter(char *str1, char *str2, int show, int replace, char *line);
int write_file_erosed (ModelConfig *cfg, ModelContext *ctx, float *total_erosion);
int write_file_grav_anom (ModelConfig *cfg, ModelContext *ctx, float *gravanom, float *geoidanom);
int write_file_maxmompoint (ModelConfig *cfg, ModelContext *ctx);
int write_file_stress (ModelConfig *cfg, ModelContext *ctx);
int write_file_Te (ModelConfig *cfg, ModelContext *ctx);
int write_file_Temperature (ModelConfig *cfg, ModelContext *ctx);
int write_file_time (ModelConfig *cfg, ModelContext *ctx, float *w, float *topo);
int write_file_Blocks(ModelConfig *cfg, ModelContext *ctx);
int write_file_resume(ModelConfig *cfg, ModelContext *ctx);
int Delete_Block(ModelContext *ctx, int i_Block);
int gradual_Block(ModelConfig *cfg, ModelContext *ctx);
int insert_new_Block(ModelConfig *cfg, ModelContext *ctx, int num_new_Block);
int Allocate_Memory(ModelConfig *cfg, ModelContext *ctx);
int Init_Stress(ModelConfig *cfg);
int make_gravi_body(ModelConfig *cfg, float *upper_hor, float *lower_hor, float *body_x, float *body_z);
int LES_matrix (ModelConfig *cfg, ModelContext *ctx, double **A, double *b, float *D, float *q, float *Dq, float *w, bool doing_visco);
float geoidanompolyg(
	float *y_pol, 			/*horizontal coordinate of line points clock-wise sorted*/
	float *z_pol, 			/*vertical coordinate of line points clock-wise sorted. z>0 downwards*/
	int numpoints, 			/*Number of line points*/
	float ym, float zm, 		/*position of measurement*/
	float dens_contrast);		/*Density of the body*/
float gravanompolyg(
	float *x_pol, 				/*x-y points of polygon clock-wise sorted*/
	float *z_pol, 				/*x-y points of polygon clock-wise sorted. z>0 downwards*/
	int numpoints, 				/*Number of polygon points*/
	float x_measure, float z_measure, 	/*position of measurement*/
	float dens_contrast);			/*Density of the body*/

float ReSort_Array (float *array, int *orden, int Nx);
float calculate_topo(ModelConfig *cfg, ModelContext *ctx, float *topo_new);
float calculate_sea_level(float current_time);
int read_file_YSE(ModelConfig *cfg);
int calculate_water_load(ModelConfig *cfg, ModelContext *ctx);
float Orographic_Precipitation (ModelConfig *cfg, ModelContext *ctx, int i, float windvel);
int Orographic_Precipitation_Evaporation_conservative (ModelConfig *cfg, ModelContext *ctx, float *precip_aux, float *evaporation, float windvel);


#include "tao+tisc.h"
