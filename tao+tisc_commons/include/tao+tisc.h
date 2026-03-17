#ifndef TAO_TISC_H
#define TAO_TISC_H

/*
	Common definitions for tAo and TISC
*/

// Enum for parameter types
typedef enum {
    PARAM_TYPE_INT,
    PARAM_TYPE_FLOAT,
    PARAM_TYPE_STRING, // For char arrays
    PARAM_TYPE_BOOL    // For bool variables
} ParamType;

// Structure to hold parameter metadata
typedef struct {
    const char *name;
    ParamType type;
    void *ptr; // Pointer to the actual variable
    size_t size; // Size for string types (MAXLENLINE, etc.)
    bool is_old_version; // Flag to indicate if this is an old/deprecated parameter name
} ParameterEntry;

/*
	Common definitions for tAo and TISC
*/

#define NMAXHEADERLINES	100
#define Match_Param_Replace_int(x, y, old)  if (!strcasecmp(str1, x)) {int i; nparams=1; \
		if (replace==0) {\
			y=atoi(str2); \
			if (old) {PRINT_INFO("Old-fashioned parameter: "x"\t%d", y);}\
			else if (show==1 && verbose_level>=3) fprintf(stdout, "\n"x"\t%d", y); }\
		else {	char newstr2[MAXLENLINE], *newline=NULL; \
			sprintf(newstr2, "%d", y); newline=replace_word(line, str2, newstr2); \
			fprintf(stdout, "%s", newline); free(newline);}\
	} 
#define Match_Param_Replace_flt(x, y, old)  if (!strcasecmp(str1, x)) {int i; nparams=1; \
		if (replace==0) {\
			y=atof(str2); \
			if (old) {PRINT_INFO("Old-fashioned parameter: "x"\t%e", y);}\
			else if (show==1 && verbose_level>=3) fprintf(stdout, "\n"x"\t%.4g", y); }\
		else {	char newstr2[MAXLENLINE], *newline=NULL; \
			sprintf(newstr2, "%.4g", y); newline=replace_word(line, str2, newstr2); \
			fprintf(stdout, "%s", newline); free(newline);}\
	} 
#define Match_Param_Replace_chr(x, y, old)  if (!strcasecmp(str1, x)) {int i; nparams=1; \
		if (replace==0) {\
			strcpy(y,str2); \
			if (old) {PRINT_INFO("Old-fashioned parameter: "x"\t%s", y);}\
			else if (show==1 && verbose_level>=3) fprintf(stdout, "\n"x"\t%s", y); }\
		else {	char newstr2[MAXLENLINE], *newline=NULL; \
			sprintf(newstr2, "%s", y); newline=replace_word(line, str2, newstr2); \
			fprintf(stdout, "%s", newline); free(newline);}\
	} 
#define Match_Param_int(x, y)  if (!strcasecmp(str1, x)) {y=atoi(str2);;  nparams=1; if (show==1 && verbose_level>=3) fprintf(stdout, "\n"x"\t%d", y); } 
#define Match_Param_flt(x, y)  if (!strcasecmp(str1, x)) {y=atof(str2);;  nparams=1; if (show==1 && verbose_level>=3) fprintf(stdout, "\n"x"\t%f", y); } 
#define Match_Param_char(x, y) if (!strcasecmp(str1, x)) {strcpy(y,str2); nparams=1; if (show==1 && verbose_level>=3) fprintf(stdout, "\n"x"\t%s", y); } 
#define Match_Param_int_old(x, y)  if (!strcasecmp(str1, x)) {y=atoi(str2);   nparams=1; PRINT_INFO("Old-fashioned parameter: "x"\t%d", y);}
#define Match_Param_flt_old(x, y)  if (!strcasecmp(str1, x)) {y=atof(str2);;  nparams=1; PRINT_INFO("Old-fashioned parameter: "x"\t%f", y);}
#define Match_Param_char_old(x, y) if (!strcasecmp(str1, x)) {strcpy(y,str2); nparams=1; PRINT_INFO("Old-fashioned parameter: "x"\t%s", y);}
#define Write_Open_Filename_Return(ext,type,retcond) {\
	char name[MAXLENFILE]; snprintf(name, sizeof(name), "%s"ext, projectname); remove(name); if (retcond) return (0); \
	if ((file=fopen(name,type))==NULL) {PRINT_WARNING("Could not open output file '%s'.",name);return 0;}\
	PRINT_INFO("Writing file '%s'.",name);};
#define Read_Open_Filename_Return(ext,type,txt) {\
	char name[MAXLENFILE];\
	snprintf(name, sizeof(name), "%s"ext, projectname);\
	if ((file = fopen(name,type)) == NULL) {PRINT_INFO("Cannot read "txt" input file '%s'.", name); return 0;}\
	PRINT_INFO("Reading "txt" at '%s'", name);};
#define Read_Header_File(file) {\
	char line[MAXLENLINE+200], *lineptr, str1[MAXLENLINE], str2[MAXLENLINE]; float value; int nlines=0, nparams=0, nread; bool switch_show=(verbose_level>=3)? true:false;\
	while ((lineptr=fgets(line, MAXLENLINE+200-1, file)) != NULL && nlines<NMAXHEADERLINES) {\
		nlines++; nread=sscanf(lineptr, "%s %s", str1, str2);\
		if (nread == 2) {\
			value = atof(str2);\
			Match_Param_int ( "mode_interp", 	mode_interp_local )\
	    	    	Match_Param_flt ( "add_random", 	add_random )\
			Match_Param_flt ( "z_default", 		z_default )\
		}\
		if (strcmp(str1, "thickness_distribution")==0) break;\
	}; rewind (file);} 



extern int nloads, n_sea_level_input_points, n_eros_level_input_points, n_record_times;
extern int i_first_Block_load, i_Block_insert, numBlocks, nwrotenfiles, run_type;

	
extern float Kerosdif, Keroseol, Ksedim, critical_slope, K_river_cap, erodibility, erodibility_sed, critical_stress;
extern float spl_m, spl_n, l_fluv_sedim, lost_rate, permeability, rain, Krain, relative_humidity, CXrain;
extern float rain_per, rain_amp, total_bedrock_eros_mass, total_sed_mass;



extern float zini, dt_record, sed_porosity, compact_depth, last_time_file_time, random_topo;
extern float **var_sea_level, **var_eros_level, **var_insolation, *horiz_record_time;

extern bool switch_file_out, switch_gradual, switch_topoest, switch_write_file_Blocks, deform_sed;


/* ===================================================================== */
/* DATA STRUCTURES                                                       */
/* ===================================================================== */

struct BLOCK { 		/*for TISC*/
	float	**thick;		/*Present thickness at each x point*/
	float	**detr_ratio;		/*Only used for sediment Blocks: % of detrital sediment (non carbonatic)*/
	float	**detr_grsize;		/*Only used for sediment Blocks: grain size of the detrital sediment*/
	float	**thickgypsum;		/*Only used for sediment Blocks: gypsum thickness*/
	float	**thickhalite;		/*Only used for sediment Blocks: halite thickness*/
	float	age;			/*Age of initial file reading*/
	float	density;		/*Density*/
	float	erodibility;		/*erosion parameter*/
	float	last_shift_x;		/*Previous x shift of Block*/
	float	last_shift_y;		/*Previous y shift of Block*/
	float	last_vel_time; 		/*Last time in which velocity changed*/
	float	shift_x;		/*Total x shift of Block*/
	float	shift_y;		/*Total y shift of Block*/
	float	time_stop;		/*Time to stop*/
	char 	type;			/*'T' means thin_sheet*/
	float	**vel_x;		/*Velocity in x direction*/
	float	**vel_y;		/*Velocity in y direction*/
	float	**visc;			/*Viscosity (only for thin sheet calculations)*/
	float	**viscTer;		/*Viscosity thermal term (only used for thin sheet calculations, in the first step)*/
};

struct BLOCK_1D {	/*for tAo*/
	float	*thick;			/*Present thickness at each x point*/
	float	*detr_ratio;		/*Only used for sediment Blocks: % of detrital sediment (non carbonatic)*/
	float	*detr_grsize;		/*Only used for sediment Blocks: grain size of the detrital sediment*/
	float	*thickgypsum;		/*Only used for sediment Blocks: gypsum thickness*/
	float	*thickhalite;		/*Only used for sediment Blocks: halite thickness*/
	float	age;			/*Age of initial file read*/
	float	density;		/*Density*/
	float	erodibility;		/*erosion parameter*/
	float	last_shift;		/*Espected shift (not affected by finite differences discretization)*/
	float	last_vel_time; 		/*Last time in which velocity changed*/
	float	shift;			/*Total horizontal shift of Block*/
	float	time_stop;		/*Time in wich Block will stop*/
	char 	type;			/*'T' means thin_sheet*/
	float	vel;			/*Velocity at wich Block moves*/
};

struct GRIDNODE {
	int row;
	int col;
};

struct DRAINAGE {
	int dr_row;		/*row of the node to where drains*/
	int dr_col;		/*column of the node to where drains*/
	float discharge;	/*water flow through the node [m3/s]*/
	float masstr;		/*sediment load: mass exiting the cell [kg/s]*/
	float grainsize;	/*average grain size of the sedload[m]*/
	float C_Ca;			/*concentration of Ca [kg/m3]*/
	float C_SO4;		/*concentration of SO4 [kg/m3]*/
	float C_Na;			/*concentration of Na [kg/m3]*/
	float C_Cl;			/*concentration of Cl [kg/m3]*/
	char type;		/*type (lake, river, sea, etc)*/
	int lake;		/*number of the lake: > 0 means is well defined; < 0 means is not still defined; 0 means it is not a lake*/
};

struct LAKE_INFO {		/*For lakes*/
	int n;			/*number of nodes INCLUDING SADDLES*/
	int *row;
	int *col;
	int n_sd;		/*number of saddles and transferring nodes*/
	int *row_sd;
	int *col_sd;
	float alt;		/*Altitude of the lake water level*/
	float vol;		/*Volume of the lake water body*/
	float mass_Ca;		/*dissolved mass of Ca in lake [kg]*/
	float mass_SO4;		/*dissolved mass of SO4 in lake [kg]*/
	float mass_Na;		/*dissolved mass of Na in lake [kg]*/
	float mass_Cl;		/*dissolved mass of Cl in lake [kg]*/
};

struct CS2D {
	float *horiz;
	float x;
	float y;
	float l;
};

struct DRAINAGE_1D {
	int dr;		/*row of the node to where drains*/
	float discharge;	/*water flow through the node [m3/s]*/
	float masstr;		/*sediment load: mass exiting the cell [kg/s]*/
	float grainsize;	/*average grain size of the sedload[m]*/
	float C_Ca;			/*concentration of Ca [kg/m3]*/
	float C_SO4;		/*concentration of SO4 [kg/m3]*/
	float C_Na;			/*concentration of Na [kg/m3]*/
	float C_Cl;			/*concentration of Cl [kg/m3]*/
	char type;		/*type (lake, river, sea, etc)*/
	int lake;		/*number of the lake: > 0 means is well defined; < 0 means is not still defined; 0 means it is not a lake*/
};

struct LAKE_INFO_1D {		/*For lakes*/
	int n;			/*number of nodes including saddles*/
	int *cell;
	int n_sd;		/*number of saddles and transferring nodes*/
	int *sd;
	float alt;		/*Altitude of the lake water level*/
	float vol;		/*Volume of the lake water body*/
	float mass_Ca;		/*dissolved mass of Ca in lake [kg]*/
	float mass_SO4;		/*dissolved mass of SO4 in lake [kg]*/
	float mass_Na;		/*dissolved mass of Na in lake [kg]*/
	float mass_Cl;		/*dissolved mass of Cl in lake [kg]*/
};

/*FUNCTION DECLARATIONS:*/
float 	*alloc_array		(int num_fil);
double	*alloc_array_dbl	(int num_fil);
int 	free_matrix		(float **, int);
int 	free_matrix_dbl		(double **, int);
int 	WriteAlmostDiagonalMatrix 	(double **A, double *b, int rango, char *filename, int NDs, int NDi);
int 	SolveAlmostDiagonalTriangularEquationSystem (double **A, double *b, int num_ecs, int NDsre, int NDire, float *x);
int 	TriangularizeAlmostDiagonalEquationSystem (double **A, double *b, int num_rows, int NDsre, int NDire);
float compaction(float phi0, float comp_depth, float z1, float z2);

char 	*replace_word(char *s,  char *old,  char *new);

#endif /* TAO_TISC_H */
