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
	char name[MAXLENFILE]; sprintf(name, "%s"ext, projectname); remove(name); if (retcond) return (0); \
	if ((file=fopen(name,type))==NULL) {PRINT_WARNING("Could not open output file '%s'.",name);return 0;}\
	PRINT_INFO("Writing file '%s'.",name);};
#define Read_Open_Filename_Return(ext,type,txt) {\
	char name[MAXLENFILE];\
	sprintf(name, "%s"ext, projectname);\
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
