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
	char line[MAXLENLINE+200], *lineptr, str1[MAXLENLINE], str2[MAXLENLINE]; float value; int nlines=0, nparams=0, nread; BOOL switch_show=(verbose_level>=3)? 1:0;\
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




int	nloads, 			/*Number of loads done*/
	n_sea_level_input_points,	/*Number of sea level input points*/
	n_eros_level_input_points, 	/*Number of input points of erosion level*/
	n_record_times, 		/*Number of horizons to record*/
	i_first_Block_load, 		/*Number of the first Block corresponding to a load*/
	i_Block_insert, 			/*Number of the Block in wich will be inserted the heigh of the next load*/
	numBlocks, 			/*Number of recorded Blocks*/
	nwrotenfiles, 			/*Number of wroten output files*/
	run_type; 			/*tao execution mode*/
	
float	Kerosdif, 		/*Diffusive erosion rate in m2/Ma, (e. g. D. Waltham & S. Hardy*/
	Keroseol, 		/*Background erosion rate in m/m/Ma*/
	Ksedim, 		/*Marine sedimentation rate in m/Ma and its maximum water depth of validity*/
	critical_slope, 	/*for landsliding*/
	K_river_cap, 		/*Constant of river transport capacity [kg/m3].*/
	erodibility, 		
	erodibility_sed, 	
	spl_m, spl_n, 		/*exponents of the stream power law over Q and S; m/n is ca. 0.5*/
	l_fluv_sedim, 		/*Length scale of fluvial sedimentation */
	lost_rate, 		/*Percent of lost water per unit length */
	rain, Krain, 		/*[m3/s/m2], [m3/s/m2/m] or m/s*/
	relative_humidity, 	/*Relative humidity at the upwind boundary (incoming air rel. humidity) [no units]*/
	CXrain, 		/*[m3/s/m2], [m3/s/m2/m], [m], [m]*/
	total_bedrock_eros_mass,
	total_sed_mass; 


float	zini, 				/*altitude of the initial plate position over the sea level [m]*/
	dt_record, 				/*maximum Time interval between automatically-generated sediment Blocks [s]*/
	sed_porosity, 
	compact_depth, 
	last_time_file_time, 
	random_topo=0, 			/*Maximum random variation of initial topo*/
	**var_sea_level, 		/*Array nx2 with sea level along time*/
	**var_eros_level, 		/*Variations of the erosion level (height dividing erosion and sedimentation, measured from sea_level)*/
	*horiz_record_time;		/*Array with times in wich to record horizons*/
	
BOOL	switch_dt_output=NO, 
	switch_file_out=NO, 
	switch_gradual, 		/*YES to distribute the load between the Time of reading file and the Time of the following one*/
	switch_topoest, 		/*YES if load files give topographic loads which must stay at zero level while the deflection room is filled up with 'densinfill' density material*/
	switch_write_file_Blocks, 	/*YES if Blocks (profile *.pfl) file is to be written*/
	deform_sed;			/*YES to deform sediment automatically, based on Blocks motion*/


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
