/*
INPUT/OUTPUT  FUNCTIONS  FOR  tao.c
*/
#include "tao.h"
#include <stdlib.h> // For EXIT_FAILURE
#include <string.h> // For strncpy

int read_file_horiz_record_time()
{
	/*
		Reads file containig times of horizons recording named 'projectname.REC'
	*/

	int 	i;
	FILE 	*file;
	char 	filename[MAXLENFILE];
	float	*aux1;

	Read_Open_Filename_Return(".REC", "rt", "Horizon recording times")

	n_record_times=0;
	aux1 = (float*)calloc(nmax_input_points, sizeof(float));
	if (!aux1) { PRINT_ERROR("Memory allocation failed for aux1 in %s.", __func__); return 0; }

	char line[MAXLENLINE];
	while (fgets(line, sizeof(line), file) != NULL) {
		// Skip comments and empty lines
		if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
		if (sscanf(line, "%f", &aux1[n_record_times]) != 1) continue;
		n_record_times++;

		if (n_record_times>=nmax_input_points-1 ) {
			PRINT_ERROR("Too many points (>%d) in horizon recording times file.", n_record_times-1);
			break;
		}
	}
	fclose(file);
	if (verbose_level>=1) fprintf(stdout, "\nHorizon recording times at '%s'. %d times were read.", filename, n_record_times);

	horiz_record_time = (float*)calloc(n_record_times, sizeof(float));
	if (!horiz_record_time) { PRINT_ERROR("Memory allocation failed for horiz_record_time in %s.", __func__); free(aux1); return 0; }

	for (i=0; i<n_record_times; i++) {
		horiz_record_time[i] = aux1[i]*Matosec;
	}
	
	/*desactivate automatic sediment Blocks generation*/
	if (n_record_times) dt_record=0;
	
	free(aux1); 
	return(1);
}




int read_file_initial_topo(float *h)
{
	int 	i ;
	FILE 	*file ;

	/*READS INITIAL TOPOGRAPHY IF EXISTS FILE WITH NAME 'projectname.ZINI'*/

	Read_Open_Filename_Return(".ZINI", "rt", "Initial topography")

	readinterplin(file, h, Nx, x0, xf) ;

	fclose(file);
	return (1);
}



int read_file_initial_deflection(float *w)
{
	int 	i ;
	FILE 	*file ;

	/*READS INITIAL DEFLECTION IF DOES EXIST FILE WITH NAME 'projectname.WINI'*/

	Read_Open_Filename_Return(".WINI", "rt", "Initial deflection")

	readinterplin(file, w, Nx, x0, xf) ;

	fclose(file);
	return (1);
}




int read_file_Crust_Thick(float crust_thick_default)
{
	/*CARGA EL FICHERO DE ESPESOR CORTICAL, DE NOMBRE  'projectname.CRUST'
			Interpola entre los puntos del fichero.*/


	int 	i;
	FILE 	*file;

	for (i=0;i<Nx;i++) {crust_thick[i] = crust_thick_default;}

	Read_Open_Filename_Return(".CRUST", "rt", "Crustal thickness")

	readinterplin(file, crust_thick, Nx, x0, xf) ;

	fclose(file);
	return(1);
}



int read_file_Upper_Crust_Thick(float crust_thick_supdefault)
{
	/*CARGA EL FICHERO DE ESPESOR DE CORTEZA SUPERIOR, DE NOMBRE  'projectname.UCRUST'
			Interpola entre los puntos del fileero.*/


	int 	i;
	FILE 	*file;

	for (i=0;i<Nx;i++) {upper_crust_thick[i] = upper_crust_thick_default;}

	Read_Open_Filename_Return(".UCRUST", "rt", "Upper crustal thickness")

	readinterplin(file, upper_crust_thick, Nx, x0, xf) ;

	for (i=0; i<Nx; i++) upper_crust_thick[i] = MIN_2(upper_crust_thick[i], crust_thick[i]) ;

	fclose(file);
	return(1);
}



int read_file_parameters (int show, int reformat) 
{
	int	nread, nparams=0, nline=0, verbose_level_ant=verbose_level;
	char 	*lineptr, str1[MAXLENLINE], str2[MAXLENLINE], 
		line[MAXLENLINE+200], PRMfilename[MAXLENFILE];
	FILE 	*file;
	bool	switch_matched_vers=false;

	/*
	READ THE PARAMETERS FILE NAMED  'projectname.PRM'
	You have an explanation of these parameters in the example project 
	parameters file (doc/template.PRM), where you have also the name of 
	the variable related to each parameter. The use of most variables is 
	described in the include files '*.h'.
	*/

	sprintf(PRMfilename, "%s.PRM", projectname);
	if (show && verbose_level>=3) fprintf(stdout, "\nCurrent tAo project: %s", projectname);
	if ((file = fopen(PRMfilename, "rt")) == NULL) {
		PRINT_ERROR("Can't open parameters file '%s'.\n", PRMfilename);
		return 0;
	}

	x0 = NO_DATA; xf = NO_DATA;
	if (show) fprintf(stdout, "\nParameters at '%s'.", PRMfilename); fflush(stdout); // Flush stdout to ensure message appears before potential errors
	while (fgets(line, sizeof(line), file) != NULL) { // Read line by line
		int status;
		status=0;
		// Skip comments and empty lines (already handled by fgets and sscanf)
		if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
		// Attempt to read two strings (parameter name and value)
		if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
		nread=sscanf(line, "%s %s", str1, str2);
		if (nread >= 2) {
			if (!strcmp(str1, "version")) {
				strcpy(version_input, str2);
				if (!strcmp(version_input, version)) {
					switch_matched_vers = true;
					nparams++;
				}
				if ((show && verbose_level>=3) || (reformat && verbose_level>=3)) PRINT_INFO("Input file version: %s", version_input);
				if (reformat) {fprintf(stdout, "version\t%s\n", version); }
			}
			status=match_parameter(str1, str2, show, reformat, line);
			nparams += status;
			if (!strcmp(str1, "version")) status=1;
		}
		/*If no parameter matched then just reproduce the entire line*/
		if (reformat==1 && !status) fprintf(stdout, "%s", line);
		nline++;
	} // End while
	if (verbose_level_ant >= 2) fprintf(stdout, " (%d parameters matched)", nparams);
	if (!switch_matched_vers) {
		if (verbose_level_ant >= 2)
		    fprintf(stderr, "\nInfo: Present version '%s' not matched in PRM file. See 'tisc/doc/template.PRM'.", version);
		if (nparams < 2) {
		    PRINT_ERROR("\aWrong format in parameters file '%s'. Only %d parameters were matched. Current version is '%s'. "
			"\nSee example file 'tisc/doc/template.PRM'.\nEND.\n", 
			PRMfilename, nparams, version);
		    exit(0); 
		}
	}

	if (x0==NO_DATA) x0 = xmin;
	if (xf == NO_DATA) xf = xmax;
	if (Kerosdif>1e5) PRINT_WARNING("Kerosdif expected in Blocks [m2/yr] !") ;

	fclose(file);
	return(1);
}




int read_file_resume(char *filename)
{
	int 	i, numBlocks_aux, i_Block_insert_aux, run_type_aux, end_check;
	FILE 	*file;
	char  	version_aux[LENGTHVERS], projectname_aux[MAXLENFILE];

	/*
	  READS A BINARY FILE WITH ALL THE INFORMATION 
	  REQUIRED TO RESTART THE PROGRAM
	*/

	if ((file = fopen(filename,"rb")) == NULL) {PRINT_ERROR("Cannot read ""Resume"" input file '%s'.\n", filename); exit(EXIT_FAILURE);}\
	if (verbose_level>=1) fprintf(stdout, "\n""Model-run data at '%s'", filename);

	/*Defined in universal.h:*/
	fread(&Nx, 		sizeof(int),		1, 	file);
	fread(&Nz, 		sizeof(int),		1, 	file);
	fread(&verbose_level, 	sizeof(int),		1, 	file);

	fread(&dx, 		sizeof(float),		1, 	file);
	fread(&dz, 		sizeof(float),		1, 	file);

	fread(version_aux,	sizeof(char),		LENGTHVERS, 	file);
	if (strcmp(version, version_aux)) PRINT_WARNING("restart file '%s' does not match present tAo version '%s'.", filename, version);
	fread(version_input,	sizeof(char),		LENGTHVERS, 	file);

	fread(&switch_geograph_coor, sizeof(bool),	1, 	file);
	fread(&switch_ps, 	sizeof(bool),		1, 	file);
	fread(&switch_write_file, sizeof(bool),		1, 	file);


	/*Defined in geomodel.h:*/
	fread(&grav_anom_type, 	sizeof(int),		1, 	file);
	fread(&isost_model, 	sizeof(int),		1, 	file);

	fread(&Te_default, 	sizeof(float),		1, 	file);
	fread(&crust_thick_default, sizeof(float),	1, 	file);
	fread(&upper_crust_thick_default, sizeof(float),1, 	file);
	fread(&densasthen, 	sizeof(float),		1, 	file);
	fread(&densmantle, 	sizeof(float),		1, 	file);
	fread(&denssedim, 	sizeof(float),		1, 	file);
	fread(&denscrust, 	sizeof(float),		1, 	file);
	fread(&densenv, 	sizeof(float),		1, 	file);
	fread(&densinfill, 	sizeof(float),		1, 	file);
	fread(&sea_level, 	sizeof(float),		1, 	file);
	fread(&temp_sea_level, 	sizeof(float),		1, 	file);
	fread(&Time, 		sizeof(float),		1, 	file);
	fread(&Timefinal, 	sizeof(float),		1, 	file);
	fread(&Timeini, 	sizeof(float),		1, 	file);
	fread(&dt, 		sizeof(float),		1, 	file);
	fread(&dt_eros, 		sizeof(float),		1, 	file);
	fread(&tau, 		sizeof(float),		1, 	file);

	fread(projectname_aux, 	sizeof(char),	MAXLENFILE, 	file);
	if (strcmp(projectname, projectname_aux)) fprintf(stdout, "\nERROR: restart file '%s' does not match present tAo project name '%s'.", filename, version);
	fread(gif_geom, 	sizeof(char),	MAXLENLINE, 	file);

	fread(&water_load, 	sizeof(bool),		1, 	file);


	/*Defined in tao+tisc.h:*/
	fread(&boundary_conds, 	sizeof(int),		1, 	file);
	fread(&nloads, 		sizeof(int),		1, 	file);
	fread(&nmax_input_points, sizeof(int),		1, 	file);
	fread(&n_sea_level_input_points, sizeof(int),	1, 	file);
	fread(&n_eros_level_input_points, sizeof(int),	1, 	file);
	fread(&n_record_times, 	sizeof(int),		1, 	file);
	fread(&i_first_Block_load, sizeof(int),	1, 	file);
	fread(&i_Block_insert_aux, sizeof(int),		1, 	file);
	fread(&numBlocks_aux, 	sizeof(int),		1, 	file);
	fread(&nwrotenfiles, 	sizeof(int),		1, 	file);
	fread(&run_type_aux, 	sizeof(int),		1, 	file);

	fread(&zini, 		sizeof(float),		1, 	file);
	fread(&dt_record, 	sizeof(float),		1, 	file);
	fread(&sed_porosity, 	sizeof(float),		1, 	file);
	fread(&compact_depth, 	sizeof(float),		1, 	file);
	fread(&Kerosdif, 	sizeof(float),		1, 	file);
	fread(&last_time_file_time, 	sizeof(float),		1, 	file);
	fread(&random_topo, 	sizeof(float),		1, 	file);

	fread(&switch_file_out, 	sizeof(bool),		1, 	file);
	fread(&switch_gradual, 	sizeof(bool),		1, 	file);
	fread(&switch_insert_load, 	sizeof(bool),		1, 	file);
	fread(&switch_topoest, 		sizeof(bool),		1, 	file);
	fread(&switch_write_file_Blocks, sizeof(bool),		1, 	file);
	fread(&deform_sed, sizeof(bool),		1, 	file);

	/*Defined in tao.h:*/
	fread(&imomentmax, 	sizeof(int),		1, 	file);
	fread(&nx_temp_input, 	sizeof(int),		1, 	file);
	fread(&nbasins, 	sizeof(int),		1, 	file);
	fread(&nlakes,	 	sizeof(int),		1, 	file);
	fread(&n_image, 	sizeof(int),		1, 	file);
	fread(&hydro_model,	sizeof(int),		1, 	file);
	fread(&erosed_model,	sizeof(int),		1, 	file);
	fread(&eros_bound_cond,	sizeof(char),		2, 	file);

	fread(&x0, 		sizeof(float),		1, 	file);
	fread(&xf, 		sizeof(float),		1, 	file);
	fread(&xmin, 		sizeof(float),		1, 	file);
	fread(&xmax, 		sizeof(float),		1, 	file);
	fread(&zmin, 		sizeof(float),		1, 	file);
	fread(&zmax, 		sizeof(float),		1, 	file);
	fread(&horz_force, 	sizeof(float),		1, 	file);
	fread(&vert_force, 	sizeof(float),		1, 	file);
	fread(&appmoment, 	sizeof(float),		1, 	file);
	fread(&Keroseol, 	sizeof(float),		1, 	file);
	fread(&Ksedim, 		sizeof(float),		1, 	file);
	fread(&critical_slope, 	sizeof(float),		1, 	file);
	fread(&K_river_cap, 	sizeof(float),		1, 	file);
	fread(&erodibility, 	sizeof(float),		1, 	file);
	fread(&erodibility_sed, sizeof(float),		1, 	file);
	fread(&l_fluv_sedim, 	sizeof(float),		1, 	file);
	fread(&lost_rate, 	sizeof(float),		1, 	file);
	fread(&evaporation_ct, 	sizeof(float),		1, 	file);
	fread(&riverbasinwidth, 	sizeof(float),		1, 	file);
	fread(&rain, 		sizeof(float),		1, 	file);
	fread(&Krain,	 	sizeof(float),		1, 	file);
	fread(&CXrain, 		sizeof(float),		1, 	file);

	fread(&switch_strs_history, 	sizeof(bool),		1, 	file);
	fread(&switch_YSE_file, 	sizeof(bool),		1, 	file);


	/*Arrays:*/
	ModelConfig cfg = {0}; ModelContext ctx = {0};
	cfg.Nx = Nx; cfg.Nz = Nz; cfg.dx = dx; cfg.dz = dz;
	cfg.x0 = x0; cfg.xf = xf; cfg.xmin = xmin; cfg.xmax = xmax;
	cfg.hydro_model = hydro_model;
	cfg.erosed_model = erosed_model;
	cfg.verbose_level = verbose_level;
	cfg.densenv = densenv; cfg.denssedim = denssedim;
	cfg.denscrust = denscrust; cfg.denswater = denswater; 
	cfg.densasthen = densasthen; cfg.densmantle = densmantle;
	Allocate_Memory(&cfg, &ctx);
	fread(w, 		sizeof(float),	Nx, 	file);
	fread(D, 		sizeof(float),	Nx, 	file);
	fread(q, 		sizeof(float),	Nx, 	file);
	fread(Dw, 		sizeof(float),	Nx, 	file);
	fread(Dq, 		sizeof(float),	Nx, 	file);
	fread(h_water, 		sizeof(float),	Nx, 	file);
	fread(h_last_unit, 	sizeof(float),	Nx, 	file);
	fread(Te, 		sizeof(float),	Nx, 	file);
	fread(crust_thick, 	sizeof(float),	Nx, 	file);
	fread(upper_crust_thick,sizeof(float),	Nx, 	file);
	fread(topo, 		sizeof(float),	Nx, 	file);
	fread(Blocks_base, 	sizeof(float),	Nx, 	file);

	horiz_record_time = calloc(n_record_times, sizeof(float));
	fread(horiz_record_time, sizeof(float), n_record_times, file);

	var_sea_level = calloc(n_sea_level_input_points, sizeof(float *));
	for (i=0; i<n_sea_level_input_points; i++) {
		var_sea_level[i] = calloc(2, sizeof(float));
		fread(var_sea_level[i], sizeof(float), 2, file);
	}
	var_eros_level = calloc(n_eros_level_input_points, sizeof(float *));
	for (i=0; i<n_eros_level_input_points; i++) {
		var_eros_level[i] = calloc(2, sizeof(float));
		fread(var_eros_level[i], sizeof(float), 2, file);
	}

	if (isost_model>=3 && !switch_YSE_file) {
		Temperature = alloc_matrix(Nx, Nz);
		if (fread(Temperature[0], sizeof(float), Nx * Nz, file) != (size_t)(Nx * Nz)) { PRINT_ERROR("Failed to read Temperature array."); goto error_read_resume; }
	}
	if (isost_model>=3) {
		stress = alloc_matrix(Nx, Nz);
		yieldcompres = alloc_matrix(Nx, Nz);
		yieldextens = alloc_matrix(Nx, Nz);
		if (fread(stress[0], sizeof(float), Nx * Nz, file) != (size_t)(Nx * Nz)) { PRINT_ERROR("Failed to read stress array."); goto error_read_resume; }
		if (fread(yieldcompres[0], sizeof(float), Nx * Nz, file) != (size_t)(Nx * Nz)) { PRINT_ERROR("Failed to read yieldcompres array."); goto error_read_resume; }
		if (fread(yieldextens[0], sizeof(float), Nx * Nz, file) != (size_t)(Nx * Nz)) { PRINT_ERROR("Failed to read yieldextens array."); goto error_read_resume; }
	}


	for (i=0; i<numBlocks_aux; i++) {
		float *ptr1, *ptr2, *ptr3;
		ModelConfig tmp_cfg = {0}; ModelContext tmp_ctx = {0};
		tmp_cfg.Nx = Nx; tmp_ctx.numBlocks = numBlocks; tmp_ctx.Time = Time;
		insert_new_Block(&tmp_cfg, &tmp_ctx, numBlocks);
		numBlocks = tmp_ctx.numBlocks;
		ptr1 = Blocks[numBlocks-1].thick; 
		fread(&Blocks[numBlocks-1], sizeof(struct BLOCK_1D), 1, file);
		Blocks[numBlocks-1].thick = ptr1;
	}
	if (numBlocks_aux != numBlocks) PRINT_ERROR("in '%s', %d Blocks?!", filename, numBlocks_aux);
	i_Block_insert=i_Block_insert_aux;
	for (i=0; i<numBlocks; i++) { // Blocks[i].thick is a 1D array, so read directly
		if (fread(Blocks[i].thick, sizeof(float), Nx, file) != (size_t)Nx) { PRINT_ERROR("Failed to read Blocks[%d].thick array.", i); goto error_read_resume; } // Check return value
	}
	for (i=0; i<numBlocks; i++) {
	    if (Blocks[i].type == 'S') {
		Blocks[i].detr_ratio = calloc(Nx, sizeof(float));
		Blocks[i].detr_grsize = calloc(Nx, sizeof(float));
		fread(Blocks[i].detr_ratio, 	sizeof(float),	Nx, 	file);
		fread(Blocks[i].detr_grsize, 	sizeof(float),	Nx, 	file);
	    }
	}

	if (erosed_model) { // These are 1D arrays, so read directly
		if (fread(eros_now, sizeof(float), Nx, file) != (size_t)Nx) { PRINT_ERROR("Failed to read eros_now array."); goto error_read_resume; } // Check return value
		if (fread(total_erosion, sizeof(float), Nx, file) != (size_t)Nx) { PRINT_ERROR("Failed to read total_erosion array."); goto error_read_resume; } // Check return value
	}
	if (hydro_model) { // These are 1D arrays, so read directly
		if (fread(precipitation, sizeof(float), Nx, file) != (size_t)Nx) { PRINT_ERROR("Failed to read precipitation array."); goto error_read_resume; } // Check return value
		if (fread(evaporation, sizeof(float), Nx, file) != (size_t)Nx) { PRINT_ERROR("Failed to read evaporation array."); goto error_read_resume; } // Check return value
		Lake = calloc (nlakes+1, sizeof(struct LAKE_INFO_1D));
		fread(Lake, sizeof(struct LAKE_INFO_1D), nlakes+1, file);
		for (int j=1; j<=nlakes; j++) {
			Lake[j].cell = calloc (Lake[j].n, sizeof(int));
			fread(Lake[j].cell, sizeof(int), Lake[j].n, file);
			Lake[j].sd = calloc (Lake[j].n_sd, sizeof(int));
			fread(Lake[j].sd, sizeof(int), Lake[j].n_sd, file);
		}
	}

	fread(&end_check,	sizeof(int),		1, 	file);
	if (end_check != 12345) {
		PRINT_ERROR("\achecking the end of resume file (%d).\n", end_check);
		goto error_read_resume;
	}
	else {
		PRINT_INFO("Check of resume file '%s' is ok.", filename);
	}

	if (switch_file_out){
		char filename[MAXLENLINE];
		FILE *file2;
		sprintf(filename, "%s.out", projectname);
		if ((file2 = fopen(filename, "a")) == NULL) {
			PRINT_ERROR("Cannot open standard output file %s.\n", filename);
		}
		else {
			PRINT_INFO("standard output redirected to %s.\n", filename);
		}
		stdout=file2;
	}

	fclose(file);
	return 1;
error_read_resume:
	return(1);
}



int read_file_sea_level()
{
	/*
		Reads file with sea level along time named 'projectname.SLV'
	*/

	int 	i, j;
	FILE 	*file;
	float	*aux1, *aux2, *aux3;

	sea_level = 0;

	Read_Open_Filename_Return(".SLV", "rt", "Sea level")

	n_sea_level_input_points=n_eros_level_input_points=0;
	aux1 = (float*)calloc(nmax_input_points, sizeof(float));
	if (!aux1) { PRINT_ERROR("Memory allocation failed for aux1 in read_file_sea_level."); return 0; }
	aux2 = (float*)calloc(nmax_input_points, sizeof(float));
	if (!aux2) { PRINT_ERROR("Memory allocation failed for aux2 in read_file_sea_level."); free(aux1); return 0; }
	aux3 = (float*)calloc(nmax_input_points, sizeof(float));
	if (!aux3) { PRINT_ERROR("Memory allocation failed for aux3 in read_file_sea_level."); free(aux1); free(aux2); return 0; }
	for (i=0; i<nmax_input_points; i++) aux3[i]=NO_DATA;

	char line[MAXLENLINE];
	while (fgets(line, sizeof(line), file) != NULL) {
		// Skip comments and empty lines
		if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
		// Use a temporary variable for nfields to avoid issues with macro expansion
		// if (sscanf(line, "%f %f %f", &aux1[n_sea_level_input_points], &aux2[n_sea_level_input_points], &aux3[n_sea_level_input_points]) != 3) continue;
		int nfields = sscanf(line, "%f %f %f", &aux1[n_sea_level_input_points], &aux2[n_sea_level_input_points], &aux3[n_sea_level_input_points]);
		if (nfields < 2) continue;
		n_sea_level_input_points++;
		if (aux3[n_sea_level_input_points-1]!=NO_DATA) 
			n_eros_level_input_points++;
		if (n_sea_level_input_points>=nmax_input_points-1 ) {
			PRINT_ERROR("Too many points (%d) in sea level file.", n_sea_level_input_points);
			break;
		}
	}
	fclose(file); 
	if (verbose_level >= 1) fprintf(stdout, "\nSea level variations file contains %d points.", n_sea_level_input_points);
	var_sea_level = (float**)calloc(n_sea_level_input_points, sizeof(float *));
	if (!var_sea_level) { PRINT_ERROR("Memory allocation failed for var_sea_level."); free(aux1); free(aux2); free(aux3); return 0; }
	var_eros_level = (float**)calloc(n_eros_level_input_points, sizeof(float *));
	if (!var_eros_level) { PRINT_ERROR("Memory allocation failed for var_eros_level."); free(aux1); free(aux2); free(aux3); free(var_sea_level); return 0; }
	for (i=0, j=0; i<n_sea_level_input_points; i++) {
		var_sea_level[i] = (float*)calloc(2, sizeof(float));
		if (!var_sea_level[i]) { PRINT_ERROR("Memory allocation failed for var_sea_level[%d].", i); free(aux1); free(aux2); free(aux3); free(var_sea_level); free(var_eros_level); return 0; }
		var_sea_level[i][0] = aux1[i]*Matosec;
		var_sea_level[i][1] = aux2[i];
		if (aux3[i]!=NO_DATA) {
			var_eros_level[j] = (float*)calloc(2, sizeof(float));
			if (!var_eros_level[j]) { PRINT_ERROR("Memory allocation failed for var_eros_level[%d].", j); free(aux1); free(aux2); free(aux3); free(var_sea_level); free(var_eros_level); return 0; }
			var_eros_level[j][0] = aux1[i]*Matosec;
			var_eros_level[j][1] = aux3[i];
			j++;
		}
	}
	free(aux1); free(aux2); free(aux3);
	return(1);
}



int read_file_Te()
{
	/*
		Reads elastic thickness file named 'projectname.EET'
		Linearly interpolates between given points.
	*/

	int 	i;
	FILE 	*file;
	float 	Dref=ET2RIG(Te_default) ;
	char 	filename[MAXLENFILE];

	sprintf(filename, "%s.EET", projectname);
	if ((file = fopen(filename, "rt")) == NULL) {
		PRINT_WARNING("Cannot read EET file '%s'. Constant rigidity = %3.4e", filename, Dref);
		if (isost_model>=3 && Te_default==-1) {
			/*Read the eet file as initial value for Te in Rheoiter*/
			sprintf(filename, "%s.eeth", projectname);
			if ((file = fopen(filename, "rt")) == NULL) {
				PRINT_WARNING("EET file '%s' not found. Constant rigidity = %3.4e", filename, Dref);
			}
			else {
				/*Note that x in this file is in km!*/
				if (verbose_level>=1) fprintf(stdout, "\nInitial elastic thickness read from '%s'.", filename);
				readinterplin(file, Te, Nx, x0/1000, xf/1000);
				fclose(file);
				for (i=0;i<Nx;i++) {D[i]=ET2RIG(Te[i]);}
				return(1);
			}
		}
		for (i=0;i<Nx;i++) {D[i]=Dref; Te[i]=Te_default; }
		return(0);
	}
	if (verbose_level>=1) fprintf(stdout, "\nElastic thickness at '%s'.", filename);
	readinterplin(file, Te, Nx, x0, xf);
	for (i=0;i<Nx;i++) {D[i]=ET2RIG(Te[i]);}
	fclose(file);
	return(1);
}



int read_file_Temperature(ModelConfig *cfg, ModelContext *ctx)
{
	/*READS TEMPERATURE FILE 'projectname.TMP'*/

	int 	i, n_fld, ix, iz, i_x_temp, i_zinp, 
		nz_temp_input, n_max_input_temp=100;
	FILE 	*file;
	bool	last_line=false;
	char 	filename[MAXLENFILE], *lin, 
		linea[MAXLENLINE];
	float	a, b, x_first, x_last, 
		*temp, *z_temp, x, z,
		**temp_input,		/*Temperature at input points of .TMP file [�C]*/
		*x_temp_input; 		/*positions of the given geotherms of .TMP file*/

	if (isost_model<3 || switch_YSE_file) return(0);

	Temperature = alloc_matrix(cfg->Nx, cfg->Nz);

	sprintf(filename, "%s.TMP", projectname) ;
	if ((file = fopen(filename, "rt")) == NULL) {
		PRINT_WARNING("Temperature file '%s' not found. Plasticity unswitched.", filename);
		isost_model = 1 ;
		return(0);
	}
	else if (verbose_level>=1) fprintf(stdout, "\nTemperature at '%s'.", filename);

	x_temp_input =	(float *) calloc(n_max_input_temp, sizeof(float));
	temp_input =	(float **)calloc(n_max_input_temp, sizeof(float *));
	z_temp = (float *) calloc(500, sizeof(float));
	temp =   (float *) calloc(500, sizeof(float));
	nz_temp_input = 0; 
	nx_temp_input = 0; 
	do {
		lin=fgets(linea, MAXLENLINE-1, file);
		if (!lin) last_line=true;
		n_fld = sscanf(linea, "%f %f", &a, &b);
		if ((n_fld==2 && nx_temp_input==0)) { 
			nx_temp_input=1;
			x_temp_input[0]=xmin;
			temp_input[nx_temp_input-1] =   (float *) calloc(cfg->Nz, sizeof(float));
			if (verbose_level>=1) fprintf(stdout, " Temperature laterally constant.");
		}
		if ((n_fld==1 && nx_temp_input>0) || !lin) {
			/*Interpolates vertically the geotherm (temperature-z)*/
			for (iz=0; iz<cfg->Nz; iz++) {
				z = iz*cfg->dz;
				for (i_zinp=0; i_zinp<nz_temp_input-1; i_zinp++)
					if (z<=z_temp[i_zinp+1]) break;
				temp_input[nx_temp_input-1][iz] = temp[i_zinp] + (z-z_temp[i_zinp]) * 
					(temp[i_zinp+1]-temp[i_zinp]) / (z_temp[i_zinp+1]-z_temp[i_zinp]);
				if (z <= z_temp[0]) 	temp_input[nx_temp_input-1][iz]=temp[0]; 
				if (z >= z_temp[nz_temp_input-1]) temp_input[nx_temp_input-1][iz]=temp[nz_temp_input-1]; 
			}
		}
		if (n_fld==1) { 
			nx_temp_input++; 
			nz_temp_input = 0; 
			x_temp_input[nx_temp_input-1] = a;
			temp_input[nx_temp_input-1] =   (float *) calloc(cfg->Nz, sizeof(float));
		}
		if (n_fld==2) { 
			nz_temp_input++; 
			z_temp[nz_temp_input-1] = a; 
			temp[nz_temp_input-1]   = b;
		}
	} while(!last_line);
	fclose(file);
	
	x_first = x_temp_input[0];
	x_last = x_temp_input[nx_temp_input-1];
	if (verbose_level>=1 && nx_temp_input>1) fprintf(stdout, " %d geotherms read from x=%.1f to %.1f km", nx_temp_input, x_first/1000, x_last/1000);

	/*Interpolates horizontally the geotherm (temperature-z)*/
	for (ix=0; ix<cfg->Nx; ix++) {
		x=ix*cfg->dx+cfg->x0;
		for (i_x_temp=0; i_x_temp<nx_temp_input-1; i_x_temp++)
			 if (x<=x_temp_input[i_x_temp+1]) break;
		for (iz=0; iz<cfg->Nz; iz++) {
			if (nx_temp_input>1) {
			    if (x>x_temp_input[0] && x<x_temp_input[nx_temp_input-1]) {
				Temperature[ix][iz] = temp_input[i_x_temp][iz] +
					(x-x_temp_input[i_x_temp]) *
					(temp_input[i_x_temp+1][iz]-temp_input[i_x_temp][iz]) /
					(x_temp_input[i_x_temp+1]-x_temp_input[i_x_temp]);
			    }
			    if (x<=x_temp_input[0]) 		
				Temperature[ix][iz] = temp_input[0][iz];
			    if (x>=x_temp_input[nx_temp_input-1]) 	
				Temperature[ix][iz] = temp_input[nx_temp_input-1][iz];
			}
			else 	Temperature[ix][iz]=temp_input[0][iz];
		}
	}
	
	write_file_Temperature_initial(cfg, ctx);
	
	free(z_temp);
	free(temp);
	free(x_temp_input);
	free(temp_input);
	return(1);
}





/******************************************************************************/






int write_file_erosed (ModelConfig *cfg, ModelContext *ctx, float *total_erosion)
{
	int 	i, ix_min, ix_max;
	float	x ;
	FILE 	*file;

	/*PRINTS A FILE with surface transport results*/

	Write_Open_Filename_Return (".eros", "wt", !cfg->erosed_model);

	fprintf(file, "#Time: %.2fMy\n#x[km]\terosion[m]\teros_rate[m/My]\ttopo[m]", ctx->Time/Matosec);
		if (cfg->erosed_model>=2) 
			fprintf(file, "\tdischarge[m3/s]\tsedload[kg/s]\ttype\tprecipt[l/m2/yr]\tevap[l/m2/yr]");
	ix_min = MAX_2((cfg->xmin-cfg->x0-.1*cfg->dx)/cfg->dx, 0) ;	ix_max = MIN_2(floor((cfg->xmax-cfg->x0+.1*cfg->dx)/cfg->dx) + 2, cfg->Nx);
	for(i=ix_min; i<ix_max; i++) {
		float drainage_aux;
		x=cfg->x0+cfg->dx*i;
		fprintf(file, "\n%7.2f\t%8.1f\t%8.2f\t%7.1f", 
			x/1000, total_erosion[i] / cfg->dx/cfg->riverbasinwidth/cfg->denscrust, 
			eros_now[i]/(ctx->dt/Matosec) / cfg->dx/cfg->riverbasinwidth/cfg->denscrust, ctx->topo[i] );
		if (cfg->erosed_model>=2) {
			drainage_aux = drainage[i].discharge;
			if (drainage[i].lake) if (Lake[drainage[i].lake].n_sd) 
				drainage_aux = drainage[Lake[drainage[i].lake].sd[0]].discharge;
			fprintf(file, "\t%8.6f\t%8.6f\t%c\t%8.1f\t%8.1f", 
				drainage_aux, 
				drainage[i].masstr, drainage[i].type, 
				precipitation[i]*secsperyr*1e3, evaporation[i]*secsperyr*1e3);
		}
	}
	fprintf(file, "\n");
	fclose(file);
	return(1);
}


int write_file_grav_anom (ModelConfig *cfg, ModelContext *ctx, float *gravanom, float *geoidanom)
{
	int 	i, ix_min, ix_max;
	float	x ;
	FILE 	*file;

	/*PRINTS A FILE WITH GRAVITY ANOMALY IN mGal*/

	Write_Open_Filename_Return (".xg", "wt", !grav_anom_type);

	fprintf(file, "#Time: %.2fMy\n#x[km]\tdg[mGal]  N[m]", ctx->Time/Matosec);
	ix_min = MAX_2((cfg->xmin-cfg->x0-.1*cfg->dx)/cfg->dx, 0) ;	ix_max = MIN_2(floor((cfg->xmax-cfg->x0+.1*cfg->dx)/cfg->dx) + 2, cfg->Nx);
	for(i=ix_min; i<ix_max; i++) {
		x=cfg->x0+cfg->dx*i; 
		fprintf(file, "\n%6.2f\t%6.1f\t%6.2f", 
			x/1000, (gravanom[i]-gravanom[ix_max-2]) * 1e5, geoidanom[i]-geoidanom[ix_max-2]);
	}
	fprintf(file, "\n");
	fclose(file);
	return(1);
}


int write_file_maxmompoint (ModelConfig *cfg, ModelContext *ctx)
{
	int	i ;
	FILE 	*file ;

	Write_Open_Filename_Return (".ysen", "wt", isost_model<3);

	fprintf(file, 	"#MMP_x= %.2f km\n"
			"#z(km)  Compress(MPa) Extens(MPa)  Temp(C)  Stress [MPa]\n",
			(imomentmax*cfg->dx+cfg->x0)/1e3);
	for (i=0; i<cfg->Nz; i++) {
		fprintf(file, "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", 
			i*cfg->dz/1e3, 
			yieldcompres[imomentmax][i]/1e6, yieldextens[imomentmax][i]/1e6, 
			(switch_YSE_file)? NO_DATA:Temperature[imomentmax][i], 
			stress[imomentmax][i]/1e6);
	}
	fclose(file);
	return(1);
}



int write_file_stress (ModelConfig *cfg, ModelContext *ctx)
{
	int	ix, iz, ix_min, ix_max ;
	FILE 	*file ;

	Write_Open_Filename_Return (".strs", "wt", isost_model<3);

	fprintf(file, "#2D Stress grid (x-z) distribution.\n") ;
	fprintf(file, "#x(km)\tz(km)\tstress[MPa]\ttemp[C]\tyieldcomp\tyieldext[MPa]\n");
	/*Write stresses in the grid file*/
	ix_min = MAX_2((cfg->xmin-cfg->x0-.1*cfg->dx)/cfg->dx, 0) ;	ix_max = MIN_2(floor((cfg->xmax-cfg->x0+.1*cfg->dx)/cfg->dx) + 2, cfg->Nx);
	for (ix=ix_min; ix<ix_max; ix++) for (iz=0; iz<cfg->Nz; iz++) {
		fprintf (file, "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", 
			(ix*cfg->dx+cfg->x0)/1e3, (iz*cfg->dz+w[ix]-zini)/1e3, stress[ix][iz]/1e6, (switch_YSE_file)? NO_DATA:Temperature[ix][iz], 
			yieldcompres[ix][iz]/1e6, yieldextens[ix][iz]/1e6);
	}
	fclose(file);
	return (1);
}



int write_file_Te (ModelConfig *cfg, ModelContext *ctx) 
{
	int	i;
	float	x;
	FILE 	*file ;

	Write_Open_Filename_Return (".eeth", "wt", isost_model<3);

	/*Write a file with the calculated EET*/
	fprintf(file, "# %s: Calculated EET.\n", projectname) ;
	for (i=0; i<cfg->Nx; i++) {
		x=cfg->x0+cfg->dx*i;
		if (x>cfg->xmin-cfg->dx && x<cfg->xmax+cfg->dx) {
 			fprintf(file, "%.2f\t%.2f\n", (i*cfg->dx+cfg->x0)/1e3, Te[i]);
		}
	}
	fclose(file);
	return(1);
}



int write_file_Temperature (ModelConfig *cfg, ModelContext *ctx)
{
	int	ix, iz, ix_min, ix_max ;
	FILE 	*file ;

	Write_Open_Filename_Return (".temp", "wt", isost_model<3 || switch_YSE_file);

	fprintf(file, "# 2D Temperature grid (x-z) distribution.\n") ;
	fprintf(file, "# x (km)\t z (km)\t temperature (�C)\n") ;
	/*Write temperatures in the grid file*/
	ix_min = MAX_2((cfg->xmin-cfg->x0-.1*cfg->dx)/cfg->dx, 0) ;	ix_max = MIN_2(floor((cfg->xmax-cfg->x0+.1*cfg->dx)/cfg->dx) + 2, cfg->Nx);
	for (ix=ix_min; ix<ix_max; ix++) for (iz=0; iz<cfg->Nz; iz++) {
		fprintf (file, "%.2f\t%.2f\t%.2f\n", 
			(ix*cfg->dx+cfg->x0)/1e3, (ctx->topo[ix]-iz*cfg->dz)/1e3, Temperature[ix][iz]);
	}
	fclose(file);
	return (1);
}



int write_file_Temperature_initial (ModelConfig *cfg, ModelContext *ctx)
{
	int	ix, iz, ix_min, ix_max ;
	FILE 	*file ;

	Write_Open_Filename_Return (".tempini", "wt", isost_model<3 || verbose_level<1);

	calculate_topo(cfg, ctx, ctx->topo);
	ix_min = MAX_2((cfg->xmin-cfg->x0-.1*cfg->dx)/cfg->dx, 0) ;	ix_max = MIN_2(floor((cfg->xmax-cfg->x0+.1*cfg->dx)/cfg->dx) + 2, cfg->Nx);
	fprintf(file, "# 2D Temperature grid (x-z) distribution.\n") ;
	fprintf(file, "# x (km)\t z (km)\t temperature (�C)\n") ;
	/*Write temperatures in the grid file*/
	for (ix=ix_min; ix<ix_max; ix++) for (iz=0; iz<cfg->Nz; iz++) {
		fprintf (file, "%.2f\t%.2f\t%.2f\n", 
			(ix*cfg->dx+cfg->x0)/1e3, (ctx->topo[ix]-iz*cfg->dz)/1e3, Temperature[ix][iz]);
	}
	fclose(file);
	return (1);
}



int write_file_time (ModelConfig *cfg, ModelContext *ctx, float *w, float *topo)
{
	/*
	  WRITES deflection and elevation along time file
	*/

	int 	i;
	FILE 	*file;
	char 	filename[MAXLENLINE], filename1[MAXLENLINE], filename2[MAXLENLINE], 
		command[MAXLENLINE];
	bool	return_cond;
	float	youngest_age=-1e16;

	for (i=0; i<ctx->numBlocks; i++) youngest_age = MAX_2(Blocks[i].age, youngest_age);
	return_cond = (!switch_write_file && (ctx->Timefinal-ctx->Time) >= ctx->dt) 
		|| !isost_model 
		|| (((ctx->Time-last_time_file_time) < dt_record*.9999 || (!dt_record && youngest_age!=ctx->Time)) && (ctx->Timefinal-ctx->Time) >= ctx->dt) 
		|| (ctx->Time-last_time_file_time) == 0;

	calculate_topo(cfg, ctx, topo);
	if (!switch_write_file) nwrotenfiles=0;
	if (nwrotenfiles==0) {
		Write_Open_Filename_Return (".xzt", "wt", return_cond);
		fprintf(file, "#Time: \t%.3f My\n#x(km)\tw(m)\th(m)", ctx->Time/Matosec);
	}
	else {
		if (return_cond) return (0);
		sprintf(filename, "%s.xzt", projectname);
		if (cfg->verbose_level>=3) fprintf(stdout, "\nInfo: Writing file '%s' (%d times).", filename, nwrotenfiles+1);
		sprintf(filename1, "%s.aux1.xzt.tao.tmp", projectname);
		sprintf(filename2, "%s.aux2.xzt.tao.tmp", projectname);
		if ((file = fopen(filename1, "wt")) == NULL) {
			PRINT_ERROR("Cannot open auxiliar output file %s.\n", filename1);
			return (0);
		}
		fprintf(file,      "%.3f My\nw(m)\th(m)", ctx->Time/Matosec);
	}

	for(i=0; i<cfg->Nx; i++) {
		float x;
		x=cfg->x0+cfg->dx*i; 
		if (x > cfg->xmin-cfg->dx && x < (cfg->xmax+cfg->dx)) {
			if (nwrotenfiles==0)	fprintf(file, "\n%1.2f\t%1.3f\t%1.3f",
							x/1000, w[i], topo[i] );
			else			fprintf(file, "\n%1.3f\t%1.3f", 
							w[i], topo[i] );
		}
	}
	fprintf(file, "\n");
	fclose(file);
	if (nwrotenfiles > 0) {
		/*Paste the new columns into the .xzt file.*/
		sprintf(command, 
			"paste %s %s > %s", filename, filename1, filename2);
		system(command);
		rename(filename2, filename); remove(filename1);
	}

	switch_write_file_Blocks=true;
	nwrotenfiles++;
	last_time_file_time = ctx->Time;

	return(1);
}


int write_file_Blocks(ModelConfig *cfg, ModelContext *ctx)
{
	FILE 	*file ;

	/*PRINTS A FILE WITH HORIZON ALTITUDES IN COLUMNS*/

	Write_Open_Filename_Return (".pfl", "wt", !switch_write_file_Blocks);

	fprintf(file, "#tAo output file of project '%s'. t= %.2f My", projectname, ctx->Time/Matosec);
	fprintf(file, "\n#Densities:\t%8.0f", cfg->denscrust);
	for (int i=0; i<ctx->numBlocks; i++) fprintf(file, "\t%8.0f", Blocks[i].density);
	if (cfg->erosed_model>=2) fprintf(file, "\t%8.0f", cfg->denswater);
	fprintf(file, "\n#x(km),Ages->\t%8.2f", ctx->Timeini/Matosec);
	for (int i=0; i<ctx->numBlocks; i++) fprintf(file, "\t%8.2f", Blocks[i].age/Matosec);
	if (cfg->erosed_model>=2) fprintf(file, "\t   water");
	for (int i=0; i<cfg->Nx; i++) {
		float x;
		x=cfg->x0+cfg->dx*i; 
		if (x > cfg->xmin-cfg->dx && x < cfg->xmax+cfg->dx) {
			float thickness_above=0, top_block;
			fprintf(file, "\n%8.2f", x/1000);
			for (int i_Block=0; i_Block<ctx->numBlocks; i_Block++) 
				thickness_above += Blocks[i_Block].thick[i];
			top_block = Blocks_base[i]-w[i];
			fprintf(file, "\t%8.1f",  top_block);
			for (int i_Block=0; i_Block<ctx->numBlocks; i_Block++) {
				thickness_above -= Blocks[i_Block].thick[i];
				top_block += Blocks[i_Block].thick[i];
				if (Blocks[i_Block].density==cfg->denssedim) top_block -= compaction(cfg->sed_porosity, compact_depth, thickness_above, thickness_above+Blocks[i_Block].thick[i]);
				fprintf(file, "\t%8.1f",  top_block);
			}
			if (cfg->erosed_model>=2) {
				float top_water=top_block;
				if (drainage[i].lake) top_water = Lake[drainage[i].lake].alt;
				fprintf(file, "\t%8.1f", top_water);
			}
		}
	}
	fclose(file);
	return(1);
}



int write_file_resume(ModelConfig *cfg, ModelContext *ctx)
{
	int 	i, j, end_check=12345;
	FILE 	*file ;

	/*
	  PRINTS A BINARY FILE WITH ALL THE INFORMATION 
	  REQUIRED TO RESTART THE PROGRAM
	*/

	Write_Open_Filename_Return (".all", "wt", !switch_write_file_Blocks);

	/*Defined in universal.h:*/
	fwrite(&cfg->Nx, 		sizeof(int),		1, 	file);
	fwrite(&cfg->Nz, 		sizeof(int),		1, 	file);
	fwrite(&cfg->verbose_level, 	sizeof(int),		1, 	file);

	fwrite(&cfg->dx, 		sizeof(float),		1, 	file);
	fwrite(&cfg->dz, 		sizeof(float),		1, 	file);

	fwrite(version,		sizeof(char),		LENGTHVERS, 	file);
	fwrite(version_input,	sizeof(char),		LENGTHVERS, 	file);

	fwrite(&switch_geograph_coor, sizeof(bool),	1, 	file);
	fwrite(&switch_ps, 	sizeof(bool),		1, 	file);
	fwrite(&switch_write_file, sizeof(bool),		1, 	file);


	/*Defined in geomodel.h:*/
	fwrite(&grav_anom_type, 	sizeof(int),		1, 	file);
	fwrite(&isost_model, 	sizeof(int),		1, 	file);

	fwrite(&Te_default, 	sizeof(float),		1, 	file);
	fwrite(&crust_thick_default, sizeof(float),	1, 	file);
	fwrite(&upper_crust_thick_default, sizeof(float),1, 	file);
	fwrite(&cfg->densasthen, 	sizeof(float),		1, 	file);
	fwrite(&cfg->densmantle, 	sizeof(float),		1, 	file);
	fwrite(&cfg->denssedim, 	sizeof(float),		1, 	file);
	fwrite(&cfg->denscrust, 	sizeof(float),		1, 	file);
	fwrite(&cfg->densenv, 	sizeof(float),		1, 	file);
	fwrite(&cfg->densinfill, 	sizeof(float),		1, 	file);
	fwrite(&ctx->sea_level, 	sizeof(float),		1, 	file);
	fwrite(&temp_sea_level, 	sizeof(float),		1, 	file);
	fwrite(&ctx->Time, 		sizeof(float),		1, 	file);
	fwrite(&ctx->Timefinal, 	sizeof(float),		1, 	file);
	fwrite(&ctx->Timeini, 	sizeof(float),		1, 	file);
	fwrite(&ctx->dt, 		sizeof(float),		1, 	file);
	fwrite(&ctx->dt_eros, 		sizeof(float),		1, 	file);
	fwrite(&tau, 		sizeof(float),		1, 	file);

	fwrite(projectname, 	sizeof(char),	MAXLENFILE, 	file);
	fwrite(gif_geom, 	sizeof(char),	MAXLENLINE, 	file);

	fwrite(&water_load, 	sizeof(bool),		1, 	file);


	/*Defined in tao+tisc.h:*/
	fwrite(&boundary_conds, 	sizeof(int),		1, 	file);
	fwrite(&nloads, 		sizeof(int),		1, 	file);
	fwrite(&nmax_input_points, sizeof(int),		1, 	file);
	fwrite(&n_sea_level_input_points, sizeof(int),	1, 	file);
	fwrite(&n_eros_level_input_points, sizeof(int),	1, 	file);
	fwrite(&n_record_times, 	sizeof(int),		1, 	file);
	fwrite(&i_first_Block_load, sizeof(int),	1, 	file);
	fwrite(&i_Block_insert, sizeof(int),		1, 	file);
	fwrite(&ctx->numBlocks, 	sizeof(int),		1, 	file);
	fwrite(&nwrotenfiles, 	sizeof(int),		1, 	file);
	fwrite(&run_type, 	sizeof(int),		1, 	file);

	fwrite(&zini, 		sizeof(float),		1, 	file);
	fwrite(&dt_record, 	sizeof(float),		1, 	file);
	fwrite(&cfg->sed_porosity, 	sizeof(float),		1, 	file);
	fwrite(&compact_depth, 	sizeof(float),		1, 	file);
	fwrite(&Kerosdif, 	sizeof(float),		1, 	file);
	fwrite(&last_time_file_time, 	sizeof(float),		1, 	file);
	fwrite(&random_topo, 	sizeof(float),		1, 	file);

	fwrite(&switch_file_out, 	sizeof(bool),		1, 	file);
	fwrite(&switch_gradual, 	sizeof(bool),		1, 	file);
	fwrite(&switch_insert_load, 	sizeof(bool),		1, 	file);
	fwrite(&switch_topoest, 		sizeof(bool),		1, 	file);
	fwrite(&switch_write_file_Blocks, sizeof(bool),		1, 	file);
	fwrite(&deform_sed, sizeof(bool),		1, 	file);

	/*Defined in tao.h:*/
	fwrite(&imomentmax, 	sizeof(int),		1, 	file);
	fwrite(&nx_temp_input, 	sizeof(int),		1, 	file);
	fwrite(&nbasins, 	sizeof(int),		1, 	file);
	fwrite(&ctx->nlakes,	 	sizeof(int),		1, 	file);
	fwrite(&n_image, 	sizeof(int),		1, 	file);
	fwrite(&cfg->hydro_model,	sizeof(int),		1, 	file);
	fwrite(&cfg->erosed_model,	sizeof(int),		1, 	file);
	fwrite(&eros_bound_cond,	sizeof(char),		2, 	file);

	fwrite(&cfg->x0, 		sizeof(float),		1, 	file);
	fwrite(&cfg->xf, 		sizeof(float),		1, 	file);
	fwrite(&cfg->xmin, 		sizeof(float),		1, 	file);
	fwrite(&cfg->xmax, 		sizeof(float),		1, 	file);
	fwrite(&zmin, 		sizeof(float),		1, 	file);
	fwrite(&zmax, 		sizeof(float),		1, 	file);
	fwrite(&horz_force, 	sizeof(float),		1, 	file);
	fwrite(&vert_force, 	sizeof(float),		1, 	file);
	fwrite(&appmoment, 	sizeof(float),		1, 	file);
	fwrite(&Keroseol, 	sizeof(float),		1, 	file);
	fwrite(&Ksedim, 		sizeof(float),		1, 	file);
	fwrite(&critical_slope, 	sizeof(float),		1, 	file);
	fwrite(&K_river_cap, 	sizeof(float),		1, 	file);
	fwrite(&erodibility, 	sizeof(float),		1, 	file);
	fwrite(&erodibility_sed, sizeof(float),		1, 	file);
	fwrite(&l_fluv_sedim, 	sizeof(float),		1, 	file);
	fwrite(&lost_rate, 	sizeof(float),		1, 	file);
	fwrite(&evaporation_ct, 	sizeof(float),		1, 	file);
	fwrite(&cfg->riverbasinwidth, 	sizeof(float),		1, 	file);
	fwrite(&rain, 		sizeof(float),		1, 	file);
	fwrite(&Krain,	 	sizeof(float),		1, 	file);
	fwrite(&CXrain, 		sizeof(float),		1, 	file);

	fwrite(&switch_strs_history, 	sizeof(bool),		1, 	file);
	fwrite(&switch_YSE_file, 	sizeof(bool),		1, 	file);


	/*Arrays:*/
	fwrite(w, 		sizeof(float),	cfg->Nx, 	file);
	fwrite(D, 		sizeof(float),	cfg->Nx, 	file);
	fwrite(q, 		sizeof(float),	cfg->Nx, 	file);
	fwrite(Dw, 		sizeof(float),	cfg->Nx, 	file);
	fwrite(Dq, 		sizeof(float),	cfg->Nx, 	file);
	fwrite(h_water, 		sizeof(float),	cfg->Nx, 	file);
	fwrite(h_last_unit, 	sizeof(float),	cfg->Nx, 	file);
	fwrite(Te, 		sizeof(float),	cfg->Nx, 	file);
	fwrite(crust_thick, 	sizeof(float),	cfg->Nx, 	file);
	fwrite(upper_crust_thick,sizeof(float),	cfg->Nx, 	file);
	fwrite(ctx->topo, 		sizeof(float),	cfg->Nx, 	file);
	fwrite(Blocks_base, 	sizeof(float),	cfg->Nx, 	file);
	fwrite(horiz_record_time, sizeof(float), n_record_times, file);

	for (i=0; i<n_sea_level_input_points; i++) 
		fwrite(var_sea_level[i], sizeof(float), 2, file);
	for (i=0; i<n_eros_level_input_points; i++) 
		fwrite(var_eros_level[i], sizeof(float), 2, file);

	if (isost_model>=3 && !switch_YSE_file) {
		if (fwrite(Temperature[0], sizeof(float), cfg->Nx * cfg->Nz, file) != (size_t)(cfg->Nx * cfg->Nz)) { PRINT_ERROR("Failed to write Temperature array."); return 0; }
	}
	if (isost_model>=3) {
		if (fwrite(stress[0], sizeof(float), cfg->Nx * cfg->Nz, file) != (size_t)(cfg->Nx * cfg->Nz)) { PRINT_ERROR("Failed to write stress array."); return 0; }
		if (fwrite(yieldcompres[0], sizeof(float), cfg->Nx * cfg->Nz, file) != (size_t)(cfg->Nx * cfg->Nz)) { PRINT_ERROR("Failed to write yieldcompres array."); return 0; }
		if (fwrite(yieldextens[0], sizeof(float), cfg->Nx * cfg->Nz, file) != (size_t)(cfg->Nx * cfg->Nz)) { PRINT_ERROR("Failed to write yieldextens array."); return 0; }
	}

	fwrite(Blocks, 		sizeof(struct BLOCK_1D),	ctx->numBlocks, file);
	for (i=0; i<ctx->numBlocks; i++) {
		fwrite(Blocks[i].thick, 		sizeof(float),	cfg->Nx, 	file);
	}
	for (i=0; i<ctx->numBlocks; i++) {
	    if (Blocks[i].type == 'S') {
		fwrite(Blocks[i].detr_ratio, 	sizeof(float),	cfg->Nx, 	file);
		fwrite(Blocks[i].detr_grsize, 	sizeof(float),	cfg->Nx, 	file);
	    }
	}

	if (cfg->erosed_model) {
		fwrite(eros_now, 	sizeof(float),	cfg->Nx, 	file);
		fwrite(total_erosion, 	sizeof(float),	cfg->Nx, 	file);
	}
	if (cfg->hydro_model) {
		fwrite(precipitation, 	sizeof(float),	cfg->Nx, 	file);
		fwrite(evaporation, 	sizeof(float),	cfg->Nx, 	file);
		fwrite(Lake, sizeof(struct LAKE_INFO_1D), ctx->nlakes+1, file);
		for (j=1; j<=ctx->nlakes; j++) {
			fwrite(Lake[j].cell, sizeof(int), Lake[j].n, file);
			fwrite(Lake[j].sd, sizeof(int), Lake[j].n_sd, file);
		}
	}

	fwrite(&end_check,	sizeof(int),		1, 	file);

	fclose(file);
	return(1);
}
