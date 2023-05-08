/*
INPUT/OUTPUT  SUBROUTINES  FOR  tisc.c
Daniel Garcia-Castellanos
*/


int find_up_river (int row, int col, int *level, int *count, float *length, float *chi, FILE *file, BOOL **done, float ref_discharge);





/*******************************  INPUT  **************************************/

int reformat_file_thin_sheet_BC(char *tmpTSBCfilename)
{
	/*
	  READS THE THIN SHEET BOUNDARY CONDITIONS FILE *.TSBC AND TRANSLATES 
	  IT INTO A TEMPORARY FILE READABLE BY thin_sheet.f, which reads lines 
	  of this type:
	  	col, row, typeBC, velx, vely 
	  Different to TISC, here row=0 is at the smallest y (y=ymin, the southernmost row) and row=Ny-1 is at ymax.
	*/
	int	i, j;
	FILE 	*file, *filetmp;
	char	boundary=0, *lin, auxstr[MAXLENFILE];
	float 	x,u1,u2, x0,u01,u02;

	Read_Open_Filename_Return(".TSBC", "rt", "Thin sheet boundary conditions")
	if ((filetmp=fopen(tmpTSBCfilename,"wt"))==NULL) {PRINT_WARNING("Could not open output file '%s'.", tmpTSBCfilename); return 0;}
	PRINT_INFO("Writing intermediate file '%s' for thin_sheet.",tmpTSBCfilename);
	fprintf(filetmp, "#temporary file for thin sheet boundary conditions (produced by TISC).");

	while (1) {
	    char  word[MAXLENLINE];
	    int   num_fields=0, type;
	    lin = fgets(auxstr, MAXLENLINE-1, file);
	    if (lin) num_fields = sscanf(lin, "%s", word);
	    if (num_fields || !lin) {
	    	if (word[0]=='N' || word[0]=='S' || word[0]=='E' || word[0]=='W' || !lin) {
		    switch (boundary) {
		      float xi;
		      case 'N':
		      case 'S':
		    	  for (i=0; i<Nx; i++) {
		    	      xi = xmin + i*dx;
		    	      if (xi>x)
		    		  fprintf(filetmp, "\n%d\t%d\t%d\t%E\t%E",
		    		      i, (boundary=='S')? 0:Ny-1, type, u1, u2);
		    	  }
		    	  break;
		      case 'E':
		      case 'W':
		    	  for (i=1; i<Ny-1; i++) {
		    	      xi = ymin + i*dy;
		    	      if (xi>x)
		    		  fprintf(filetmp, "\n%d\t%d\t%d\t%E\t%E",
		    		      (boundary=='W')? 0:Nx-1, i, type, u1, u2);
		    	  }
		    	  break;
		    }
		    if (!lin) break;
		    boundary = word[0]; type = atoi(&word[1]);
		    x0=u01=u02=SIGNAL;
		}
	    }
	    if ((num_fields = sscanf(lin, "%f %f %f", &x, &u1, &u2)) == 3) {
	      if (type==1) {u1*=1e3/Matosec; u2*=1e3/Matosec;};
	      switch (boundary) {
	        float xi;
		case 'N':
		case 'S':
		    for (i=0; i<Nx; i++) {
			xi = xmin + i*dx;
			if (x0==SIGNAL) { 
			    if (xi<=x) 
			    	fprintf(filetmp, "\n%d\t%d\t%d\t%E\t%E", 
			    	    i, (boundary=='S')? 0:Ny-1, type, u1, u2);
			}
			else {
			    if (xi>x0 && xi<=x) 
			    	fprintf(filetmp, "\n%d\t%d\t%d\t%E\t%E", 
			    	    i, (boundary=='S')? 0:Ny-1, type, LININTERP(xi,x0,x,u01,u1), LININTERP(xi,x0,x,u02,u2));
			}
		    }
		    break;
		case 'E':
		case 'W':
		    for (i=1; i<Ny-1; i++) {
		    	xi = ymin + i*dy;
			if (x0==SIGNAL) {
			    if (xi<=x) 
			    	fprintf(filetmp, "\n%d\t%d\t%d\t%E\t%E", 
			    	    (boundary=='W')? 0:Nx-1, i, type, u1, u2);
			}
			else {
			    if (xi>x0 && xi<=x) 
			    	fprintf(filetmp, "\n%d\t%d\t%d\t%E\t%E", 
			    	    (boundary=='W')? 0:Nx-1, i, type, LININTERP(xi,x0,x,u01,u1), LININTERP(xi,x0,x,u02,u2));
			}
		    }
		    break;
	      }
	      x0=x; u01=u1; u02=u2; 
	    }
	}

    	fprintf(filetmp, "\n"); 
	PRINT_INFO("Thin sheet BC in '%s.TSBC'", projectname);

	fclose(file); fclose(filetmp);

	return 1;
}




int read_file_horiz_record_time()
{
	/*
	  Reads file containing the ages of horizons to be recorded: 'projectname.REC'
	*/

	int 	i, nmax_input_points=5000;
	FILE 	*file;
	float	*aux1;

	Read_Open_Filename_Return(".REC", "rt", "Horizon recording times")

	n_record_times=0;
	aux1 = calloc(nmax_input_points, sizeof(float));
	for (;;) {
		TAKE_LINE_1(aux1[n_record_times]);
		n_record_times++;
		if (n_record_times>=nmax_input_points-1 ) {
			PRINT_ERROR("Too many points (>%d) in horizon recording times file.", n_record_times-1);
			break;
		}
	}
	fclose(file);
	PRINT_INFO("%d times read.", n_record_times);

	horiz_record_time = calloc(n_record_times, sizeof(float));
	for (i=0; i<n_record_times; i++) {
		horiz_record_time[i] = aux1[i]*Matosec;
	}
	
	/*desactivate automatic sediment Blocks generation*/
	/*if (n_record_times) dt_record=0;*/
	
	free(aux1); 
	return(1);
}




int read_file_initial_deflection (float **w)
{
	/*
	  READS INITIAL DEFLECTION from 'projectname.WINI'
	*/

	FILE 	*file;
	float	z_default=0;
	int	mode_interp_local=mode_interp;

	Read_Open_Filename_Return(".WINI", "rt", "Initial deflection")

	readinterp2D(file, w, mode_interp_local, z_default, xmin, xmax, ymin, ymax, Nx, Ny);

	fclose(file);
	return (1);
}



int read_file_initial_topo(float **h)
{
	/*
	  READS FILE CONTAINING INTIAL TOPOGRAPHY
	*/

	float	z_default=0;
	int	mode_interp_local=mode_interp;
	FILE 	*file;

	Read_Open_Filename_Return(".ZINI", "rt", "Initial topography")

	readinterp2D(file, h, mode_interp_local, z_default, xmin, xmax, ymin, ymax, Nx, Ny);

	fclose(file);

	return (1) ;
}




int read_file_initial_rivers()
{
	/* 
	  Reads an input file with the rivers that must exist at Tini.
	  Initial topography is adjusted to fit this rivers.
	*/

	FILE 	*file;
	int	num_fields, n_riv_pts=0, n_rivers=0, 
		i, j, ia, ja, ii, jj, iia, jja;
	float	x, y, z, xa, ya, topoant, toporiv=SIGNAL,
		min_descent = .1;

#ifdef SURFACE_TRANSPORT
	Read_Open_Filename_Return(".RIV", "rt", "Initial rivers")

	while (1) {
		char *lin;
		/*Reads a river x,y, point*/
		while (1) {
			char auxstr[MAXLENLINE];
			lin=fgets(auxstr, MAXLENLINE-1, file);
			if (lin==NULL) break;
			if ((num_fields = sscanf(lin, "%f %f %f", &x, &y, &z)) >= 2) {
				if (num_fields>2) toporiv=z; else toporiv=SIGNAL;
				if (n_riv_pts==0) {
					n_riv_pts = 1;
					xa=x; ya=y;
				}
				else if (sqrt((x-xa)*(x-xa)+(y-ya)*(y-ya)) > dxy) break;
			}
			if (lin[0] == '>') {n_riv_pts=0;}
		}
		if (lin==NULL) break;
		n_riv_pts++;
		if (n_riv_pts==2) n_rivers++;

		ja = floor((xa-xmin)/dx + .49999999);
		j  = floor((x -xmin)/dx + .49999999);
		ia = floor((ymax-ya)/dy + .49999999);
		i  = floor((ymax-y )/dy + .49999999);
		iia=ia; jja=ja;
		if (IN_DOMAIN(ia,ja)) topoant=topo[ia][ja];
		if (fabs((float) (j-ja)) >= fabs((float) (i-ia))) {
		    if (j>ja)  for (jj=ja+1; jj<=j; jj++) {
			    ii = ia+(i-ia)*(jj-ja)/(j-ja);
			    if (IN_DOMAIN(ii,jj) && IN_DOMAIN(ia,ja)) {
				topo[ii][jj] =
					MIN_2(topoant-min_descent, topo[ii][jj]);
				Damn_River_Node (iia, jja, ii, jj);
				topoant = topo[ii][jj];
				iia=ii; jja=jj;
			    }
		    }
		    if (j<ja)  for (jj=ja-1; jj>=j; jj--) {
			    ii = ia+(i-ia)*(jj-ja)/(j-ja);
			    if (IN_DOMAIN(ii,jj) && IN_DOMAIN(ia,ja)) {
				topo[ii][jj] =
					MIN_2(topoant-min_descent, topo[ii][jj]);
				Damn_River_Node (iia, jja, ii, jj);
				topoant = topo[ii][jj];
				iia=ii; jja=jj;
			    }
		    }
		}
		else {
		    if (i>ia)  for (ii=ia+1; ii<=i; ii++) {
			    jj = ja+(j-ja)*(ii-ia)/(i-ia);
			    if (IN_DOMAIN(ii,jj) && IN_DOMAIN(ia,ja)) {
				topo[ii][jj] =
					MIN_2(topoant-min_descent, topo[ii][jj]);
				Damn_River_Node (iia, jja, ii, jj);
				topoant = topo[ii][jj];
				iia=ii; jja=jj;
			    }
		    }
		    if (i<ia)  for (ii=ia-1; ii>=i; ii--) {
			    jj = ja+(j-ja)*(ii-ia)/(i-ia);
			    if (IN_DOMAIN(ii,jj) && IN_DOMAIN(ia,ja)) {
				topo[ii][jj] =
					MIN_2(topoant-min_descent, topo[ii][jj]);
				Damn_River_Node (iia, jja, ii, jj);
				topoant = topo[ii][jj];
				iia=ii; jja=jj;
			    }
		    }
		}
		xa=x;
		ya=y;
	}

	PRINT_INFO("%d rivers (last defined by %d x-y points)", n_rivers, n_riv_pts);

	fclose(file);
#endif

	return 1;
}



int read_file_parameters (int show, int reformat) 
{
	int	nread=0, nparams=0, nline=0, verbose_level_ant=verbose_level;
	char 	*lineptr=NULL, str1[MAXLENLINE]="", str2[MAXLENLINE]="", 
		line[MAXLENLINE+200], PRMfilename[MAXLENFILE];
	FILE 	*file;
	BOOL	switch_matched_vers=NO;

	/*
	READ THE PARAMETERS FILE NAMED  'projectname.PRM'
	You have an explanation of these parameters in the example (and default) 
	parameters file (doc/template.PRM), in function match_parameter you have 
	the name of the variable related to each parameter. The use of most 
	variables is described in the include files '*.h'.
	*/

	sprintf(PRMfilename, "%s.PRM", projectname);
	if (show && verbose_level>=3) fprintf(stdout, "\nCurrent TISC project: %s", projectname);
	if ((file = fopen(PRMfilename, "rt")) == NULL) {
		PRINT_INFO("Can't open parameters file '%s'.\n", PRMfilename);
		return(0);
	}

	if (show) fprintf(stdout, "\nParameters at '%s'.", PRMfilename);
	if (reformat==1) fprintf(stdout, "\n\n");
	while ((lineptr=fgets(line, MAXLENLINE+200-1, file)) != NULL) {
		int status;
		status=0;
		nread=sscanf(lineptr, "%s %s", str1, str2);
		if (nread >= 2) {
			if (!strcmp(str1, "version")) {
				strcpy(version_input, str2);
				if (!strcmp(version_input, version)) {
					switch_matched_vers = YES;
					nparams++;
				}
				if ((show && verbose_level>=3) || (reformat && verbose_level>=3)) fprintf(stderr, "\nInput file version: %s\n", version_input);
				if (reformat) {fprintf(stdout, "version\t%s\n", version); }
			}
			status=match_parameter(str1, str2, show, reformat, line);
			nparams += status; 
			if (!strcmp(str1, "version")) status=1;
		}
		/*If no parameter matched then just reproduce the entire line*/
		if (reformat==1 && !status) fprintf(stdout, "%s", line);
		nline++;
	}
	if (verbose_level_ant>=2) fprintf(stdout, " (%d parameters matched)", nparams);
	if (!switch_matched_vers) {
		if (verbose_level_ant>=2) 
		    fprintf(stderr,
			"\nInfo: Present version name '%s' not matched in PRM file. See 'tisc/doc/template.PRM'.", version);
		if (nparams<2) {
		    PRINT_ERROR("\aWrong format in parameters file '%s'. Only %d parameters were matched. Current version is '%s'. "
			"\nSee example file 'tisc/doc/template.PRM'.\nEND.\n", 
			PRMfilename, nparams, version);
		    exit(0); 
		}
	}

	fclose(file);
	return(1);
}





int match_parameter (char *str1, char *str2, int show, int replace, char *line)
{
	BOOL switch_debug=NO;
	int nparams=0;

	Match_Param_Replace_int ( "Nx", 	Nx,  	0 )
	Match_Param_Replace_int ( "Ny", 	Ny,  	0 )
	Match_Param_Replace_flt ( "xmin", 	xmin,  	0 )
	Match_Param_Replace_flt ( "xmax", 	xmax,  	0 )
	Match_Param_Replace_flt ( "ymin", 	ymin,  	0 )
	Match_Param_Replace_flt ( "ymax", 	ymax,  	0 )
	Match_Param_Replace_flt ( "Te",		Te_default,  	0 )
	Match_Param_Replace_flt ( "zini",	zini,  	0 )
	Match_Param_Replace_flt ( "random_topo",	random_topo,  	0 )
	Match_Param_Replace_int ( "mode_interp",	mode_interp,  	0 )
	Match_Param_Replace_flt ( "densasthen",	densasthen,  	0 )
	Match_Param_Replace_flt ( "densmantle",	densmantle,  	0 )
	Match_Param_Replace_flt ( "denscrust",	denscrust,  	0 )
	Match_Param_Replace_flt ( "densinfill",	densinfill,  	0 )
	Match_Param_Replace_flt ( "denssedim",	denssedim,  	0 )
	Match_Param_Replace_flt ( "densenv",	densenv,  	0 )
	Match_Param_Replace_flt ( "sed_porosity",	sed_porosity,  	0 )
	Match_Param_Replace_flt ( "compact_depth",	compact_depth, 	0 )
	Match_Param_Replace_chr ( "boundary_conds",	boundary_conds,  	0 )
	Match_Param_Replace_flt ( "Px", 	Px,  	0 )
	Match_Param_Replace_flt ( "Py", 	Py,  	0 )
	Match_Param_Replace_flt ( "Pxy",	Pxy,  	0 )
	Match_Param_Replace_int ( "hydro_model",	hydro_model,  	0 )
	Match_Param_Replace_flt ( "rain",	rain,  	0 )
	Match_Param_Replace_flt ( "Krain",	Krain,  	0 )
	Match_Param_Replace_flt ( "relhumid",	relative_humidity,  	0 )
	Match_Param_Replace_flt ( "windazimut",	windazimut,  	0 )
	Match_Param_Replace_flt ( "CXrain",	CXrain,  	0 )
	Match_Param_Replace_flt ( "CYrain",	CYrain,  	0 )
	Match_Param_Replace_flt ( "evaporation",	evaporation_ct,  	0 )
	Match_Param_Replace_flt ( "lost_rate",	lost_rate,  	0 )
	Match_Param_Replace_flt ( "permeability",	permeability,  	0 )
	Match_Param_Replace_int ( "erosed_model",	erosed_model,  	0 )
	Match_Param_Replace_flt ( "Kerosdif",	Kerosdif,  	0 )
	Match_Param_Replace_flt ( "Keroseol",	Keroseol, 	0 )
	Match_Param_Replace_flt ( "Ksedim",	Ksedim, 	0 )
	Match_Param_Replace_flt ( "critical_slope",	critical_slope,  	0 )
	Match_Param_Replace_flt ( "K_river_cap",	K_river_cap,  	0 )
	Match_Param_Replace_flt ( "erodibility",	erodibility,  	0 )
	Match_Param_Replace_flt ( "erodibility_sed",	erodibility_sed,  	0 )
	Match_Param_Replace_flt ( "critical_stress",	critical_stress,  	0 )
	Match_Param_Replace_flt ( "l_fluv_sedim",	l_fluv_sedim,  	0 )
	Match_Param_Replace_flt ( "temp_sea_level",	temp_sea_level,  	0 )
	Match_Param_Replace_int ( "deform_sed",  	deform_sed, 	0 )
	Match_Param_Replace_flt ( "K_ice_eros",	K_ice_eros,  	0 )
	Match_Param_Replace_flt ( "dt_ice",	dt_ice,  	0 )
	Match_Param_Replace_int ( "n_ice_flow",	n_ice_flow,  	0 )
	Match_Param_Replace_flt ( "A_ice_rheo",	A_ice_rheo,  	0 )
	Match_Param_Replace_flt ( "A_ice_slide",	A_ice_slide,  	0 )
	Match_Param_Replace_chr ( "eros_bound_cond",	eros_bound_cond,  	0 )
	Match_Param_Replace_flt ( "Timeini",	Timeini,  	0 )
	Match_Param_Replace_flt ( "Timefinal",	Timefinal,  	0 )
	Match_Param_Replace_flt ( "tau",	tau,  	0 )
	Match_Param_Replace_flt ( "dt",		dt,  	0 )
	Match_Param_Replace_flt ( "dt_eros",	dt_eros,  	0 )
	Match_Param_Replace_flt ( "dt_record",	dt_record,  	0 )
	Match_Param_Replace_int ( "isost_model",	isost_model,  	0 )
	Match_Param_Replace_int ( "water_load", 	water_load,  	0 )
	Match_Param_Replace_int ( "switch_topoest", 	switch_topoest,  	0 )
	Match_Param_Replace_int ( "switch_files",	switch_write_file,  	0 )
	Match_Param_Replace_int ( "switch_ps",  	switch_ps,  	0 )
	Match_Param_Replace_int ( "verbose_level",	verbose_level,  	0 )

	/*Old versions:*/
	Match_Param_Replace_flt ( "erodability",	erodibility,  	1 )
	Match_Param_Replace_flt ( "erodability_sed",	erodibility_sed,  	1 )
	Match_Param_Replace_int ( "switch_verbose",	verbose_level,  	1 )
	Match_Param_Replace_int ( "switch_debug",	switch_debug,  	1 )
	if (switch_debug) verbose_level=3;
	Match_Param_Replace_flt ( "alt0",	zini,  	1 )
	Match_Param_Replace_int ( "lith_type",	isost_model,  	1 )
	Match_Param_Replace_int ( "erosed_type",	erosed_model,  	1 )
	Match_Param_Replace_int ( "switch_hydro",	hydro_model,  	1 )
	Match_Param_Replace_flt ( "leng_fluv_eros",	erodibility,  	1 )
	Match_Param_Replace_flt ( "leng_fluv_sedim",	l_fluv_sedim,  	1 )
	Match_Param_Replace_int ( "switch_erosed",	erosed_model,  	1 )
	Match_Param_Replace_int ( "switch_sea", 	water_load,  	1 )
	Match_Param_Replace_flt ( "l_fluv_eros",	erodibility,  	1 )
	Match_Param_Replace_flt ( "l_fluv_eros_sed",	erodibility_sed,  	1 )
	Match_Param_Replace_flt ( "dtmemounit",	dt_record,  	1 )

	return (nparams);
}




int read_file_rain(float **precipitation_file)
{
	/*
	  READ THE RUNOFF FILE *.RAIN
	*/

	FILE 	*file;
	float	total_rain_file=0;
	int	i, j, mode_interp_local=mode_interp;

#ifdef SURFACE_TRANSPORT
	if (!hydro_model) return(0);

	for (i=0;i<Ny;i++) for (j=0; j<Nx; j++) precipitation_file[i][j] = -9999;

	Read_Open_Filename_Return(".RAIN", "rt", "Precipitation distribution")

	readinterp2D(file, precipitation_file, mode_interp_local, -9999, xmin, xmax, ymin, ymax, Nx, Ny);

	fclose(file);

	for (i=0;i<Ny;i++) for (j=0; j<Nx; j++) {
		precipitation_file[i][j] /= (secsperyr*1e3);
		if (precipitation_file[i][j]>0) total_rain_file += precipitation_file[i][j]*dx*dy;
	}
	PRINT_INFO("Total precipitation in file: %.2f m3/s", total_rain_file);
#endif

	return (1) ;
}




int read_file_resume(char *filename)
{
	/*
	  PRINTS A BINARY FILE WITH ALL THE INFORMATION 
	  REQUIRED TO RESUME A PROJECT RUN
	*/

	int 	i, j, numBlocks_aux, i_Block_insert_aux, run_type_aux, end_check;
	float 	**auxptr1;
	FILE 	*file;
	char  	version_aux[LENGTHVERS], projectname_aux[MAXLENFILE];

	if ((file = fopen(filename,"rb")) == NULL) {PRINT_ERROR("Cannot read ""Resume"" input file '%s'.\n", filename); exit(0);}
	if (verbose_level>=1) fprintf(stdout, "\n""Reading Model-run data at '%s'", filename);

	/*Defined in universal.h:*/
	fread(&Nx, 		sizeof(int),		1, 	file);
	fread(&Ny, 		sizeof(int),		1, 	file);
	fread(&Nz, 		sizeof(int),		1, 	file);
	fread(&verbose_level,	sizeof(int),		1, 	file);

	fread(&dx, 		sizeof(float),		1, 	file);
	fread(&dy, 		sizeof(float),		1, 	file);
	fread(&dz, 		sizeof(float),		1, 	file);

	fread(version_aux,	sizeof(char),		LENGTHVERS, 	file);
	if (strcmp(version, version_aux)) PRINT_WARNING("version in resume file '%s' does not match present TISC version '%s'.", version_aux, version);
	fread(version_input,	sizeof(char),		LENGTHVERS, 	file);

	fread(&switch_ps, 	sizeof(BOOL),		1, 	file);
	fread(&switch_write_file, sizeof(BOOL),		1, 	file);


	/*Defined in geomodel.h:*/
	fread(&grav_anom_type, 	sizeof(int),		1, 	file);
	fread(&isost_model, 	sizeof(int),		1, 	file);
	fread(&water_load, 	sizeof(int),		1, 	file);

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
	fread(&dt_eros, 	sizeof(float),		1, 	file);
	fread(&tau, 		sizeof(float),		1, 	file);

	fread(projectname_aux, 	sizeof(char),	MAXLENFILE, 	file);
	if (strlen(projectname)<1) {
		/*Don't use Warnings cause verbose_level is not set yet*/
		fprintf(stdout, "\nAttention: project name in resume file '%s' will be used.", projectname_aux);
		strcpy(projectname,projectname_aux);
	}
	if (strcmp(projectname, projectname_aux)) {
		/*Don't use Warnings cause verbose_level is not set yet*/
		fprintf(stdout, "\nAttention: project name in resume file '%s' will be ignored, '%s' wil be used instead.", projectname_aux, projectname);
		strcpy(projectname_aux,projectname);
	}

	fread(&switch_geograph_coor, sizeof(BOOL),	1, 	file);


	/*Defined in tao+tisc.h:*/
	fread(&nloads, 		sizeof(int),		1, 	file);
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

	fread(&switch_file_out, 	sizeof(BOOL),		1, 	file);
	fread(&switch_gradual, 	sizeof(BOOL),		1, 	file);
	fread(&switch_topoest, 		sizeof(BOOL),		1, 	file);
	fread(&switch_write_file_Blocks, sizeof(BOOL),		1, 	file);
	fread(&deform_sed, sizeof(BOOL),		1, 	file);


	/*Defined in tisc.h:*/
	fread(&erosed_model, 		sizeof(int),		1, 	file);
	fread(&mode_interp, 	sizeof(int),		1, 	file);
	fread(&nbasins, 	sizeof(int),		1, 	file);
	fread(&nlakes, 		sizeof(int),		1, 	file);
	fread(&n_image, 	sizeof(int),		1, 	file);

	fread(&xmin, 		sizeof(float),		1, 	file);
	fread(&xmax, 		sizeof(float),		1, 	file);
	fread(&ymin, 		sizeof(float),		1, 	file);
	fread(&ymax, 		sizeof(float),		1, 	file);
	fread(&dxy, 		sizeof(float),		1, 	file);
	fread(&Px, 		sizeof(float),		1, 	file);
	fread(&Py, 	sizeof(float),		1, 	file);
	fread(&Pxy, 	sizeof(float),		1, 	file);
	fread(&Keroseol, 	sizeof(float),		1, 	file);
	fread(&Ksedim, 		sizeof(float),		1, 	file);
	fread(&critical_slope, 	sizeof(float),		1, 	file);
	fread(&K_river_cap, 	sizeof(float),		1, 	file);
	fread(&K_ice_eros, 	sizeof(float),		1, 	file);
	fread(&dt_ice, 	sizeof(float),		1, 	file);
	fread(&n_ice_flow, 	sizeof(int),		1, 	file);
	fread(&A_ice_rheo, 	sizeof(float),		1, 	file);
	fread(&A_ice_slide, 	sizeof(float),		1, 	file);
	fread(&erodibility, 	sizeof(float),		1, 	file);
	fread(&erodibility_sed, 	sizeof(float),		1, 	file);
	fread(&critical_stress, 	sizeof(float),		1, 	file);
	fread(&l_fluv_sedim, 	sizeof(float),		1, 	file);
	fread(&lost_rate, 	sizeof(float),		1, 	file);
	fread(&permeability, 	sizeof(float),		1, 	file);
	fread(&evaporation_ct, 	sizeof(float),		1, 	file);
	fread(&rain, 	sizeof(float),		1, 	file);
	fread(&Krain, 	sizeof(float),		1, 	file);
	fread(&relative_humidity, 	sizeof(float),		1, 	file);
	fread(&windazimut, 	sizeof(float),		1, 	file);
	fread(&CXrain, 	sizeof(float),		1, 	file);
	fread(&CYrain, 	sizeof(float),		1, 	file);
	fread(&total_bedrock_eros_mass, 	sizeof(float),		1, 	file);
	fread(&total_sed_mass, 	sizeof(float),		1, 	file);

	fread(&hydro_model, 	sizeof(int),		1, 	file);
	fread(&lake_instant_fill, 	sizeof(int),		1, 	file);

	fread(boundary_conds, 	sizeof(char),	5, 	file);
	fread(eros_bound_cond, 	sizeof(char),	5, 	file);
	fread(&solver_type, 	sizeof(char),	1, 	file);
	fread(gif_geom, 	sizeof(char),	MAXLENLINE, 	file);

	/*Arrays:*/
	Allocate_Memory();
	for (i=0; i<Ny; i++) fread(w[i], sizeof(float), Nx, file);
	for (i=0; i<Ny; i++) fread(D[i], sizeof(float), Nx, file);
	for (i=0; i<Ny; i++) fread(q[i], sizeof(float), Nx, file);
	for (i=0; i<Ny; i++) fread(Dw[i], sizeof(float), Nx, file);
	for (i=0; i<Ny; i++) fread(Dq[i], sizeof(float), Nx, file);
	for (i=0; i<Ny; i++) fread(h_water[i], sizeof(float), Nx, file);
	for (i=0; i<Ny; i++) fread(h_last_unit[i], sizeof(float), Nx, file);
	for (i=0; i<Ny; i++) fread(EET[i], sizeof(float), Nx, file);
	for (i=0; i<Ny; i++) fread(topo[i], sizeof(float), Nx, file);
	for (i=0; i<Ny; i++) fread(Blocks_base[i], sizeof(float), Nx, file);

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

	if (hydro_model) {
		fread(sortcell, sizeof(struct GRIDNODE), Nx*Ny, file);
		for (i=0; i<Ny; i++) fread(drainage[i], sizeof(struct DRAINAGE), Nx, file);
		for (i=0; i<Ny; i++) fread(lake_former_step[i], sizeof(int), Nx, file);
		if (K_ice_eros) {
			for (i=0; i<Ny; i++) fread(ice_thickness[i], sizeof(float), Nx, file);
			for (i=0; i<Ny; i++) fread(ice_sedm_load[i], sizeof(float), Nx, file);
			for (i=0; i<Ny; i++) fread(ice_velx_sl[i],  sizeof(float), Nx, file);
			for (i=0; i<Ny; i++) fread(ice_vely_sl[i],  sizeof(float), Nx, file);
			for (i=0; i<Ny; i++) fread(ice_velx_df[i],  sizeof(float), Nx, file);
			for (i=0; i<Ny; i++) fread(ice_vely_df[i],  sizeof(float), Nx, file);
		}
	}

	/*numBlocks and i_Block_insert are special*/
	for (j=0; j<numBlocks_aux; j++) {
		insert_new_Block(numBlocks);
		auxptr1 = Blocks[numBlocks-1].thick;
		fread(&Blocks[numBlocks-1], sizeof(struct BLOCK), 1, file);
		Blocks[numBlocks-1].thick = 	auxptr1;
	}
	if (numBlocks_aux != numBlocks) PRINT_ERROR("%d Blocks?!", numBlocks_aux);
	i_Block_insert=i_Block_insert_aux;
	for (j=0; j<numBlocks; j++) {
		for (i=0; i<Ny; i++) 	fread(Blocks[j].thick[i], 	sizeof(float),	Nx, 	file);
	}
	for (j=0; j<numBlocks; j++) {
	    if (Blocks[j].type == 'V') {
		Blocks[j].vel_x = alloc_matrix(Ny, Nx);
		Blocks[j].vel_y = alloc_matrix(Ny, Nx);
		Blocks[j].visc  = alloc_matrix(Ny, Nx);
		Blocks[j].viscTer = alloc_matrix(Ny, Nx);
		for (i=0; i<Ny; i++) {
			fread(Blocks[j].vel_x[i], 	sizeof(float),	Nx, 	file);
			fread(Blocks[j].vel_y[i], 	sizeof(float),	Nx, 	file);
			fread(Blocks[j].visc[i], 	sizeof(float),	Nx, 	file);
			fread(Blocks[j].viscTer[i], 	sizeof(float),	Nx, 	file);
		}
	    }
	    else {
		Blocks[j].vel_x = alloc_matrix(1, 1);
		Blocks[j].vel_y = alloc_matrix(1, 1);
		fread(&Blocks[j].vel_x[0][0], 	sizeof(float),	1, 	file);
		fread(&Blocks[j].vel_y[0][0], 	sizeof(float),	1, 	file);
	    }
	    if (Blocks[j].type == 'S') {
		Blocks[j].detr_ratio = alloc_matrix(Ny, Nx);
		Blocks[j].detr_grsize = alloc_matrix(Ny, Nx);
		for (i=0; i<Ny; i++) {
			fread(Blocks[j].detr_ratio[i], 	sizeof(float),	Nx, 	file);
			fread(Blocks[j].detr_grsize[i], 	sizeof(float),	Nx, 	file);
		}
	    }
	}

	if (erosed_model) {
		for (i=0; i<Ny; i++) fread(eros_now[i], sizeof(float), Nx, file);
		for (i=0; i<Ny; i++) fread(accumul_erosion[i], sizeof(float), Nx, file);
	}
	if (hydro_model) {
		for (i=0; i<Ny; i++) fread(evaporation[i], sizeof(float), Nx, file);
		for (i=0; i<Ny; i++) fread(precipitation[i], sizeof(float), Nx, file);
		for (i=0; i<Ny; i++) fread(precipitation_snow[i], sizeof(float), Nx, file);
		for (i=0; i<Ny; i++) fread(precipitation_file[i], sizeof(float), Nx, file);
		Lake = calloc (nlakes+1, sizeof(struct LAKE_INFO));
		fread(Lake, sizeof(struct LAKE_INFO), nlakes+1, file);
		for (j=1; j<=nlakes; j++) {
			Lake[j].row = calloc (Lake[j].n, sizeof(int));
			Lake[j].col = calloc (Lake[j].n, sizeof(int));
			fread(Lake[j].row, 	sizeof(int), Lake[j].n, file);
			fread(Lake[j].col, 	sizeof(int), Lake[j].n, file);
			Lake[j].row_sd = calloc (Lake[j].n_sd, sizeof(int));
			Lake[j].col_sd = calloc (Lake[j].n_sd, sizeof(int));
			fread(Lake[j].row_sd, 	sizeof(int), Lake[j].n_sd, file);
			fread(Lake[j].col_sd, 	sizeof(int), Lake[j].n_sd, file);
		}
	}

	fread(&end_check, 	sizeof(int),		1, 	file);
	if (end_check != 123456) {
		PRINT_ERROR("\aChecking the end of resume file (%d) failed!\n", end_check);
		exit(0);
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
			PRINT_DEBUG("standard output redirected to %s.\n", filename);
		}
		stdout=file2;
	}

	fclose(file);
	return 1;
}



int read_file_sea_level()
{
	/*
	  Reads file with sea level along time named 'projectname.SLV'
	*/

	int 	i, j, nmax_input_points=5000;
	FILE 	*file;
	float	*aux1, *aux2, *aux3;

	sea_level = 0;

	Read_Open_Filename_Return(".SLV", "rt", "Sea level")

	n_sea_level_input_points=n_eros_level_input_points=0;
	aux1 = calloc(nmax_input_points, sizeof(float));
	aux2 = calloc(nmax_input_points, sizeof(float));
	aux3 = calloc(nmax_input_points, sizeof(float));
	for (i=0; i<nmax_input_points; i++) aux3[i]=NO_DATA;
	
	for (;;) {
		{char auxstr[MAXLENLINE], *lin; int nfields=0; while (nfields<2) {lin=fgets(auxstr, MAXLENLINE-1, file); if (lin==NULL) break; nfields=sscanf(lin, "%f %f %f", &aux1[n_sea_level_input_points], &aux2[n_sea_level_input_points], &aux3[n_sea_level_input_points]);}; if (lin==NULL) break;}
		n_sea_level_input_points++;
		if (aux3[n_sea_level_input_points-1]!=NO_DATA) 
			n_eros_level_input_points++;
		if (n_sea_level_input_points>=nmax_input_points-1 ) {
			PRINT_ERROR("Too many points (%d) in sea level file.", n_sea_level_input_points);
			break;
		}
	}
	fclose(file); 
	var_sea_level = calloc(n_sea_level_input_points, sizeof(float *));
	var_eros_level = calloc(n_eros_level_input_points, sizeof(float *));
	for (i=0, j=0; i<n_sea_level_input_points; i++) {
		var_sea_level[i] = calloc(2, sizeof(float));
		var_sea_level[i][0] = aux1[i]*Matosec;
		var_sea_level[i][1] = aux2[i];
		if (aux3[i]!=NO_DATA) {
			var_eros_level[j] = calloc(2, sizeof(float));
			var_eros_level[j][0] = aux1[i]*Matosec;
			var_eros_level[j][1] = aux3[i];
			j++;
		}
	}
	free(aux1); free(aux2); free(aux3);
	PRINT_INFO("Rows in sea level file '%d'.", n_sea_level_input_points);
	return(1);
}



int read_file_node_defs(float dt_st)
{
	FILE 	*file;
	char 	filename[MAXLENFILE];

	/* 
	  READ THE FILE *.NDEF containing changes to be applied to 
	  certain nodes, e.g. in water flow.
	*/
	
	if (!erosed_model) return(0);

	sprintf(filename, "%s.NDEF", projectname);
	if ((file = fopen(filename, "rt")) == NULL) {
		PRINT_DEBUG("Cannot find file '%s' with node specifications", filename);
		return(0);
	}
	PRINT_DEBUG("Reading Node definitions at '%s'", filename);

	for (;;) {
		int i,j, type;
		float x, y, type_aux, value;
		TAKE_LINE_4(x, y, type_aux, value)
		i = floor((ymax-y)/dy + .49999999);
		j = floor((x-xmin)/dx + .49999999);
		type = (int) type_aux;
		DOMAIN_LIMIT(i,j);
		PRINT_DEBUG("Node specified: %d,%d, %d, %f", i,j, type, value);
		switch (type) {
		    case 1:
			drainage[i][j].discharge += value;
			total_rain += value;
			break;
		    case 2:
			drainage[i][j].masstr += value;
			total_bedrock_eros_mass += value*dt_st;
			break;
		}
	}
	fclose(file);

	return 1;
}




int read_file_Te()
{
	int	i, j, mode_interp_local=mode_interp;
	FILE 	*file;
	BOOL 	switch_EET_file=NO;
	char 	filename[MAXLENFILE];

	/* 
	  READ THE ELASTIC THICKNESS FILE *.EET and writes 
	  interpolated Te in *.eeth
	*/
	
	if (isost_model<=0) return(0);

	for (i=0;i<Ny;i++) for (j=0; j<Nx; j++) {
		EET[i][j]=Te_default;
	}

	sprintf(filename, "%s"".EET", projectname);
	if ((file = fopen(filename,"rt")) == NULL) {
		PRINT_DEBUG("Cannot read Elastic Thickness input file '%s'.", filename);
	}
	else {
		float	z_default=Te_default;
		readinterp2D(file, EET, mode_interp_local, z_default, xmin, xmax, ymin, ymax, Nx, Ny);
		switch_EET_file=YES;
		fclose(file);
	}

	for (i=Te_default=0; i< Ny; i++) {
		for (j=0; j<Nx; j++) {
			EET[i][j] = MAX_2 (EET[i][j], 0);
			D[i][j] = ET2RIG(EET[i][j]) ;
			Te_default += EET[i][j];
		}
	}
	Te_default /= (Nx*Ny);

	/*Writes Te output file*/
	sprintf(filename, "%s.eeth", projectname); 
	remove(filename);
	if (switch_EET_file && switch_ps && switch_write_file_Blocks) {
		Write_Open_Filename_Return (".eeth", "wt", 0);
		for (i=0; i<Ny; i++) {
			for (j=0; j<Nx; j++) {
				fprintf(file, "%7.1f\t%7.1f\t%.1f\n", (xmin+j*dx)/1000, (ymin+(Ny-i-1)*dy)/1000, EET[i][j]);
			}
		}
		fclose(file);
	}

	return 1;
}




int read_file_output_Blocks ()
{
	int 	i, j, k, num_fields;
	FILE 	*file ;
	float 	longperfil, hori, horiant, dens[160], x, y, xa, ya;
	char	auxstr[MAXLENLINE], *lin, word[MAXLENLINE];

	/*
	  This routine is not used by tisc, but by related programs such 
	  as cuthrz or gravanom_3D	  
	  READS A FILE WITH ELEVATION OF HORIZONS IN COLUMNS.

	  WONT WORK WITH COPACTION!!
	*/

	/*Horizons file*/
	Read_Open_Filename_Return(".hrz", "rt", "3D Blocks")

	fgets(auxstr, MAXLENLINE-1, file);
	fgets(auxstr, MAXLENLINE-1, file);
	/*Counts number of Blocks and Nx,Ny, and finds densities in line #2*/
	numBlocks = -3 + sscanf(auxstr, 
		"%s %s %f "
		"%f %f %f %f %f %f %f %f %f %f   %f %f %f %f %f %f %f %f %f %f   %f %f %f %f %f %f %f %f %f %f   %f %f %f %f %f %f %f %f %f %f   %f %f %f %f %f %f %f %f %f %f " 
		"%f %f %f %f %f %f %f %f %f %f   %f %f %f %f %f %f %f %f %f %f   %f %f %f %f %f %f %f %f %f %f   %f %f %f %f %f %f %f %f %f %f   %f %f %f %f %f %f %f %f %f %f " 
		"%f %f %f %f %f %f %f %f %f %f   %f %f %f %f %f %f %f %f %f %f   %f %f %f %f %f %f %f %f %f %f   %f %f %f %f %f %f %f %f %f %f   %f %f %f %f %f %f %f %f %f %f " 
		"%f %f %f %f %f %f %f %f %f %f ", 
		word, word, &denscrust, 
		&dens[0],  &dens[1],  &dens[2],  &dens[3],  &dens[7],  &dens[6],  &dens[6],  &dens[7],  &dens[8],  &dens[9], 
		&dens[10], &dens[11], &dens[12], &dens[13], &dens[14], &dens[15], &dens[16], &dens[17], &dens[18], &dens[19], 
		&dens[20], &dens[21], &dens[22], &dens[23], &dens[24], &dens[25], &dens[26], &dens[27], &dens[28], &dens[29], 
		&dens[30], &dens[31], &dens[32], &dens[33], &dens[34], &dens[35], &dens[36], &dens[37], &dens[38], &dens[39], 
		&dens[40], &dens[41], &dens[42], &dens[43], &dens[44], &dens[45], &dens[46], &dens[47], &dens[48], &dens[49], 
		&dens[50], &dens[51], &dens[52], &dens[53], &dens[54], &dens[55], &dens[56], &dens[57], &dens[58], &dens[59], 
		&dens[60], &dens[61], &dens[62], &dens[63], &dens[64], &dens[65], &dens[66], &dens[67], &dens[68], &dens[69], 
		&dens[70], &dens[71], &dens[72], &dens[73], &dens[74], &dens[75], &dens[76], &dens[77], &dens[78], &dens[79], 
		&dens[80], &dens[81], &dens[82], &dens[83], &dens[84], &dens[85], &dens[86], &dens[87], &dens[88], &dens[89], 
		&dens[90], &dens[91], &dens[92], &dens[93], &dens[94], &dens[95], &dens[96], &dens[97], &dens[98], &dens[99], 
		&dens[100], &dens[101], &dens[102], &dens[103], &dens[104], &dens[105], &dens[106], &dens[107], &dens[108], &dens[109], 
		&dens[110], &dens[111], &dens[112], &dens[113], &dens[114], &dens[115], &dens[116], &dens[117], &dens[118], &dens[119], 
		&dens[120], &dens[121], &dens[122], &dens[123], &dens[124], &dens[125], &dens[126], &dens[127], &dens[128], &dens[129], 
		&dens[130], &dens[131], &dens[132], &dens[133], &dens[134], &dens[135], &dens[136], &dens[137], &dens[138], &dens[139], 
		&dens[140], &dens[141], &dens[142], &dens[143], &dens[144], &dens[145], &dens[146], &dens[147], &dens[148], &dens[149], 
		&dens[150], &dens[151], &dens[152], &dens[153], &dens[154], &dens[155], &dens[156], &dens[157], &dens[158], &dens[159] 
		);
	fgets(auxstr, MAXLENLINE-1, file);
	sscanf(auxstr, "%f %f", &xmin, &ymax);
	xa=xmin; ya=ymax;
	Nx=Ny=1;
	while (1) {
		TAKE_LINE_3(x, y, hori)
		if (y==ymax) Nx++;
		if (x<xa) Ny++;
		xa=x; ya=y;
	}
	xmin*=1e3; ymax*=1e3; 
	xmax=x*1e3; ymin=y*1e3;
	dx = (xmax-xmin)/(Nx-1);
	dy = (ymax-ymin)/(Ny-1);
	Blocks_base = alloc_matrix(Ny, Nx);
	topo = alloc_matrix(Ny, Nx);
	Blocks =  (struct BLOCK *) calloc(NmaxBlocks, sizeof(struct BLOCK));
	for (k=0; k<numBlocks; k++) {
		Blocks[k].thick = alloc_matrix(Ny, Nx);
		Blocks[k].density = dens[k];
	}
	if (verbose_level>=1) {
		fprintf (stderr, 
			"\nnumBlocks=%d;"
			"\nxmin/xmax/ymin/ymax = %.1f/%.1f/%.1f/%.1f m"
			"\nNx=%d, Ny=%d"
			"\ndenscrust = %.1f kg/m3",
			numBlocks, xmin,xmax,ymin,ymax, Nx,Ny, denscrust);
		for (i=0; i<numBlocks; i++) 
			fprintf (stderr, "\ndens[%d] = %.1f kg/m3", i, Blocks[i].density);
	}

	rewind(file);

	fgets(auxstr, MAXLENLINE-1, file);
	fgets(auxstr, MAXLENLINE-1, file);
	for (i=0; i<Ny; i++)  for (j=0; j<Nx; j++) {
		fscanf(file, "%f %f %f", &x, &y, &Blocks_base[i][j]);
		for (k=0,horiant=Blocks_base[i][j]; k<numBlocks; k++) {
			fscanf(file, "%f", &hori);
			Blocks[k].thick[i][j] = hori-horiant;
			if (Blocks[k].thick[i][j]<0) PRINT_WARNING("negative thickness in Block %d at %.1f,%.1f.", k, x, y) ;
			horiant=hori;
		}
		topo[i][j] = hori;
	}
	fclose(file);
	return 1;
}




int read_file_2D_CS (struct BLOCK *Blocks, struct CS2D *CrossSection, int Nx2D)
{
	/* 
	  READ THE CROSS SECTION POLIGON FILE *.PRFL 
	*/
	int	i, j, i2D, num_fields, nmax_input_points=10000, 
		n_2D_CS_pol;			/*Sides of the cross section poligon*/ 
	FILE 	*file;
	char	auxstr[MAXLENFILE], *lin ;
	float	x, y, dl, l=0, pol_length=0, 
		*x_2D_CS_pol,			/*Cross section poligon (x)*/ 
		*y_2D_CS_pol;			/*Cross section poligon (y)*/ 

	{
	    char name[MAXLENFILE];
	    sprintf(name, "%s"".PRFL", projectname);
	    if ((file = fopen(name,"rt")) == NULL) {
	    	PRINT_DEBUG("Cannot read ""Cross-section poligon"" input file '%s'.", name);
	    	return 0;
	    }
	    PRINT_INFO("Reading Cross section poligon"" at '%s'", name);
	}

	x_2D_CS_pol = alloc_array(nmax_input_points);
	y_2D_CS_pol = alloc_array(nmax_input_points);

	n_2D_CS_pol = 0;

	while (1) {
		while (1) {
			lin=fgets(auxstr, MAXLENLINE-1, file);
			if (lin==NULL) break;
			if ((num_fields = sscanf(lin, "%f %f", &x, &y)) >= 2) break;
		}
		if (lin==NULL)	break;
		x_2D_CS_pol[n_2D_CS_pol]=x;
		y_2D_CS_pol[n_2D_CS_pol]=y;
		n_2D_CS_pol++;
	}
	if (n_2D_CS_pol<2) {
		PRINT_WARNING("Need at least 2 x-y points in '%s.PRFL'", projectname); 
		return(0);
	}

	/*Define the points of the cross-section poligon in which to interpolate*/
	for (i=0; i<n_2D_CS_pol-1; i++) 
		pol_length += sqrt(SQUARE(x_2D_CS_pol[i+1]-x_2D_CS_pol[i]) + SQUARE(y_2D_CS_pol[i+1]-y_2D_CS_pol[i]));
	dl = pol_length / (Nx2D-1);
	for (i2D=0; i2D<Nx2D; i2D++) {
		float l_sides, ls=0;
		int is;
		for (is=0, l_sides=0; is<n_2D_CS_pol-1; is++) {
			ls = sqrt(SQUARE(x_2D_CS_pol[is+1]-x_2D_CS_pol[is]) + SQUARE(y_2D_CS_pol[is+1]-y_2D_CS_pol[is]));
			if (l_sides+ls>l) break;
			l_sides += ls;
		}
		CrossSection[i2D].x = x_2D_CS_pol[is] + (l - l_sides) * (x_2D_CS_pol[is+1]-x_2D_CS_pol[is])/ls;
		CrossSection[i2D].y = y_2D_CS_pol[is] + (l - l_sides) * (y_2D_CS_pol[is+1]-y_2D_CS_pol[is])/ls;
		CrossSection[i2D].l = l;
		l += dl;
	}

	PRINT_INFO("Poligon in '%s.PRFL' has %d points and is %.2f km long.", projectname, n_2D_CS_pol, pol_length/1e3);

	fclose(file);
	free(x_2D_CS_pol);
	free(y_2D_CS_pol);

	return 1;
}






/*******************************  OUTPUT  *************************************/


int Calculate_2D_Cross_Section (struct BLOCK *Blocks, struct CS2D *CrossSection, int Nx2D)
{
	float	**hori_aux, **thickness_above;

	hori_aux = alloc_matrix(Ny, Nx);
	thickness_above = alloc_matrix(Ny, Nx);

	/*Blocks_base horizon:*/
	for (int i=0; i<Ny; i++)  for (int j=0; j<Nx; j++) {
		hori_aux[i][j] = Blocks_base[i][j] - w[i][j];
	}
	for (int i2D=0; i2D<Nx2D; i2D++) {
		CrossSection[i2D].horiz[0] =
			interpol_point_in_mesh (hori_aux, Nx, Ny, xmin, dx, ymin, dy, CrossSection[i2D].x, CrossSection[i2D].y) ;
	}

/*
		float thickness_above=0, top_block;
		fprintf(file, "\n%7.2f\t%7.2f", (xmin+j*dx)/1000, (ymax-i*dy)/1000);
		for (int i_Block=0; i_Block<numBlocks; i_Block++) 
			thickness_above += Blocks[i_Block].thick[i][j];
		top_block = Blocks_base[i][j]-w[i][j];
		fprintf(file, "\t%.1f",  top_block);
		for (int i_Block=0; i_Block<numBlocks; i_Block++) {
			comp
			thickness_above -= Blocks[i_Block].thick[i][j];
			top_block += Blocks[i_Block].thick[i][j];
			if (Blocks[i_Block].density==denssedim) top_block -= compaction(sed_porosity, compact_depth, thickness_above, thickness_above+Blocks[i_Block].thick[i][j]);
			fprintf(file, "\t%8.1f",  top_block);
		}
*/


	/*Block horizons:*/
	for (int i_Block=0; i_Block<numBlocks; i_Block++) 
		for (int i=0; i<Ny; i++)  for (int j=0; j<Nx; j++) 
			thickness_above[i][j] += Blocks[i_Block].thick[i][j];
	for (int i_Block=0; i_Block<numBlocks; i_Block++) {
		for (int i=0; i<Ny; i++)  for (int j=0; j<Nx; j++) {
			thickness_above[i][j] -= Blocks[i_Block].thick[i][j];
			hori_aux[i][j] += Blocks[i_Block].thick[i][j];
			if (Blocks[i_Block].density==denssedim) hori_aux[i][j] -= compaction(sed_porosity, compact_depth, thickness_above[i][j], thickness_above[i][j]+Blocks[i_Block].thick[i][j]);
		}
		for (int i2D=0; i2D<Nx2D; i2D++) {
			CrossSection[i2D].horiz[i_Block+1] =
				interpol_point_in_mesh (hori_aux, Nx, Ny, xmin, dx, ymin, dy, CrossSection[i2D].x, CrossSection[i2D].y) ;
		}
	}
    /*========Added By Michael Berry, adding lakes to cross section ===============*/
    if (hydro_model){
            for (int i=0; i<Ny; i++)  for (int j=0; j<Nx; j++){
                    hori_aux[i][j] += h_water[i][j];
            }
            for (int i2D=0; i2D<Nx2D; i2D++) {
                    CrossSection[i2D].horiz[numBlocks+1] =
                            interpol_point_in_mesh (hori_aux, Nx, Ny, xmin, dx, ymin, dy, CrossSection[i2D].x, CrossSection[i2D].y) ;
            }
    }
    /*==========end of addition ===================================================*/
	/*ice horizon:*/
	if (K_ice_eros) {
		for (int i=0; i<Ny; i++)  for (int j=0; j<Nx; j++) {
			hori_aux[i][j] += ice_thickness[i][j];
		}
		for (int i2D=0; i2D<Nx2D; i2D++) {
			CrossSection[i2D].horiz[numBlocks+2] =
				interpol_point_in_mesh (hori_aux, Nx, Ny, xmin, dx, ymin, dy, CrossSection[i2D].x, CrossSection[i2D].y) ;
		}
	}

	free_matrix(hori_aux, Ny);
	free_matrix(thickness_above, Ny);
	return(1);
}




int write_file_cross_section ()
{
	int 	i, j,
		Nx2D=1001; 	/*Number of points of the 2D selected profile*/
	FILE 	*file ;
	BOOL 	switch_CrossSection;
	struct CS2D *CrossSection;
	/*
	  CALCULATES AND WRITES 2D CROSS SECTION FILE
	*/

	CrossSection = (struct CS2D *) calloc(Nx2D, sizeof(struct CS2D));
	for (i=0; i<Nx2D; i++)  CrossSection[i].horiz = (float *) calloc(numBlocks+3, sizeof(float));

	switch_CrossSection = read_file_2D_CS(Blocks, CrossSection, Nx2D);

	if (!switch_CrossSection) {for (i=0; i<Nx2D; i++)  free(CrossSection[i].horiz); free(CrossSection);}
	Write_Open_Filename_Return (".pfl", "wt", !switch_CrossSection);
	Calculate_2D_Cross_Section (Blocks, CrossSection, Nx2D);

	fprintf(file, "# x(km)  \t y(km)  \t  long(km)\t   z(m)-->\n#\t\t\t\t  Densities->\t%8.0f", denscrust) ;
	for (i=0; i<numBlocks; i++) fprintf(file, "\t%8.0f", Blocks[i].density);

        /*====Added by Michael Berry====*/
        if (hydro_model) fprintf(file,"\t%8.0f", denswater);
        /*===end of added===============*/

	if (K_ice_eros) fprintf(file, "\t%8.0f", densice);
	fprintf(file, "\n#\t\t\t\t  Ages->\t%8.2f", Timeini/Matosec);
	for (i=0; i<numBlocks; i++) fprintf(file, "\t%8.2f", Blocks[i].age/Matosec);
	for (i=0; i<Nx2D; i++) {
		fprintf(file, "\n%8.1f\t%8.1f\t%8.1f  ",
		CrossSection[i].x/1000,
		CrossSection[i].y/1000,
		CrossSection[i].l/1000) ;
		for (j=0; j<numBlocks+1; j++) fprintf(file, "\t%8.1f", CrossSection[i].horiz[j] );
            /*===========Added by Michael Berry====*/
            if (hydro_model) fprintf(file,"\t%8.1f", CrossSection[i].horiz[numBlocks+1]);
            /*===========end of addition===========*/
		if (K_ice_eros) fprintf(file, "\t%8.1f", CrossSection[i].horiz[numBlocks+2] );
	}
	fclose(file);

	if (verbose_level>=1) {
		float max2D=-1e9, min2D=1e9, *basam2D;
		basam2D = alloc_array(Nx2D);
		for(i=0; i<Nx2D; i++)  basam2D[i] = CrossSection[i].horiz[numBlocks];
		Perfil_info(basam2D, Nx2D, &max2D, &min2D);
		fprintf(stdout, "\n  2D prof. :  max = %9.1f m     min = %9.1f m   ", max2D, min2D);
		free(basam2D);
	}

	for (i=0; i<Nx2D; i++)  free(CrossSection[i].horiz);
	free (CrossSection);
	return 1;
}




int write_file_Blocks ()
{
	FILE 	*file ;

	/*
	  WRITES A FILE WITH ELEVATION OF ELEVATION OF SURFACES BETWEEN BLOCKS IN COLUMNS
	*/

	Write_Open_Filename_Return (".hrz", "wt", !switch_write_file_Blocks);

	fprintf(file, "# x(km)\t\t y(km)\t\t z(m)-->  \t\t(t=%.2f My)\n#    \tDens:\t%.0f", Time/Matosec, denscrust) ;
	for (int k=0; k<numBlocks; k++) {
		fprintf(file, "\t%.0f", Blocks[k].density);
	}
	if (erosed_model>=2) fprintf(file, "\t   water");
	for (int i=0; i<Ny; i++)  for (int j=0; j<Nx; j++) {
		float thickness_above=0, top_block;
		fprintf(file, "\n%9.3f\t%9.3f", (xmin+j*dx)/1000, (ymax-i*dy)/1000);
		for (int i_Block=0; i_Block<numBlocks; i_Block++) 
			thickness_above += Blocks[i_Block].thick[i][j];
		top_block = Blocks_base[i][j]-w[i][j];
		fprintf(file, "\t%.1f",  top_block);
		for (int i_Block=0; i_Block<numBlocks; i_Block++) {
			thickness_above -= Blocks[i_Block].thick[i][j];
			top_block += Blocks[i_Block].thick[i][j];
			if (Blocks[i_Block].density==denssedim) top_block -= compaction(sed_porosity, compact_depth, thickness_above, thickness_above+Blocks[i_Block].thick[i][j]);
			fprintf(file, "\t%8.1f",  top_block);
		}
/*		if (erosed_model>=2) {
			float top_water=top_block;
			if (drainage[i][j].lake) top_water = Lake[drainage[i][j].lake].alt;
			fprintf(file, "\t%8.1f", top_water);
		}
*/	}
	fclose(file);
	return 1;
}



#define NDERS 8	

int write_file_drainage ()
{
	int 	i, j, k;
	FILE 	*file;
	int 	ro[NDERS], co[NDERS];

	/*
	  WATER and MASS discharge FILE
	*/

#ifdef SURFACE_TRANSPORT
	Write_Open_Filename_Return (".xyw", "wt", !switch_write_file_Blocks || Time==Timeini || !hydro_model);

	fprintf(file, "#TISC output: drainage.  sea_level: %.1f m\n# x(km) y(km) water(m3/s) sed[kg/s] type topo[m] x-to y-to topo-to precipt[mm/y] evapora[mm/y] ice_thick[m] ice_sed_load[m] swim_dist[km]\n", sea_level);
	for (i=0; i<Ny; i++) for (j=0; j<Nx; j++) {
		int il, dcol=drainage[i][j].dr_col, drow=drainage[i][j].dr_row, switch_mouth, ik, jk;
		char dr_type;
		float dist, maxdist;
		float **done;
		done = (float **) calloc(Ny, sizeof(float *));
		for (int i=0; i<Ny; i++) done[i] = (float *) calloc(Nx, sizeof(float));
		maxdist=0;
		if (drow==SIGNAL && dcol==SIGNAL) {drow=i; dcol=j;}
		dr_type = drainage[i][j].type;
		if ((il=drainage[i][j].lake)) 
		  IF_LAKE_IS_SEA(il) dr_type = 'S';
		if (!dr_type) 
		  dr_type = '-';
		/*Find the maximum swimming distance with the neigbouring cells:*/
		ro[0]=i-1, ro[1]=i,   ro[2]=i+1, ro[3]=i,   ro[4]=i-1, ro[5]=i+1, ro[6]=i+1, ro[7]=i-1;
		co[0]=j,   co[1]=j+1, co[2]=j,   co[3]=j-1, co[4]=j+1, co[5]=j+1, co[6]=j-1, co[7]=j-1;
		/*Follow down the i,j river and mark the path to a mouth*/
		ik=i; jk=j; 
		dist=0;
		for (int k=1; ; k++) {
			int id, jd, ild;
//printf("\n&&&&&%d    %d %d  %d %d\n", k, i, j, ik, jk);
			id=drainage[ik][jk].dr_row; jd=drainage[ik][jk].dr_col;
			il=drainage[ik][jk].lake;
			if (IN_DOMAIN(id,jd))	    ild=drainage[id][jd].lake;
			else			    ild=0;
			if (IN_DOMAIN(id,jd)) dist += sqrt((id-ik)*(id-ik)*dy*dy+(jd-jk)*(jd-jk)*dx*dx);
			done[ik][jk]=dist; /*mark this node with the distance to i,j*/
			switch_mouth=NO; 
			/*If this is not a lake and it drains to the sea, then it is a river mouth*/
			if (!il) IF_LAKE_IS_SEA(ild)
			    switch_mouth=YES;
			/*If it drains outside, it is the end of a river*/
			if (OUT_DOMAIN(id,jd))
				    switch_mouth=YES;
			/*If it drains to an endorheic lake, it is the end of a river*/
			if (!il && ild && !Lake[ild].n_sd)
				    switch_mouth=YES;
			if (switch_mouth) break;
/*WHY?!!*/		if (k>Nx*Ny) {break;}
			ik=id; jk=jd;
		}
		/*Follow down from the neigbours and take the max swimming distance*/
		for (int l=0; l<NDERS; l++) {
			float distneighb;
			ik=ro[l]; jk=co[l]; 
			distneighb=0;
			if (IN_DOMAIN(ik,jk)) for (int k=1; ; k++) {
				int id, jd, ild;
				id=drainage[ik][jk].dr_row; jd=drainage[ik][jk].dr_col;
				il=drainage[ik][jk].lake;
				if (IN_DOMAIN(id,jd))	    ild=drainage[id][jd].lake;
				else			    ild=0;
				if (IN_DOMAIN(id,jd)) distneighb += sqrt((id-ik)*(id-ik)*dy*dy+(jd-jk)*(jd-jk)*dx*dx);
				switch_mouth=NO; 
				/*If this is not a lake and it drains to the sea, then it is a river mouth*/
				if (!il) IF_LAKE_IS_SEA(ild)
				    switch_mouth=YES;
				/*If it drains outside, it is the end of a river*/
				if (OUT_DOMAIN(id,jd))
					    switch_mouth=YES;
				/*If it drains to an endorheic lake, it is the end of a river*/
					if (!il && ild && !Lake[ild].n_sd)
					    switch_mouth=YES;
				if (switch_mouth) {distneighb+=dist; break;}
				if (done[id][jd]) {distneighb+=done[id][jd]; break;}
/*WHY?!!*/		if (k>Nx*Ny) {break;}
				ik=id; jk=jd;
			}
			if (distneighb>maxdist) maxdist = distneighb;
		}
		
		fprintf(file, "%7.2f\t%7.2f\t%.2f\t%.2f\t%c\t%.1f\t%7.2f\t%7.2f\t%.1f\t%7.1f\t%7.1f\t%6.1f\t%6.2f\t%6.2f\n",
			(xmin+j*dx)/1000, (ymax-i*dy)/1000,
			drainage[i][j].discharge,
			drainage[i][j].masstr,
			dr_type, 
			topo[i][j], 
			(xmin+dcol*dx)/1000, (ymax-drow*dy)/1000, topo[drow][dcol], 
			precipitation[i][j]*secsperyr*1e3, evaporation[i][j]*secsperyr*1e3, 
			(K_ice_eros)? ice_thickness[i][j] : 0, 
			(K_ice_eros)? ice_sedm_load[i][j] : 0,
			maxdist/1e3 );
		for (int i=0; i<Ny; i++) free(done[i]);
		free(done);
	}
	fclose(file);
#endif
	return (1);
}




int find_up_river (int row, int col, int *level, int *count, float *length, float *chi, FILE *file, BOOL **done, float ref_discharge)
{
	int 	n_above=1;
	float	tp, dtp;

	/*Just in case there is a feedback bucle by mistake, track the checked nodes*/
	if (done[row][col]) {
		fprintf(file, "\tERROR: drainage loop in [%d][%d]", row, col); 
		PRINT_WARNING("drainage loop in [%d][%d]", row, col); 
		return(0);
	}
	done[row][col] = YES;

	(*level) ++;
	for (int i=0; i<Ny; i++)  for (int j=0; j<Nx; j++) {
		if (drainage[i][j].dr_row == row && drainage[i][j].dr_col == col) {
			float 	Dl = sqrt((col-j)*dx*(col-j)*dx + (row-i)*dy*(row-i)*dy), 
				weight = pow(drainage[row][col].discharge/ref_discharge, -spl_m/spl_n);
			(*length) += Dl;
			(*chi)    += Dl*weight;
			n_above   += find_up_river (i, j, level, count, length, chi, file, done, ref_discharge);
			(*length) -= Dl;
			(*chi)    -= Dl*weight;
		}
	}
	(*count) ++;
	{
	int drow, dcol; float dl, weight;
	drow=drainage[row][col].dr_row;
	dcol=drainage[row][col].dr_col;
	if (OUT_DOMAIN(drow,dcol)) {drow=row; dcol=col;}
	dl = sqrt((dcol-col)*dx*(dcol-col)*dx + (drow-row)*dy*(drow-row)*dy);
	tp  = topo[ row][ col];	if (drainage[ row][ col].lake)  tp = Lake[drainage[ row][ col].lake].alt;
	dtp = topo[drow][dcol];	if (drainage[drow][dcol].lake) dtp = Lake[drainage[drow][dcol].lake].alt;
	weight=pow(drainage[drow][dcol].discharge/ref_discharge, -spl_m/spl_n);
	fprintf(file, "\n%7.2f\t%7.2f\t%.2f\t%.2f\t%c\t%.1f",
		(xmin+col*dx)/1000,
		(ymax-row*dy)/1000,
		drainage[row][col].discharge,
		drainage[row][col].masstr,
		drainage[row][col].type,
		tp);
	fprintf(file, "\t%.2f\t%.2f\t%7.2f\t%7.2f\t%.1f",
		(*length)/1000,
		(*chi)/1000,
		(xmin+dcol*dx)/1000,
		(ymax-drow*dy)/1000,
		dtp); 
	fprintf(file, "\t%.1f\t%.1f",
		(*length - dl)/1000, 
		(*chi - dl*weight)/1000); 
	fprintf(file, "\t%.1f\t%.1f",
		(erosed_model)? eros_now[row][col]/(dt/Matosec)/dx/dy/denscrust  :  0, 
		(erosed_model)? accumul_erosion[row][col]/dx/dy/denscrust  :  0); 
	fprintf(file, "\t%d\t%d",
		*level,
		n_above);
	}
	fflush(file);
	(*level) --;
	return (n_above);
}




int write_file_river_basins ()
{
	int 	id, jd, il, ild, 
		count=0, level, n_river=0, river_nodes;
	float 	length, chi, maxdisch=-1e9;
	FILE 	*file;
	BOOL 	switch_mouth, **done;

	/*
	  WRITES INFORMATION ABOUT HYDROLOGICAL BASINS
	*/

	Write_Open_Filename_Return (".bas", "wt", !hydro_model || !switch_write_file || Time==Timeini);

	done = (BOOL **) calloc(Ny, sizeof(BOOL *));
	for (int i=0; i<Ny; i++) done[i] = (BOOL *) calloc(Nx, sizeof(BOOL));

	for (int i=0; i<Ny; i++) for (int j=0; j<Nx; j++) maxdisch=MAX_2(maxdisch, drainage[i][j].discharge);

	fprintf(file, "#TISC output: river basins.\n#x[km] y[km] dischg[m3/s] masstr[kg/s] type topo[m] length[km] chi[km] x-to y-to topo-to[m] length-to[km] chi-to[km] eros_rate[m/My] accumul_erosion[m] level nodes_above");
	if (verbose_level>=3) fprintf(stdout, "  %3.0f%%", count*100.0/Nx/Ny);

	/*Look for river mouths in sea, closed lakes or boundaries*/
	for (int i=0; i<Ny; i++) for (int j=0; j<Nx; j++) {
	    switch_mouth=NO;
	    id=drainage[i][j].dr_row; jd=drainage[i][j].dr_col;
	    il=drainage[i][j].lake;      
	    if (IN_DOMAIN(id,jd)) 	ild=drainage[id][jd].lake;
	    else 			ild=0;
 	    /*If this is not a lake and it drains to the sea, then it is a river mouth*/
 	    if (!il) IF_LAKE_IS_SEA(ild) 
	    		switch_mouth=YES;
   	    /*If it drains outside (but is not an endorheic lake) and is above sea_level, it is the end of a river*/
	    if (!il && OUT_DOMAIN(id,jd) && topo[i][j] >= sea_level)
 	    		switch_mouth=YES;
 	    /*If it drains to an endorheic lake, it is the end of a river*/
	    if (!il && ild && !Lake[ild].n_sd)
  	    		switch_mouth=YES;
	    if (switch_mouth) {
	    	if (verbose_level>=3) {fprintf(stdout, "\b\b\b\b%3.0f%%", count*100.0/Nx/Ny); fflush(stdout);}
	    	n_river++; level=0; length=chi=0;
	    	fprintf(file, "\n> begin river basin %d", n_river); fflush(file);
	    	river_nodes = find_up_river (i, j, &level, &count, &length, &chi, file, done, maxdisch);
	    	fprintf(file, "\n# END river basin %d of %.0f km2: %.2f m3/s %.2f kg/s at [%.0f,%.0f] km", 
			n_river, river_nodes*dx*dy/1e6, drainage[i][j].discharge, drainage[i][j].masstr, 
			(xmin+j*dx)/1000,  (ymax-i*dy)/1000); 
		fflush(file);
 	    }
	}
	if (verbose_level>=3) {fprintf(stdout, "\b\b\b\b     "); fflush(stdout);}
	fprintf(file, "\n");
	fclose(file);
	for (int i=0; i<Ny; i++) free(done[i]);
	free(done);
	return 1;
}





int write_file_lakes ()
{
	/*
	  WRITES INFORMATION ABOUT LAKES
	*/

	int 	i, j, k;
	FILE 	*file;

#ifdef SURFACE_TRANSPORT
	Write_Open_Filename_Return (".lak", "wt", !hydro_model || !switch_write_file || !nlakes || Time==Timeini);

	fprintf(file, "#TISC output: \n#Lakes: %d", nlakes);
	for (i=1; i<=nlakes; i++) {
		float vol, height_lake, lake_evaporation;
		for (j=0, height_lake=sea_level; j<Lake[i].n; j++)
		    height_lake = MAX_2(height_lake, topo[Lake[i].row[j]][Lake[i].col[j]]);
		for (j=0, vol=0; j<Lake[i].n; j++)
		    vol += height_lake - topo[Lake[i].row[j]][Lake[i].col[j]];
		vol *= dx*dy;
		lake_evaporation = Lake[i].n * dx*dy * evaporation_ct;
		fprintf(file, "\n>Lake %d:  %d nodes", i, Lake[i].n);
		if (Lake[i].n)    fprintf(file, "  %.3f km3  %.1f m  %2d outl.  %.2f m3/s inp", vol/1e9, height_lake, Lake[i].n_sd, Lake_Input_Discharge(i));
		if (Lake[i].n_sd) fprintf(file, "   %.2f m3/s evap.", lake_evaporation);
		for (j=0; j<Lake[i].n; j++) {
			int 	row=Lake[i].row[j], col=Lake[i].col[j], 
				dcol, drow;
			char 	dr_type=drainage[row][col].type;
			dcol=drainage[row][col].dr_col; 
			drow=drainage[row][col].dr_row;
			if (drow==SIGNAL && dcol==SIGNAL) {drow=row; dcol=col;}
			IF_LAKE_IS_SEA(i) dr_type = 'S';
			fprintf(file, "\n%7.2f\t%7.2f\t%c",
				(xmin+col*dx)/1000, (ymax-row*dy)/1000,
				dr_type);
			if (drainage[row][col].type == 'E')
				fprintf(file, "\t%.2f m3/s\t%.2f\t%.2f",
					drainage[row][col].discharge,
					(xmin+dcol*dx)/1000, (ymax-drow*dy)/1000  );
		}
	}
	fclose(file);
#endif
	return 1;
}




int write_file_ice ()
{
	int 	i, j;
	FILE 	*file;

	/*
	  Ice flow file
	*/

#ifdef SURFACE_TRANSPORT
	Write_Open_Filename_Return (".ice", "wt", !switch_write_file_Blocks || !hydro_model || !K_ice_eros || Time==Timeini);

	fprintf(file, "#TISC output: ice flow.  sea_level: %.1f m\n# x(km) y(km) topo[m] ice_thick[m]  vx_df vy_df[m/y]  vx_sl vy_sl  sol_prec[mm/y] ice_sed_load[m]\n", sea_level);
	for (i=0; i<Ny; i++) for (j=0; j<Nx; j++) {
		fprintf(file, "%7.2f\t%7.2f\t%.1f\t%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.1f\t%.1f\n",
			(xmin+j*dx)/1000, (ymax-i*dy)/1000,
			topo[i][j], 
			ice_thickness[i][j], 
			ice_velx_df[i][j]*secsperyr, ice_vely_df[i][j]*secsperyr, 
			ice_velx_sl[i][j]*secsperyr, ice_vely_sl[i][j]*secsperyr, 
			precipitation[i][j]*secsperyr*1e3, 
			ice_sedm_load[i][j]);
	}
	fclose(file);
#endif
	return (1);
}




int write_file_resume()
{
	int 	i, j, i_hori, end_check=123456;
	FILE 	*file ;

	/*
	  WRITES A BINARY FILE WITH ALL THE INFORMATION 
	  REQUIRED TO RESUME A PROJECT RUN
	*/

	Write_Open_Filename_Return (".all", "wb", !switch_write_file || !switch_write_file_Blocks);

	/*Defined in universal.h:*/
	fwrite(&Nx, 		sizeof(int),		1, 	file);
	fwrite(&Ny, 		sizeof(int),		1, 	file);
	fwrite(&Nz, 		sizeof(int),		1, 	file);
	fwrite(&verbose_level,	sizeof(int),		1, 	file);

	fwrite(&dx, 		sizeof(float),		1, 	file);
	fwrite(&dy, 		sizeof(float),		1, 	file);
	fwrite(&dz, 		sizeof(float),		1, 	file);

	fwrite(version, 	sizeof(char),		LENGTHVERS, 	file);
	fwrite(version_input,	sizeof(char),		LENGTHVERS, 	file);

	fwrite(&switch_ps, 	sizeof(BOOL),		1, 	file);
	fwrite(&switch_write_file, sizeof(BOOL),		1, 	file);


	/*Defined in geomodel.h:*/
	fwrite(&grav_anom_type, 	sizeof(int),		1, 	file);
	fwrite(&isost_model, 	sizeof(int),		1, 	file);
	fwrite(&water_load, 	sizeof(int),		1, 	file);

	fwrite(&Te_default, 	sizeof(float),		1, 	file);
	fwrite(&crust_thick_default, sizeof(float),	1, 	file);
	fwrite(&upper_crust_thick_default, sizeof(float),1, 	file);
	fwrite(&densasthen, 	sizeof(float),		1, 	file);
	fwrite(&densmantle, 	sizeof(float),		1, 	file);
	fwrite(&denssedim, 	sizeof(float),		1, 	file);
	fwrite(&denscrust, 	sizeof(float),		1, 	file);
	fwrite(&densenv, 	sizeof(float),		1, 	file);
	fwrite(&densinfill, 	sizeof(float),		1, 	file);
	fwrite(&sea_level, 	sizeof(float),		1, 	file);
	fwrite(&temp_sea_level, 	sizeof(float),		1, 	file);
	fwrite(&Time, 		sizeof(float),		1, 	file);
	fwrite(&Timefinal, 	sizeof(float),		1, 	file);
	fwrite(&Timeini, 	sizeof(float),		1, 	file);
	fwrite(&dt, 		sizeof(float),		1, 	file);
	fwrite(&dt_eros, 	sizeof(float),		1, 	file);
	fwrite(&tau, 		sizeof(float),		1, 	file);

	fwrite(projectname, 	sizeof(char),	MAXLENFILE, 	file);

	fwrite(&switch_geograph_coor, sizeof(BOOL),	1, 	file);

	/*Defined in tao+tisc.h:*/
	fwrite(&nloads, 		sizeof(int),		1, 	file);
	fwrite(&n_sea_level_input_points, sizeof(int),	1, 	file);
	fwrite(&n_eros_level_input_points, sizeof(int),	1, 	file);
	fwrite(&n_record_times, 	sizeof(int),		1, 	file);
	fwrite(&i_first_Block_load, sizeof(int),	1, 	file);
	fwrite(&i_Block_insert, sizeof(int),		1, 	file);
	fwrite(&numBlocks, 	sizeof(int),		1, 	file);
	fwrite(&nwrotenfiles, 	sizeof(int),		1, 	file);
	fwrite(&run_type, 	sizeof(int),		1, 	file);

	fwrite(&zini, 		sizeof(float),		1, 	file);
	fwrite(&dt_record, 	sizeof(float),		1, 	file);
	fwrite(&sed_porosity, 	sizeof(float),		1, 	file);
	fwrite(&compact_depth, 	sizeof(float),		1, 	file);
	fwrite(&Kerosdif, 	sizeof(float),		1, 	file);
	fwrite(&last_time_file_time, 	sizeof(float),		1, 	file);
	fwrite(&random_topo, 	sizeof(float),		1, 	file);

	fwrite(&switch_file_out, 	sizeof(BOOL),		1, 	file);
	fwrite(&switch_gradual, 	sizeof(BOOL),		1, 	file);
	fwrite(&switch_topoest, 		sizeof(BOOL),		1, 	file);
	fwrite(&switch_write_file_Blocks, sizeof(BOOL),		1, 	file);
	fwrite(&deform_sed, sizeof(BOOL),		1, 	file);


	/*Defined in tisc.h:*/
	fwrite(&erosed_model, 	sizeof(int),		1, 	file);
	fwrite(&mode_interp, 	sizeof(int),		1, 	file);
	fwrite(&nbasins, 	sizeof(int),		1, 	file);
	fwrite(&nlakes, 	sizeof(int),		1, 	file);
	fwrite(&n_image, 	sizeof(int),		1, 	file);

	fwrite(&xmin, 		sizeof(float),		1, 	file);
	fwrite(&xmax, 		sizeof(float),		1, 	file);
	fwrite(&ymin, 		sizeof(float),		1, 	file);
	fwrite(&ymax, 		sizeof(float),		1, 	file);
	fwrite(&dxy, 		sizeof(float),		1, 	file);
	fwrite(&Px, 	sizeof(float),		1, 	file);
	fwrite(&Py, 	sizeof(float),		1, 	file);
	fwrite(&Pxy, 	sizeof(float),		1, 	file);
	fwrite(&Keroseol, 	sizeof(float),		1, 	file);
	fwrite(&Ksedim, 		sizeof(float),		1, 	file);
	fwrite(&critical_slope, 	sizeof(float),		1, 	file);
	fwrite(&K_river_cap, 	sizeof(float),		1, 	file);
	fwrite(&K_ice_eros, 	sizeof(float),		1, 	file);
	fwrite(&dt_ice, 	sizeof(float),		1, 	file);
	fwrite(&n_ice_flow, 	sizeof(int),		1, 	file);
	fwrite(&A_ice_rheo, 	sizeof(float),		1, 	file);
	fwrite(&A_ice_slide, 	sizeof(float),		1, 	file);
	fwrite(&erodibility, 	sizeof(float),		1, 	file);
	fwrite(&erodibility_sed, 	sizeof(float),		1, 	file);
	fwrite(&critical_stress, 	sizeof(float),		1, 	file);
	fwrite(&l_fluv_sedim, 	sizeof(float),		1, 	file);
	fwrite(&lost_rate, 	sizeof(float),		1, 	file);
	fwrite(&permeability, 	sizeof(float),		1, 	file);
	fwrite(&evaporation_ct, 	sizeof(float),		1, 	file);
	fwrite(&rain, 	sizeof(float),		1, 	file);
	fwrite(&Krain, 	sizeof(float),		1, 	file);
	fwrite(&relative_humidity, 	sizeof(float),		1, 	file);
	fwrite(&windazimut, 	sizeof(float),		1, 	file);
	fwrite(&CXrain, 	sizeof(float),		1, 	file);
	fwrite(&CYrain, 	sizeof(float),		1, 	file);
	fwrite(&total_bedrock_eros_mass, 	sizeof(float),		1, 	file);
	fwrite(&total_sed_mass, 	sizeof(float),		1, 	file);

	fwrite(&hydro_model, 	sizeof(int),		1, 	file);
	fwrite(&lake_instant_fill, 	sizeof(int),		1, 	file);

	fwrite(boundary_conds, 	sizeof(char),	5, 	file);
	fwrite(eros_bound_cond,	sizeof(char),	5, 	file);
	fwrite(&solver_type, 	sizeof(char),	1, 	file);
	fwrite(gif_geom, 	sizeof(char),	MAXLENLINE, 	file);


	/*Arrays:*/
	for (i=0; i<Ny; i++) fwrite(w[i], sizeof(float), Nx, file);
	for (i=0; i<Ny; i++) fwrite(D[i], sizeof(float), Nx, file);
	for (i=0; i<Ny; i++) fwrite(q[i], sizeof(float), Nx, file);
	for (i=0; i<Ny; i++) fwrite(Dw[i], sizeof(float), Nx, file);
	for (i=0; i<Ny; i++) fwrite(Dq[i], sizeof(float), Nx, file);
	for (i=0; i<Ny; i++) fwrite(h_water[i], sizeof(float), Nx, file);
	for (i=0; i<Ny; i++) fwrite(h_last_unit[i], sizeof(float), Nx, file);
	for (i=0; i<Ny; i++) fwrite(EET[i], sizeof(float), Nx, file);
	for (i=0; i<Ny; i++) fwrite(topo[i], sizeof(float), Nx, file);
	for (i=0; i<Ny; i++) fwrite(Blocks_base[i], sizeof(float), Nx, file);

	fwrite(horiz_record_time, sizeof(float), n_record_times, file);

	for (i=0; i<n_sea_level_input_points; i++) 
		fwrite(var_sea_level[i], sizeof(float), 2, file);
	for (i=0; i<n_eros_level_input_points; i++) 
		fwrite(var_eros_level[i], sizeof(float), 2, file);

	if (hydro_model) {
		fwrite(sortcell, sizeof(struct GRIDNODE), Nx*Ny, file);
		for (i=0; i<Ny; i++) fwrite(drainage[i], sizeof(struct DRAINAGE), Nx, file);
		for (i=0; i<Ny; i++) fwrite(lake_former_step[i], sizeof(int), Nx, file);
		if (K_ice_eros) {
			for (i=0; i<Ny; i++) fwrite(ice_thickness[i], sizeof(float), Nx, file);
			for (i=0; i<Ny; i++) fwrite(ice_sedm_load[i], sizeof(float), Nx, file);
			for (i=0; i<Ny; i++) fwrite(ice_velx_sl[i],  sizeof(float), Nx, file);
			for (i=0; i<Ny; i++) fwrite(ice_vely_sl[i],  sizeof(float), Nx, file);
			for (i=0; i<Ny; i++) fwrite(ice_velx_df[i],  sizeof(float), Nx, file);
			for (i=0; i<Ny; i++) fwrite(ice_vely_df[i],  sizeof(float), Nx, file);
		}
	}

	fwrite(Blocks, 	sizeof(struct BLOCK),	numBlocks, file);
	for (j=0; j<numBlocks; j++) {
		for (i=0; i<Ny; i++) fwrite(Blocks[j].thick[i], sizeof(float), Nx, file);
	}
	for (j=0; j<numBlocks; j++) {
		if (Blocks[j].type == 'V') {
			for (i=0; i<Ny; i++) {
				fwrite(Blocks[j].vel_x[i], sizeof(float), Nx, file);
				fwrite(Blocks[j].vel_y[i], sizeof(float), Nx, file);
				fwrite(Blocks[j].visc[i],  sizeof(float), Nx, file);
				fwrite(Blocks[j].viscTer[i],  sizeof(float), Nx, file);
			}
		}
		else {
				fwrite(&Blocks[j].vel_x[0][0], sizeof(float), 1, file);
				fwrite(&Blocks[j].vel_y[0][0], sizeof(float), 1, file);
		}
		if (Blocks[j].type == 'S') {
			for (i=0; i<Ny; i++) {
				fwrite(Blocks[j].detr_ratio[i], sizeof(float), Nx, file);
				fwrite(Blocks[j].detr_grsize[i], sizeof(float), Nx, file);
			}
		}
	}

	if (erosed_model) {
		for (i=0; i<Ny; i++) fwrite(eros_now[i], sizeof(float), Nx, file);
		for (i=0; i<Ny; i++) fwrite(accumul_erosion[i], sizeof(float), Nx, file);
	}
	if (hydro_model) {
		for (i=0; i<Ny; i++) fwrite(evaporation[i], sizeof(float), Nx, file);
		for (i=0; i<Ny; i++) fwrite(precipitation[i], sizeof(float), Nx, file);
		for (i=0; i<Ny; i++) fwrite(precipitation_snow[i], sizeof(float), Nx, file);
		for (i=0; i<Ny; i++) fwrite(precipitation_file[i], sizeof(float), Nx, file);
		fwrite(Lake, sizeof(struct LAKE_INFO), nlakes+1, file);
		for (j=1; j<=nlakes; j++) {
			fwrite(Lake[j].row, 	sizeof(int), Lake[j].n, file);
			fwrite(Lake[j].col, 	sizeof(int), Lake[j].n, file);
			fwrite(Lake[j].row_sd, 	sizeof(int), Lake[j].n_sd, file);
			fwrite(Lake[j].col_sd, 	sizeof(int), Lake[j].n_sd, file);
		}
	}

	fwrite(&end_check, 	sizeof(int),		1, 	file);

	fclose(file);
	return 1;
}



int write_file_surftransp ()
{
	FILE 	*file;

	/*
	  EROSION and SEDIMENTATION file
	*/

	Write_Open_Filename_Return (".st", "wt", !erosed_model || !switch_write_file_Blocks || Time==Timeini);

	fprintf(file, "#TISC output drainage.  sea_level: %.1f m\n# x(km)  y(km) topo[m]  accumul_erosion[m] eros_rate[m/My]\n", sea_level);
	for (int i=0; i<Ny; i++) for (int j=0; j<Nx; j++) {
		fprintf(file, "%7.2f\t%7.2f\t%.1f\t%.1f\t%.1f\n",
			(xmin+j*dx)/1000, (ymax-i*dy)/1000, topo[i][j], 
			accumul_erosion[i][j] / dx/dy/denscrust, eros_now[i][j]/(dt/Matosec) / dx/dy/denscrust);
	}
	fclose(file);
	return (1);
}




int write_file_time (float **xarxa1, float **xarxa2)
{
	/*
	  WRITES deflection and elevation along time file
	*/

	int 	i, j ;
	FILE 	*file;
	char 	filename[MAXLENLINE], filename1[MAXLENLINE], filename2[MAXLENLINE], 
		command[MAXLENLINE];
	BOOL	return_cond;
	float	youngest_age=-1e16;

	for (i=0; i<numBlocks; i++) youngest_age = MAX_2(Blocks[i].age, youngest_age);
	return_cond = !switch_write_file_Blocks 
		|| !isost_model 
		|| (((Time-last_time_file_time) < dt_record || (!dt_record && youngest_age!=Time)) && (Timefinal-Time) >= dt) 
		||  (Time-last_time_file_time) == 0;

	if (nwrotenfiles==0) {
		Write_Open_Filename_Return (".xyzt", "wt", return_cond);
		fprintf(file, "#\t\tTime: %.2f My\n#x(km)\ty(km)\tw(m)\th_water(m)\n", Time/Matosec);
	}
	else {
		if (return_cond) return (0);
		sprintf(filename, "%s.xyzt", projectname);
		PRINT_INFO("Writing file '%s' (%d times).", filename, nwrotenfiles+1);
		sprintf(filename1, "%s.aux1.xyzt.tisc.tmp", projectname);
		sprintf(filename2, "%s.aux2.xyzt.tisc.tmp", projectname);
		if ((file = fopen(filename1, "wt")) == NULL) {
			PRINT_ERROR("Cannot open auxiliar output file %s.\n", filename1);
			return (0);
		}
		fprintf(file,      "Time: %.2f My\nw(m)\th(m)\n", Time/Matosec);
	}

	for (i=0; i<Ny; i++) {
		for (j=0; j<Nx; j++) {
			if (nwrotenfiles==0) {
				fprintf (file, "%7.2f\t%7.2f", 
					(xmin+j*(xmax-xmin)/(Nx-1))/1000, 
					(ymax-i*(ymax-ymin)/(Ny-1))/1000
				);
				fprintf(file, "\t%.1f\t%.1f\n",  xarxa1[i][j], xarxa2[i][j]);
			}
			if (nwrotenfiles>0) 
				fprintf(file, "%.1f\t%.1f\n",    xarxa1[i][j], xarxa2[i][j]);
		}
	}
	fclose(file);
	if (nwrotenfiles>0) {
		/*Paste the new columns into the .xyzt file.*/
		sprintf(command, 
			"paste %s %s > %s", filename, filename1, filename2);
		system(command);
		rename(filename2, filename); remove(filename1);
	}

	nwrotenfiles++;
	last_time_file_time = Time;

	return (1);
}



int write_file_deflection ()
{
	/*
	  WRITES deflection and elevation
	*/
	FILE 	*file;

	Write_Open_Filename_Return (".xyzt", "wt", !isost_model || !switch_write_file_Blocks);

	fprintf(file, "#TISC output deflection. Te_default: %.1f m\n# x(km)  y(km) w[m] dw/dt[m/Myr] topo[m]\n", Te_default);
	for (int i=0; i<Ny; i++) for (int j=0; j<Nx; j++) {
		fprintf(file, "%7.2f\t%7.2f\t%.1f\t%.1f\t%.1f\n",
			(xmin+j*dx)/1000, (ymax-i*dy)/1000, w[i][j], Dw[i][j]/dt*Matosec, topo[i][j]);
	}
	fclose(file);
	return (1);
}



int write_file_velocity_field ()
{
	int 	i, j, k, thin_sheet_Blocks=0;
	FILE 	*file;

	/*
	  WRITES A FILE WITH THE VELOCITY FIELD AND THE VISCOSITY OF THE 
	  THIN SHEET BlockS
	*/

	for (k=0; k<numBlocks; k++) if (Blocks[k].type == 'V') thin_sheet_Blocks++;
	Write_Open_Filename_Return (".vel", "wt", !switch_write_file_Blocks || !thin_sheet_Blocks);

	fprintf(file, "# x(km)\ty(km) \tvel_x(km/My) vel_y(km/My)  viscosity(Pa s-1)-->  \t\t(t=%.2f My)\n#    \tDens:\t", Time/Matosec) ;
	for (k=0; k<numBlocks; k++) {
		if (Blocks[k].type == 'V') fprintf(file, "\t\t\t%.0f", Blocks[k].density);
	}
	for (i=0; i<Ny; i++)  for (j=0; j<Nx; j++) {
		fprintf(file, "\n%7.2f\t%7.2f", (xmin+j*dx)/1000, (ymax-i*dy)/1000);
		for (k=0; k<numBlocks; k++) {
			if (Blocks[k].type == 'V') fprintf(file, "\t%.2f\t%.2f\t\t%.5e", Blocks[k].vel_x[i][j]/1000*Matosec, Blocks[k].vel_y[i][j]/1000*Matosec, Blocks[k].visc[i][j]);
		}
	}
	fclose(file);
	return 1;
}


