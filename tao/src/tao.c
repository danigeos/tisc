/*******************************************************************************
***********						tAo main source file					********
********************************************************************************
	For compilation and installation check the file ../tao/README.
*******************************************************************************/

#include "tao.h"


int main(int argc, char **argv)
{
	ModelConfig cfg = {0};
	ModelContext ctx = {0};

	/*get input parameters and files*/
	inputs (&cfg, &ctx, argc, argv) ;

	fprintf(stdout, "\nT= %.4f My", Time/Matosec);

	/* Populate configuration context */
	cfg.Nx = Nx; cfg.Nz = Nz;
	cfg.dx = dx; cfg.dz = dz;
	cfg.x0 = x0; cfg.xf = xf; cfg.xmin = xmin; cfg.xmax = xmax;
	cfg.hydro_model = hydro_model;
	cfg.erosed_model = erosed_model;
	cfg.verbose_level = verbose_level;
	cfg.densenv = densenv; cfg.denssedim = denssedim;
	cfg.denscrust = denscrust; cfg.denswater = denswater; 
	cfg.densasthen = densasthen; cfg.densmantle = densmantle;
	cfg.riverbasinwidth = riverbasinwidth;
	cfg.sed_porosity = sed_porosity;

	ctx.Timeini = Timeini; ctx.Timefinal = Timefinal;
	ctx.dt_eros = dt_eros;
	ctx.topo = topo;
	ctx.numBlocks = numBlocks;
	ctx.Time = Time;
	ctx.dt = dt;
	ctx.sea_level = sea_level;
	ctx.nlakes = nlakes;

	if (switch_ps>=2) {calculate_topo(&cfg, &ctx, topo); Write_Ouput(&cfg, &ctx);}

	/*MAIN LOOP: In this loop time increases from Timeini to Timefinal*/
	do { 
		/*Remember topography before tectonics and flexure*/
		calculate_topo(&cfg, &ctx, topo);

		/*Calculate tectonic deformation and tectonic load*/
		tectload(&cfg, &ctx);

		/*Sea level variations*/
		ctx.sea_level = calculate_sea_level(ctx.Time);

		/*Define & solve elastic equation*/
		Elastoplastic_Deflection(&cfg, &ctx);

		/*Define & solve viscoelastic flexure equation*/
		Viscous_Relaxation(&cfg, &ctx);

		/* Update dynamic context */
		ctx.Time = Time;
		ctx.dt = dt;
		ctx.sea_level = sea_level;

		/*Calculates surface water flow and sediment-erosion*/
		surface_processes(&cfg, &ctx);

		/* Synchronize context with changes made by the physics engine */

		Time += dt;
		ctx.Time = Time;
		fprintf(stdout, "\nT= %.4f My", Time/Matosec);

		if (switch_ps>=2) Write_Ouput(&cfg, &ctx);
	} while (Time < Timefinal-dt/10);

	The_End(&cfg, &ctx);

	return(1);
}



/**************************************************/
/****	  SUBROUTINES  IN  RANDOM  ORDER	 *****/
/**************************************************/



int inputs(ModelConfig *cfg, ModelContext *ctx, int argc, char **argv)
{
	int 	reformat=0 ;
	char 	resume_filename[MAXLENFILE], 
		command[MAXLENLINE], 
		load_file_name[MAXLENLINE];
	bool	success_def_prm=false, switch_initial_geom=false;

	run_type = 0;
	nmax_input_points = 5000;
	switch_strs_history = true;
	setbuf(stdout, NULL);

	putenv("tao_dir=" TAODIR); 
	
	/*Version of tAo will be matched against the parameters file *.PRM*/
	/*¡¡ UPDATE template.PRM !!*/
	strncpy(version, VERSION, sizeof(version) - 1);
	version[sizeof(version) - 1] = '\0';

	/*Default parameter values are read from ./tao/doc/template.PRM*/
	snprintf(projectname, sizeof(projectname), "%s/doc/template", TAODIR);
	success_def_prm = read_file_parameters(0, 0);
	projectname[0] = '\0';

	for (int iarg=1; iarg<argc; iarg++) {
		if (argv[iarg][0] == '-') {
			float 	value;
			char 	prm[MAXLENLINE];
			for (int ilet=2; ilet < strlen(argv[iarg])+2; ilet++) 
				prm[ilet-2] = argv[iarg][ilet];
			value=atof(prm);
			switch (argv[iarg][1]) {
				case 'f':
					reformat=1;
					if (argv[iarg][2]) reformat = value;
					break;
				case 'F':
					run_type=2;
					if (strlen(prm)>0) {
						strncpy(resume_filename, prm, sizeof(resume_filename) - 1);
						resume_filename[sizeof(resume_filename) - 1] = '\0';
					}
					else snprintf(resume_filename, sizeof(resume_filename), "%s.all", projectname);
					break;
				case 'h':
					switch (argv[iarg][2]) {
						case 'p':
							fprintf(stderr, "\nFile ./tao/doc/template.PRM (sample parameters file) follows in stdout:\n") ;
						snprintf(command, sizeof(command), "cat %s/doc/template.PRM", TAODIR);
						system(command) ;
						break;
						case 'c':
							fprintf(stderr, "\nFile ./tao/doc/template.PRM (sample parameters file) follows in stdout:\n") ;
						snprintf(command, sizeof(command), "cat %s/doc/template.PRM | %s/script/cleanPRM", TAODIR, TAODIR);
						system(command) ;
						break;
						case 'u':
							fprintf(stderr, "\nFile ./tao/doc/template.UNIT (sample unit file) follows in stdout:\n") ;
						snprintf(command, sizeof(command), "cat %s/doc/template.UNIT", TAODIR);
						system(command) ;
						break;
						default:
						syntax();
						AUTHORSHIP;
						break;
					}
					fprintf(stderr, "\n") ;
					exit(0);
				case '-':
					if (strcmp(argv[iarg], "--help") == 0) {
						fprintf(stderr, "\nFile ./tao/doc/tao.info.txt follows:\n") ;
						snprintf(command, sizeof(command), "more %s/doc/tao.info.txt", TAODIR);
						system(command) ;
						exit(0);
					}
					break;
				case 'Q':
					run_type=1;
					strncpy(load_file_name, prm, sizeof(load_file_name) - 1);
					load_file_name[sizeof(load_file_name) - 1] = '\0';
					break;
				case 'V':
					verbose_level = 1;
					if (argv[iarg][2]) verbose_level = value;
					break;
			}
		}
		else {
			if (run_type != 2) run_type=10;
			if (strlen(projectname)<1) {
				strncpy(projectname, argv[iarg], sizeof(projectname) - 1);
				projectname[sizeof(projectname) - 1] = '\0';
			}
		}
	}

	if (verbose_level>=1) {
	  fprintf(stdout, 	
	  	"\n**** tAo: 2D LITHOSPHERIC FLEXURE, EROSION/SEDIMENTATION, AND FORELAND BASIN MODELING ****"
	  	"\nVersion: %s", version);
	  fflush(stdout);
	}

	if (!run_type) {
		syntax();
		fprintf(stdout, "\n\nType %s -h for further information.\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	nloads=0; n_image=0; nlakes=0;
	numBlocks=0; i_first_Block_load=0; i_Block_insert=0;
	nwrotenfiles=0; switch_topoest=false;

	switch (run_type)
	{
		case 0:
			fprintf(stdout, "\n\n*** END of run *** \n\n");
			exit(0);
		case 1:
			interpr_command_line_opts(argc, argv);
			Direct_mode(cfg, ctx, load_file_name);
			exit(0);
		case 2:
			read_file_resume(resume_filename);
			interpr_command_line_opts(argc, argv);
			if (verbose_level>=1) fprintf(stdout, "\nResuming project '%s'. Timefinal=%.1f My", projectname, Timefinal/Matosec);
			if (switch_ps>=2) n_image--; /*Don't produce 2 jpg's of the same stage of restart*/
			return(1);
		case 10:
			if (!read_file_parameters(verbose_level>=1, 0)) {
				syntax();
				fprintf(stderr, "\n\tAvailable parameter files in this directory:\n");	
				system("ls *.PRM");
				if (!success_def_prm) {
					PRINT_ERROR("\t\aDefault parameters file './tao/doc/template.PRM' could not be read.\n"); 
				}
				exit(EXIT_FAILURE);
			}
			if (reformat) {
				snprintf(projectname, sizeof(projectname), "%s/doc/template", TAODIR);
				read_file_parameters(0, reformat);
				exit(0);
			}
			interpr_command_line_opts(argc, argv);
			
			cfg->Nx = Nx; cfg->Nz = Nz; cfg->dx = dx; cfg->dz = dz;
			cfg->x0 = x0; cfg->xf = xf; cfg->xmin = xmin; cfg->xmax = xmax;
			cfg->hydro_model = hydro_model;
			cfg->erosed_model = erosed_model;
			cfg->verbose_level = verbose_level;
			cfg->densenv = densenv; cfg->denssedim = denssedim;
			
			Allocate_Memory(cfg, ctx);
			break;
	}

	{
		char filename[MAXLENLINE]; FILE *file;
		snprintf(filename, sizeof(filename), "%s.out", projectname);
		if (switch_file_out) {
			if ((file = fopen(filename, "w")) == NULL) {
				PRINT_ERROR("Cannot open standard output file %s.\n", filename);
			}
			else {
				if (verbose_level>=1) fprintf(stdout, "\nInfo: standard output redirected to %s.\n", filename);
			stdout=file;
			}
		}
		else {
		remove(filename);
		}
	}

	dx = (xf-x0) / (Nx-1);
	dz = (100000) / (Nz-1);
	dt *= Matosec;
	tau *= Matosec; 
	dt_record *= Matosec;
	Timefinal *= Matosec;
	Timeini *= Matosec;
	Time = Timeini;
	last_time_file_time = Timeini - 9999*dt_record;	/*very old*/
	dt_eros *= Matosec;
	Kerosdif /= secsperyr;
	Keroseol /= Matosec;
	Ksedim   /= Matosec;
	rain /= secsperyr*1e3;
	if (hydro_model==1) Krain *= 1e6/Matosec/1e3/1e3;
	evaporation_ct /= secsperyr*1e3;
	lost_rate *= 1e-2 * 1e-3;
	temp_sea_level += TEMP_FREEZE_WATER;
	switch_write_file_Blocks = 1;


	if (verbose_level>=3) {
		fprintf(stdout, "\nPlate model: %d \tBC=%d\tNx=%d from %.2f to %.2f (dx=%f m)", isost_model, boundary_conds, Nx, x0, xf, dx);
		fprintf(stdout, "\nDensities: asthenosphere =%f; mantle=%f; crust=%.2f", densasthen, densmantle, denscrust);
		fprintf(stdout, "\nTiming: from %.2f to %.2f every %.2f My\n", Timeini/Matosec, Timefinal/Matosec, dt/Matosec);
	}


	/*Test of incompatibilities between parameters*/
	if (densenv && water_load)		{ water_load=0; 	PRINT_WARNING("Sea not possible when densenv<>0. Sea switch turned off.") ; }
	if (!erosed_model && (Ksedim || Kerosdif || Keroseol))	{ Ksedim=Kerosdif=Keroseol=0 ; if (verbose_level>=1) PRINT_WARNING("Erosion-sedimentation unswitched. Parameters Ksedim, Kerosdif & Keroseol have no effect.") ; }
	if (!water_load && Ksedim)		{ Ksedim=0;			PRINT_WARNING("Warning: Sea presence isn't switched, Ksedim has no effect.");}
	if (!erosed_model)  			{ Ksedim=Kerosdif=Keroseol=0 ; }
	if (!hydro_model && erosed_model>1)  	{ erosed_model=1 ; PRINT_WARNING("Warning: hydro_model is 0; erosed_model is switched to 1 accordingly.") ; }
	if (switch_ps && !switch_write_file)	{ if (verbose_level>=3)	PRINT_WARNING("Warning: switch_write_file needed to make a postscript. Postscript may not be done.") ; }
	if (isost_model!=2)  			{ tau=0; }
	if (tau<=0 && isost_model==2)  		{ isost_model=1; }
	if (boundary_conds == 0)		{ appmoment = 0; }
	if (boundary_conds == 2)		{ vert_force=0; appmoment = 0; }


	snprintf(command, sizeof(command), "rm -f %s.temp0 %s.mtrz %s.grv_mod", projectname, projectname, projectname);
	system(command);

	read_file_sea_level(); 
	ctx->sea_level = calculate_sea_level(ctx->Time);
	read_file_horiz_record_time();
	read_file_Te();
	read_file_Crust_Thick(crust_thick_default);
	read_file_Upper_Crust_Thick(upper_crust_thick_default);	
	read_file_YSE(cfg);
	
	cfg->Nx = Nx; cfg->Nz = Nz; cfg->dx = dx; cfg->dz = dz;
	cfg->x0 = x0; cfg->xf = xf; cfg->xmin = xmin; cfg->xmax = xmax;
	cfg->hydro_model = hydro_model;
	cfg->erosed_model = erosed_model;
	cfg->verbose_level = verbose_level;
	cfg->denscrust = denscrust;

	read_file_Temperature(cfg, ctx) ;
	Init_Stress(cfg);

	switch_initial_geom = read_file_initial_deflection(w) + read_file_initial_topo(topo) ;
	for (int i=0; i<cfg->Nx; i++)  { 
		ctx->topo[i] += random_topo * ((((float) rand()) / ((float) RAND_MAX)) -.5);
		ctx->topo[i] += zini;
		Blocks_base[i] = ctx->topo[i];
		ctx->topo[i] -= w[i];
		h_water[i] = MAX_2(ctx->sea_level-ctx->topo[i], 0);
	}
	if (switch_initial_geom && verbose_level>=1) {
		float	altmax=-1e9, altmin=1e9;
		for (int i=0; i<cfg->Nx; i++) {
			if (altmax<ctx->topo[i])	altmax=ctx->topo[i];
			if (altmin>ctx->topo[i])	altmin=ctx->topo[i];
		}
		fprintf(stdout, "\n  alt.Init.:  max = %9.1f m	 min = %9.1f m   ", altmax, altmin);
	}

	return(1);
}



int interpr_command_line_opts(int argc, char **argv) 
{
	/*Interpretates the command line options described in tao.info.txt*/

	PRINT_INFO("Enetering command line interpretation.");
	for (int iarg=1; iarg<argc; iarg++) {
		if (argv[iarg][0] == '-') {
			float 	value, value2;
			char 	prm[MAXLENLINE], prm2[MAXLENLINE], *ptr;
			for (int ilet=2; ilet < strlen(argv[iarg])+2; ilet++) 
				prm[ilet-2] = argv[iarg][ilet];
			for (int ilet=3; ilet < strlen(argv[iarg])+2; ilet++) 
				prm2[ilet-3] = argv[iarg][ilet];
			value  = atof(prm);
			value2 = atof(prm2);
			PRINT_DEBUG("\aArgument: %s", argv[iarg]);
			switch (argv[iarg][1]) {
				case 'A':
					grav_anom_type = value;
					break;
				case 'B':
					boundary_conds = value;
					break;
				case 'D':
					if (run_type!=2) {
						x0 = atof(strtok(prm, "/"));
						xf = atof(strtok(NULL, "/"));
					}
					if (xmin<x0) xmin=x0;
					if (xmax>xf) xmax=xf;
					else if (run_type!=1) fprintf(stdout, "\nWarning: Impossible to change the domain when resuming a model.");
					break;
				case 'd':
					if (run_type!=2) Nx = (int) (xf-x0)/value + 1;
					else if (run_type!=1) fprintf(stdout, "\nWarning: Imposible to change dx or Nx value when resuming a model.");
					break;
				case 'L': /*OLD*/
					if (run_type!=2) x0 = value;
					else if (run_type!=1) fprintf(stdout, "\nWarning: Imposible to change x0 value when resuming a model.");
					if (xmin<x0) xmin=x0;
					break;
				case 'M':
					isost_model = value;
					if (prm[1] == 't') switch_strs_history = false;
					break;
				case 'm':
					appmoment = value;
					break;
				case 'N':
					Nx = value;
					break;
				case 'o':
					switch_file_out=true;
					break;
				case 'P':
					switch_ps=(strlen(prm)>0)?value:1; /*default is 1 to keep old command line syntax with -Pc*/
					switch_write_file_Blocks=true;
					if (argv[iarg][2] == 'c') {
						switch_ps=2;
						gif_geom[0] = '\0';
						if (strlen(prm2)>0) {
							strncpy(gif_geom, prm2, sizeof(gif_geom) - 1);
							gif_geom[sizeof(gif_geom) - 1] = '\0';
						}
					}
					if (argv[iarg][2] == 'p') {
						switch_ps=3;
						gif_geom[0] = '\0';
						if (strlen(prm2)>0) {
							strncpy(gif_geom, prm2, sizeof(gif_geom) - 1);
							gif_geom[sizeof(gif_geom) - 1] = '\0';
						}
					}
					break;
				case 'p':
					horz_force = value;
					break;
				case 'q':
					ptr = strtok(prm, "=");
					match_parameter(ptr, strtok(NULL, "/"), 1, 0, prm);
					break;
				case 'R': /*OLD*/
					if (run_type!=2) xf = value;
					else fprintf(stdout, "\nWarning: Imposible to change xf value when restarting a model.");
					if (xmax>xf) xmax=xf;
					break;
				case 'r':
					switch (argv[iarg][2]) {
						case 'e':	densenv	= value2; 	break;
						case 'c':	denscrust  = value2; 	break;
						case 'i':	densinfill = value2; 	break;
						case 'm':	densmantle = value2; 	break;
						case 'a':	densasthen = value2; 	break;
					}
					break;
				case 'S':
					{
						int iblock, nblocks;
						struct BLOCK_1D	Block_aux;
						iblock = atoi(strtok(prm, "/"));
						nblocks = atoi(strtok(NULL, "/"));
						PRINT_INFO("Block %d will be moved by %d positions", iblock, nblocks);
						Block_aux=Blocks[iblock];
						if (nblocks>0) {
							for (int iu=iblock; iu<iblock+nblocks; iu++) {
								Blocks[iu]=Blocks[iu+1];
								PRINT_INFO("%d = %d", iu, iu+1);
							}
						}
						else {
							for (int iu=iblock; iu>iblock+nblocks; iu--) {
								Blocks[iu]=Blocks[iu-1];
								PRINT_INFO("%d = %d", iu, iu-1);
							}
						}
						Blocks[iblock+nblocks]=Block_aux;
						PRINT_INFO("%d = %d", iblock+nblocks, numBlocks-1);
					}
					break;
				case 's':
					vert_force = value;
					break;
				case 'T':
					Te_default = value;
					break;
				case 't':
					switch (argv[iarg][2]) {
						case 'i':	Timeini	= value2; 	if (run_type==2) Timeini *= Matosec;	break;
						case 'f':	Timefinal  = value2;if (run_type==2) Timefinal *= Matosec;	break;
						case 'd':	dt		 = value2; 	if (run_type==2) dt *= Matosec;			break;
						case 'e':	dt_eros	= value2; 	if (run_type==2) dt_eros *= Matosec;	break;
						case 'v':	tau		= value2; 	if (run_type==2) tau *= Matosec; 		break;
						case 'r':	dt_record = value2; if (run_type==2) dt_record *= Matosec;	break;
					}
					break;
				case 'V':
					verbose_level = 1;
					if (argv[iarg][2]) verbose_level = value;
					break;
				case 'v':
					{
						float density, velocity;
						density = atof(strtok(prm, "/"));
						velocity = atof(strtok(NULL, "/"));
						for (int iu=0; iu<numBlocks; iu++) {
							if (Blocks[iu].density==-density || iu==density) {
								Blocks[iu].vel=velocity*1e3/Matosec;
								Blocks[iu].last_vel_time=Time-dt;/*!!*/
								Blocks[iu].last_shift=0;
							}
						}
					}
					break;
			}
		}
	}
	return(1);
}



int Direct_mode(ModelConfig *cfg, ModelContext *ctx, char *load_file_name)
{
	int 	i;
	FILE 	*file;

	/*Solves flexure problem in direct mode: taking a single external load file and 
	writing deflection in standar output. There are no other input files neither 
	output files*/

	PRINT_INFO("Entering direct mode. x0,xf,Nx= %.1f,%.1f,%d\n", x0,xf, Nx);
	dx = (xf-x0) / (Nx-1) ;
	cfg->Nx = Nx; cfg->Nz = Nz; cfg->dx = dx; cfg->dz = dz;
	cfg->x0 = x0; cfg->xf = xf; cfg->xmin = xmin; cfg->xmax = xmax;
	cfg->densenv = densenv; cfg->denscrust = denscrust; cfg->densmantle = densmantle;
	ctx->Time = Time; ctx->Timeini = Timeini; ctx->dt = dt;

	Allocate_Memory(cfg, ctx);
	if (strcmp(load_file_name, "")) {
		if ((file = fopen(load_file_name, "rt")) == NULL) {
			fprintf(stderr, "\nLoad file '%s' not found.\n", load_file_name);
			exit(0);
		}
		readinterplin(file, h_last_unit, cfg->Nx, cfg->x0, cfg->xf) ;
		fclose(file);
	}
	for (i=0; i<cfg->Nx; i++) {
		D[i] = ET2RIG(Te_default); 
		Dq[i] = h_last_unit[i];
	}
	Elastoplastic_Deflection(cfg, ctx);
	fprintf(stdout, "\n\n#x[km]\t\tw[m]\t\tpressure[Pa]\n"); 
	for (i=0; i<cfg->Nx; i++) 
		fprintf(stdout, "%8.1f\t%8.1f\t%8.1f\n", 
			(cfg->x0+i*cfg->dx)/1e3, w[i], h_last_unit[i]);
	fprintf(stdout, "\n"); 
	return 1;
}




int tectload(ModelConfig *cfg, ModelContext *ctx)
{
	/*
	CALCULATES NEW LOAD INCREMENT FROM UNIT FILES, Returns 1 if elastic
	flexure must be done (i.e, if changes in load  occurred), 0 otherwise.
	*/

	PRINT_ARRAY_INFO_1D(topo, "topogr.", "m", "m2") 

	/*Reads external load from file(s)*/
	while (read_file_unit(cfg, ctx));

	/*Moves Blocks*/
	move_Blocks(cfg, ctx);

	/*Interpolates loads through time*/
	gradual_Block(cfg, ctx);

	Repare_Blocks(cfg, ctx);

	return (1);
}




int Elastoplastic_Deflection(ModelConfig *cfg, ModelContext *ctx)
{
	int 	i;
	bool	load_changes=false;

	/*
	CALLS SUBROUTINES TO SOLVE ELASTIC OR ELASTIC-PLASTIC FLEXURE WITH 
	THE NEW LOAD INCREMENT AND ADDS NEW DEFLECTION INCREMENT TO THE 
	TOTAL DEFLECTION.
	*/
	
	for (i=0; i<cfg->Nx; i++) if (Dq[i]) load_changes = true;
	if (isost_model>0 && (load_changes || (ctx->Time==ctx->Timeini && (horz_force || vert_force || appmoment)))) {
		fprintf(stdout, " e");	  fflush(stdout);
		if (isost_model<3) {
			if (!Te_default) {
				/*LOCAL ISOSTASY*/
				float Krest;
				for (i=0; i<cfg->Nx; i++) {
					GET_KREST(Krest, q, i)
					Dw[i] = Dq[i] / Krest;
				}
			}
			else {
				/*REGIONAL ISOSTASY*/
				/*Pure Elastic Flexure without rheological inputs*/
				float 	*moment;
				double	**A,		/*Linear System Matrix (diagonal terms)*/
					*b;		/*Independent Column*/
				int 	NDs=3, NDi=3;
					A = alloc_matrix_dbl (cfg->Nx, NDi+1+NDs);
					b = (double *) calloc (cfg->Nx , sizeof(double));
					moment = (float *) calloc (cfg->Nx , sizeof(float));
					LES_matrix(cfg, ctx, A, b, D, q, Dq, w, false) ;
					solveLES(A, b, cfg->Nx, NDs, NDi, Dw) ;
					free_matrix_dbl(A, cfg->Nx);
				free(b); 

					for (i=0;i<cfg->Nx;i++) {
					w[i] += Dw[i] ;
						if (i != 0 && i != cfg->Nx-1) 
						moment[i] += 
							-D[i] * (Dw[i-1] - 2*Dw[i] + Dw[i+1]) 
								/ pow(cfg->dx,2) ;
					}
					if (isost_model!=2) flexural_stats(cfg, ctx, moment);
				free(moment);

			}
		}
		else {
			/*Elasto-plastic flexure with EET rheological calculation*/
			Rheo_Flex_Iter(cfg, ctx);
		}

		if (switch_topoest) {
	  		/*Defines the thickness of last infill Block*/
	  		for (i=0; i<cfg->Nx; i++)
	  			Blocks[i_first_Block_load-1].thick[i] +=  MAX_2(Dw[i], 0) ;
		}
	}

	PRINT_ARRAY_INFO_1D(Dq, "load_incrm", "N/m2", "N/m") 

	/*Adds all this load to the total load array*/
	for (int i=0; i<cfg->Nx; i++)  {
		q[i] += Dq[i];
		if (Dq[i]) load_changes = true;
	}
	if (load_changes) fprintf(stdout, " d");
	/*Resets deflection and load arrays*/
	for (i=0;i<cfg->Nx;i++)  Dq[i]=Dw[i]=0;

	fflush(stdout);
	return(1);
}



int surface_processes(ModelConfig *cfg, ModelContext *ctx)
{
	/* 
	CALCULATES EROSION AND SEDIMENTATION:
	*/
	float 	eros_level;
	bool	switch_horiz_record=false;

	PRINT_DEBUG("erosed_model=%d hydro_model=%d", erosed_model, hydro_model);

	total_sed_mass=total_bedrock_eros_mass=0;

	if (!erosed_model && !hydro_model) return (0);
	switch_topoest=false;
	

	/*Creates a new sediment Block if necessary*/
	if (erosed_model) {
		int i;
		float TimelastBlock=-9999*Matosec;
		for (int i=0; i<Nx; i++) eros_now[i]=0;
		for (int i=0; i<numBlocks; i++) 
		if (Blocks[i].age > TimelastBlock && Blocks[i].density==denssedim) TimelastBlock = Blocks[i].age;
		for (int i=0; i<n_record_times; i++) 
		if (Time>horiz_record_time[i]-dt/2 && Time<=horiz_record_time[i]+dt/2) 
			switch_horiz_record=true;
		if (Time == Timeini 
		  || ((Time-TimelastBlock)>(dt_record-.001*dt) && dt_record && !n_record_times)
		  || switch_horiz_record) {
			insert_new_Block(cfg, ctx, ctx->numBlocks);
			Blocks[ctx->numBlocks-1].type = 'S' ;
			Blocks[ctx->numBlocks-1].density = cfg->denssedim ;
			Blocks[ctx->numBlocks-1].erodibility = erodibility_sed ;
			Blocks[ctx->numBlocks-1].detr_ratio = calloc(cfg->Nx, sizeof(float));
			Blocks[ctx->numBlocks-1].detr_grsize = calloc(cfg->Nx, sizeof(float));
		}
	}
	/*Calculates height at every point*/
	calculate_topo(cfg, ctx, topo);

	/*Diffusive Erosion*/
	Diffusive_Eros_1D (cfg, ctx, Kerosdif);

	/*Adds background erosion and sea sedimentation*/
	constant_rate_eros (cfg, ctx, Keroseol, Ksedim, water_load, n_eros_level_input_points, var_eros_level, &eros_level);

	/*Adds fluvial transport*/
	Surface_Transport (cfg, ctx, topo);


	/*Calculates water column load*/
	calculate_water_load(cfg, ctx);



	if (n_eros_level_input_points)
		fprintf(stdout,
			"\teros_lvl.: %8.1f m", eros_level);

	if (cfg->verbose_level>=1) {
		if (cfg->erosed_model>1) {
			int i, i_biggest_nosea=0, n_biggest_nosea=0;
			double error_water;
			for (i=1; i<=ctx->nlakes; i++) {
				if (Lake[i].n > n_biggest_nosea) {
					if (Lake[i].n_sd) {
					if (ctx->topo[Lake[i].sd[0]] > ctx->sea_level)
						n_biggest_nosea = Lake[i].n; i_biggest_nosea = i;
					}
					else	n_biggest_nosea = Lake[i].n; i_biggest_nosea = i;
				}
			}
			for (i=1; i<=ctx->nlakes; i++) {
				float vol=0;
				if (i==i_biggest_nosea || (Lake[i].n>ceil((double) cfg->Nx/100) && cfg->verbose_level>=3) || cfg->verbose_level>=4) {
				PRINT_SUMLINE("lake %d/%d: %6.2f km2 %6.1f km%5.0f m ", 
					i, ctx->nlakes, Lake[i].vol/1e6, Lake[i].n*cfg->dx/1e3, Lake[i].alt);
				/*write the position of the lowest node*/
				if (Lake[i].n) fprintf(stdout, "bott@ %.0f km", (Lake[i].cell[0]*cfg->dx+cfg->x0)/1e3);
				if (Lake[i].n_sd) {
					if (ctx->topo[Lake[i].sd[0]]>ctx->sea_level) {
					 fprintf(stdout, " %d out@ %.0f km %8.2e m3/s", Lake[i].n_sd, (Lake[i].sd[0]*cfg->dx+cfg->x0)/1e3, drainage[Lake[i].sd[0]].discharge/*Lake_Input_Discharge(i)*/);
					}
					else fprintf(stdout, " Sea	   %8.2e m3/s", Lake_Input_Discharge(cfg, i));
				}
				else	 fprintf(stdout, " Endorheic %8.2e m3/s", Lake_Input_Discharge(cfg, i));
				}
			}
		}
	}
	if (cfg->verbose_level>=1) {
		int i, imaxerosion=SIGNAL, imaxsediment=SIGNAL, imaxaltsediment=SIGNAL;
		float top_Block, totalerosion=0, volsediment=0, maxerosion=0, maxsediment=0, max_altit_seds=-9999;
		for (i=(cfg->xmin-cfg->x0)/cfg->dx; i<cfg->Nx-(cfg->xf-cfg->xmax)/cfg->dx; i++) {
			totalerosion  -= MASS2SEDTHICK_1D(cfg, total_erosion[i])*cfg->dx;
			if (maxerosion<total_erosion[i])
				{maxerosion = total_erosion[i]; imaxerosion=i;}
			if (maxsediment<-total_erosion[i])
				{maxsediment = -total_erosion[i]; imaxsediment=i;}
			top_Block = Blocks_base[i]-w[i];
			for (int j=0; j<ctx->numBlocks; j++) {
				top_Block+=Blocks[j].thick[i];
				if (Blocks[j].density==cfg->denssedim) volsediment+=Blocks[j].thick[i]*cfg->dx;
				if (Blocks[j].density==cfg->denssedim && top_Block>max_altit_seds && Blocks[j].thick[i]>1) {
					max_altit_seds = top_Block; imaxaltsediment=i;
				}
			}
		}
		PRINT_SUMLINE("total_cumul_eros-sed= %8.1f km2\tvol.sediment= %8.1f km2",
			totalerosion/1e6, volsediment/1e6
		);
		if (imaxerosion!=SIGNAL) fprintf(stdout,
			"\n  erosion_max.   = %8.1f m  \t@x= %5.1f km", 
			MASS2SEDTHICK_1D(cfg, maxerosion),  (cfg->x0+imaxerosion*cfg->dx)/1000
		);
		if (imaxsediment!=SIGNAL) fprintf(stdout,
			"\n  sediment.max.  = %8.1f m  \t@x= %5.1f km" 
			"\n  max.altit.sed. = %8.1f m  \t@x= %5.1f km", 
			MASS2SEDTHICK_1D(cfg, maxsediment), (cfg->x0+imaxsediment*cfg->dx)/1000, 
			max_altit_seds, (cfg->x0+imaxaltsediment*cfg->dx)/1000
		);
	}

	fflush(stdout);
	return(1);
}




int gravity_anomaly(ModelConfig *cfg, ModelContext *ctx)
{
	/*
	USES gravanompolyg() TO CALCULATE GRAVITY ANOMALY ALONG HORIZONS PROFILE 
			(results in m/s2) 
	USES geoidanompolyg() TO CALCULATE GEOID ANOMALY 
			(results in m) 
		z-values relative to sea level. This means that zero 
		deflection is at an elevation of 'zini'.
	*/

	register int   	i, ix, k_Block, i_Block, np_aux_Block, ix_min, ix_max;
	float 	z_max_grav_model = -zini+crust_thick_default+40000, 
		*Block_aux_x, *Block_aux_z, 
		*upper_hori_aux, *lower_hori_aux, 
		*alt_measurement,
		aux;
	float	*gravanom, 		/*Gravity anomaly	[mgal]*/
		*geoidanom;		/*Geoid anomaly		[m]   */
	FILE 	*file;
	char 	filename[MAXLENFILE];

	sprintf(filename,"%s.xg", projectname);
	remove(filename);

	if (!grav_anom_type) return(0);

	fprintf(stdout, " a");	fflush(stdout);

	geoidanom =  (float *) calloc(cfg->Nx, sizeof(float));
	gravanom =   (float *) calloc(cfg->Nx, sizeof(float));
	Block_aux_x = (float *) calloc (cfg->Nx*2+5, sizeof(float));
	Block_aux_z = (float *) calloc (cfg->Nx*2+5, sizeof(float));
	alt_measurement = (float *) calloc (cfg->Nx, sizeof(float));
	lower_hori_aux = (float *) calloc (cfg->Nx, sizeof(float));
	upper_hori_aux = (float *) calloc (cfg->Nx, sizeof(float));
	for (ix=0; ix<cfg->Nx; ix++) {gravanom[ix] = 0;  geoidanom[ix] = 0;}

	ix_min = MAX_2((cfg->xmin-cfg->x0-.1*cfg->dx)/cfg->dx, 0) ;	ix_max = MIN_2(floor((cfg->xmax-cfg->x0+.1*cfg->dx)/cfg->dx) + 2, cfg->Nx);

	Repare_Blocks(cfg, ctx);

	sprintf(filename,"%s.grv_mod", projectname);
	remove(filename);
	if ((file = fopen(filename, "wt")) == NULL) {
		if (cfg->verbose_level>=3) fprintf(stderr, "Warning: Cannot open output file '%s'.\n", filename);
	}
	fprintf(file, "#Time: %.2fMy\tGravity model bodies: x(km)-z(m)", ctx->Time/Matosec);

	/*Calculates topography at every point*/
	calculate_topo(cfg, ctx, ctx->topo);
	for (ix=0; ix<cfg->Nx; ix++)  alt_measurement[ix] = MAX_2(ctx->topo[ix], 0) + 10;

	/*Calculates gravity atraction of the Blocks*/
	np_aux_Block=cfg->Nx*2+4;
	for (i_Block=0; i_Block<ctx->numBlocks; i_Block++) {
		/*Constructs the auxiliar Block*/
		for (ix=0; ix<cfg->Nx; ix++) {
			upper_hori_aux[ix] = w[ix]-Blocks_base[ix];
			for (k_Block=0; k_Block<=i_Block; k_Block++) upper_hori_aux[ix] -= Blocks[k_Block].thick[ix];
			lower_hori_aux[ix] = w[ix]-Blocks_base[ix];
			for (k_Block=0; k_Block<i_Block; k_Block++)  lower_hori_aux[ix] -= Blocks[k_Block].thick[ix];
		}
		make_gravi_body (cfg, upper_hori_aux, lower_hori_aux, Block_aux_x, Block_aux_z);
		fprintf(file, "\n>");
		for (i=0; i<np_aux_Block; i++) fprintf(file, "\n%9.2f %9.1f", Block_aux_x[i]/1000, -Block_aux_z[i]);

		/*Calculates anomaly for this Block*/
		for (ix=ix_min; ix<ix_max; ix++) {
			gravanom[ix]  += gravanompolyg (Block_aux_x, Block_aux_z, np_aux_Block, cfg->x0+ix*cfg->dx, -alt_measurement[ix], Blocks[i_Block].density);
			geoidanom[ix] += geoidanompolyg(Block_aux_x, Block_aux_z, np_aux_Block, cfg->x0+ix*cfg->dx, -10e3, Blocks[i_Block].density);
//PRINT_ERROR("\ngrav %.2e    geoid %.2e   iblock=%d", gravanom[ix], geoidanom[ix], i_Block);
		}
	}

	/*Calculates gravity atraction of sea water*/
	if (water_load) {
		for (ix=0; ix<cfg->Nx; ix++) {
			upper_hori_aux[ix] = 0;
			lower_hori_aux[ix] = MAX_2(-ctx->topo[ix], 0);
		}
		make_gravi_body (cfg, upper_hori_aux, lower_hori_aux, Block_aux_x, Block_aux_z);
		fprintf(file, "\n>");
		for (i=0; i<np_aux_Block; i++) fprintf(file, "\n%9.2f %9.1f", Block_aux_x[i]/1000, -Block_aux_z[i]);

		for (ix=ix_min; ix<ix_max; ix++) {
			gravanom[ix]  += gravanompolyg (Block_aux_x, Block_aux_z, np_aux_Block, cfg->x0+ix*cfg->dx, -alt_measurement[ix], cfg->denswater);
			geoidanom[ix] += geoidanompolyg(Block_aux_x, Block_aux_z, np_aux_Block, cfg->x0+ix*cfg->dx, -10e3, cfg->denswater);
		}
	}

	/*Calculates gravity atraction of crust*/
	for (ix=0; ix<cfg->Nx; ix++) {
		upper_hori_aux[ix] = w[ix]-Blocks_base[ix];
		lower_hori_aux[ix] = w[ix]+crust_thick[ix];
	}
	make_gravi_body (cfg, upper_hori_aux, lower_hori_aux, Block_aux_x, Block_aux_z);
	fprintf(file, "\n>");
	for (i=0; i<np_aux_Block; i++) fprintf(file, "\n%9.2f %9.1f", Block_aux_x[i]/1000, -Block_aux_z[i]);

	for (ix=ix_min; ix<ix_max; ix++) {
		gravanom[ix]  += gravanompolyg (Block_aux_x, Block_aux_z, np_aux_Block, cfg->x0+ix*cfg->dx, -alt_measurement[ix], cfg->denscrust);
		aux=geoidanompolyg(Block_aux_x, Block_aux_z, np_aux_Block, cfg->x0+ix*cfg->dx, -10e3, cfg->denscrust);
		geoidanom[ix] += aux;
	}

	/*Calculates gravity atraction of mantle*/
	for (ix=0; ix<cfg->Nx; ix++) {
		upper_hori_aux[ix] = w[ix]+crust_thick[ix] - zini;
		lower_hori_aux[ix] = z_max_grav_model;
	}
	make_gravi_body (cfg, upper_hori_aux, lower_hori_aux, Block_aux_x, Block_aux_z);
	fprintf(file, "\n>");
	for (i=0; i<np_aux_Block; i++) fprintf(file, "\n%9.2f %9.1f", Block_aux_x[i]/1000, -Block_aux_z[i]);

	for (ix=ix_min; ix<ix_max; ix++) {
		aux = gravanompolyg (Block_aux_x, Block_aux_z, np_aux_Block, cfg->x0+ix*cfg->dx, -alt_measurement[ix], cfg->densmantle);
		gravanom[ix]  += aux;
		geoidanom[ix] -= geoidanompolyg(Block_aux_x, Block_aux_z, np_aux_Block, cfg->x0+ix*cfg->dx, -10e3, cfg->densmantle);
	}

	/*Substracts the atraction of normal water layer*/
	if (water_load && zini<0) {
		Block_aux_x[0] = x0 - 1e8;	Block_aux_z[0] = - zini;
		Block_aux_x[1] = x0 - 1e8;	Block_aux_z[1] = 0;
		Block_aux_x[2] = xf + 1e8;	Block_aux_z[2] = 0;
		Block_aux_x[3] = xf + 1e8;	Block_aux_z[3] = - zini;
		np_aux_Block = 4;
		fprintf(file, "\n>");
		for (i=0; i<np_aux_Block; i++) fprintf(file, "\n%9.2f %9.1f", Block_aux_x[i]/1000, -Block_aux_z[i]);

		for (ix=ix_min; ix<ix_max; ix++) {
			gravanom[ix]  -= gravanompolyg (Block_aux_x, Block_aux_z, np_aux_Block, cfg->x0+ix*cfg->dx, -alt_measurement[ix], cfg->denswater);
			geoidanom[ix] -= geoidanompolyg(Block_aux_x, Block_aux_z, np_aux_Block, cfg->x0+ix*cfg->dx, -10e3, cfg->denswater);
		}
	}

	/*Substracts the atraction of normal crust*/
	Block_aux_x[0] = cfg->x0 - 1e8;	Block_aux_z[0] = - zini + crust_thick_default;
	Block_aux_x[1] = cfg->x0 - 1e8;	Block_aux_z[1] = - zini;
	Block_aux_x[2] = cfg->xf + 1e8;	Block_aux_z[2] = - zini;
	Block_aux_x[3] = cfg->xf + 1e8;	Block_aux_z[3] = - zini + crust_thick_default;
	np_aux_Block = 4;
	fprintf(file, "\n>");
	for (i=0; i<np_aux_Block; i++) {fprintf(file, "\n%9.2f %9.1f", Block_aux_x[i]/1000, -Block_aux_z[i]);}

	for (ix=ix_min; ix<ix_max; ix++) {
		aux = gravanompolyg (Block_aux_x, Block_aux_z, np_aux_Block, cfg->x0+ix*cfg->dx, -alt_measurement[ix], cfg->denscrust);
		gravanom[ix]  -= aux;
		geoidanom[ix] -= geoidanompolyg(Block_aux_x, Block_aux_z, np_aux_Block, cfg->x0+ix*cfg->dx, -10e3, cfg->denscrust);
	}

	/*Substracts the atraction of normal mantle*/
	Block_aux_x[0] = cfg->x0 - 1e8;	Block_aux_z[0] = z_max_grav_model;
	Block_aux_x[1] = cfg->x0 - 1e8;	Block_aux_z[1] = - zini + crust_thick_default;
	Block_aux_x[2] = cfg->xf + 1e8;	Block_aux_z[2] = - zini + crust_thick_default;
	Block_aux_x[3] = cfg->xf + 1e8;	Block_aux_z[3] = z_max_grav_model;
	np_aux_Block = 4;
	fprintf(file, "\n>");
	for (i=0; i<np_aux_Block; i++) fprintf(file, "\n%9.2f %9.1f", Block_aux_x[i]/1000, -Block_aux_z[i]);

	for (ix=ix_min; ix<ix_max; ix++) {
		aux = gravanompolyg (Block_aux_x, Block_aux_z, np_aux_Block, cfg->x0+ix*cfg->dx, -alt_measurement[ix], cfg->densmantle);
		gravanom[ix]  -= aux;
		geoidanom[ix] -= geoidanompolyg(Block_aux_x, Block_aux_z, np_aux_Block, cfg->x0+ix*cfg->dx, -10e3, cfg->densmantle);
	}


	/* Now we add Bouguer correction:
		Free-air correction is not applied since in a 2D model there are no
		vertical variations of the regional gravity field related to the 
		distance to the center of earth (Turcotte & Schubert, 1982).
	*/
	if (grav_anom_type==1) {
		for (ix=0; ix<cfg->Nx; ix++) if (cfg->x0+ix*cfg->dx >= cfg->xmin && cfg->x0+ix*cfg->dx <= cfg->xmax) {
			/*Assumed standar 2.67 g/cm3*/
			if (ctx->topo[ix]>ctx->sea_level) 
				gravanom[ix] -= 2*pi * 2670 * CGU * ctx->topo[ix];
			if (ctx->topo[ix]<ctx->sea_level) 
				gravanom[ix] -= 2*pi*(2670 - ((water_load)? 1000:0))* CGU * ctx->topo[ix];
		}
	}

	write_file_grav_anom(cfg, ctx, gravanom, geoidanom);

	free(Block_aux_x);
	free(Block_aux_z);
	free(alt_measurement);
	free(upper_hori_aux);
	free(lower_hori_aux);
	fprintf(file, "\n>");
	fclose(file);

	return(1);
}




int move_Blocks(ModelConfig *cfg, ModelContext *ctx)
{
	int	*nshift;
	float	*new_thick;

	/*
	  Moves the Blocks and calculates the isostatic load and thickness change
	  Deforms the sediment.
	*/

	new_thick = calloc(cfg->Nx, sizeof(float));
	nshift = calloc(ctx->numBlocks, sizeof(int));

	/*CRAWL UP along all Blocks to find the relevant moving ones first*/
	for (int iu=0; iu<ctx->numBlocks; iu++) {
		PRINT_DEBUG("move_Block:  Block %d; nshift=%d", iu, nshift[iu]);
		if (Blocks[iu].density == cfg->denssedim) {
		/*DEFORM SEDIMENT UNITS*/
		for (int i=0; i<cfg->Nx; i++) new_thick[i] = Blocks[iu].thick[i];
		for (int i=0; i<cfg->Nx; i++) {
			float sedthick;
			/*CRAWL DOWN Blocks to find the uppermost moving Block below this point*/
			for (int ju=iu-1, sedthick=Blocks[iu].thick[i]; ju>=0; ju--) {
			/*Calculate the thickness of sediments between the top of this sed. Block and the moving Block*/
			if (Blocks[ju].density == cfg->denssedim) {
				sedthick += Blocks[ju].thick[i];
			}
			else {
   			  /*Amount of cells to propagate the deformation: ~20 deg assumed.*/
			  int nprop = SIGN(nshift[ju]) * (int) ceil(sedthick*5/cfg->dx);
			  int i_unprop = i-nprop;
			  DOMAIN_LIMIT_1D(i_unprop);
			  if (Blocks[ju].thick[i_unprop]>.1) {
				if (!nshift[ju]) 
					break;
				else {
					int i_shift =   i+nshift[ju];
					/*If block ju is moving below [i] then shift seds.*/
					if (deform_sed && IN_DOMAIN_1D(i_shift)) 
					new_thick[i_shift] += Blocks[iu].thick[i];
					if (deform_sed && IN_DOMAIN_1D(i))  
					new_thick[i]	   -= Blocks[iu].thick[i];
					break;
				}
			  }
				}
			}
		}
		}
		else 
		{
		/*MOVE BLOCK UNITS and define nshift[]*/
		float theor_shift = Blocks[iu].vel * (ctx->Time-Blocks[iu].last_vel_time);
		nshift[iu] = floor((theor_shift - Blocks[iu].last_shift) /cfg->dx +.5);
		if (ctx->Time > Blocks[iu].time_stop + .1*ctx->dt) {nshift[iu]=0;}
		Blocks[iu].shift += nshift[iu]*cfg->dx;
		Blocks[iu].last_shift += nshift[iu]*cfg->dx;
		for (int i=0; i<cfg->Nx; i++) {
			int i_unshifted = i-nshift[iu];	
			DOMAIN_LIMIT_1D(i_unshifted);
			new_thick[i] = Blocks[iu].thick[i_unshifted];
		}
		}
		for (int i=0; i<cfg->Nx; i++) {
			if (new_thick[i]<0) PRINT_ERROR("\aBlock %d has a negative thickness: %.2f m", iu, new_thick[i]);
				Dq[i] += g * (new_thick[i] - Blocks[iu].thick[i]) * Blocks[iu].density;
				Blocks[iu].thick[i] = new_thick[i];
		}
	}
	fflush(stdout);
	free(new_thick);
	free(nshift);
	return(1);
}




int read_file_unit(ModelConfig *cfg, ModelContext *ctx)
{
	/*
	  READS UNIT FILE NAMED 'projectnameNUM.UNIT' WHERE 'NUM' IS 1 FOR THE
	  FIRST UNIT, 2 FOR THE SECOND, ETC. Interpolates this input.
	  Creates new Block to store its properties and cuts sediment
	  Blocks when file contains fault depth rather than a thickness itself.
	*/

	int 	cut_Block, nparams=0;
	float	time_stop=9999/*My*/, time_unit, 
		erodibility_aux=NO_DATA, fill_up_to=NO_DATA, 
		vel=0, density=NO_DATA;
	bool 	insert, cut_seds, cut_Blocks, cut_all, top, fault, switch_move,
		ride, hidden, z_absol;
	FILE 	*file;
	char 	filename[MAXLENFILE];

	/*Read the next unit age*/
	sprintf(filename, "%s%d.UNIT", projectname, nloads+1);
	if ((file = fopen(filename, "rt")) == NULL) {
		PRINT_WARNING("Cannot read unit file '%s'.", filename);
		return (0);
	}
	time_unit = ctx->Timeini/Matosec;
	{
		int nlines=0, nread, show, replace=0;
		char str1[MAXLENLINE], str2[MAXLENLINE], line[MAXLENLINE+200], *lineptr;
		show=(cfg->verbose_level>=3)? 1 : 0;
		rewind(file);
		while ((lineptr=fgets(line, MAXLENLINE+200-1, file)) != NULL && nlines<NMAXHEADERLINES) {
		nlines++; nread=sscanf(lineptr, "%s %s", str1, str2);
		if (nread == 2) {
			Match_Param_Replace_flt ( "time",  	time_unit, 0 )
			/*Old versions:*/
		}
		}
	}
	time_unit *= Matosec;
	/*Return if it isn't time yet to read the new unit file*/
	if (time_unit>ctx->Time+.1*ctx->dt || time_unit<ctx->Timeini) return(0);

	PRINT_INFO("Reading '%s'", filename);
	switch_move = fault = switch_gradual = false; // Reset flags for each unit file
	insert = hidden = cut_seds = cut_Blocks = cut_all = false; // Reset flags
	top = ride = z_absol = false; // Reset flags
	i_Block_insert = ctx->numBlocks;
	cut_Block = 0;

	/*READS AND INTERPOLATES UNIT/LOAD FILE*/
	{
		int nlines=0, nread, show, replace=0;
		char str1[MAXLENLINE], str2[MAXLENLINE], line[MAXLENLINE+200], *lineptr;
		rewind(file); 
		while ((lineptr=fgets(line, MAXLENLINE+200-1, file)) != NULL && nlines<NMAXHEADERLINES) {
				nlines++; nread=sscanf(lineptr, "%s %s", str1, str2);
				if (nread == 2) {
			Match_Param_Replace_flt ( "vel",		vel, 	0 )
			Match_Param_Replace_flt ( "time_stop",  	time_stop, 	0 )
			Match_Param_Replace_flt ( "density", 	 	density, 	0 )
			Match_Param_Replace_flt ( "erodibility",	erodibility_aux, 	0 )
			Match_Param_Replace_int ( "gradual",		switch_gradual, 	0 )
			Match_Param_Replace_int ( "hidden",		hidden, 	0 )
			Match_Param_Replace_int ( "ride",		ride, 	0 )
			Match_Param_Replace_int ( "insert",		insert, 	0 )
			Match_Param_Replace_int ( "top",		top, 	0 )
			Match_Param_Replace_int ( "move",  		switch_move, 	0 )
			Match_Param_Replace_int ( "fault",	  	fault, 	0 )
			Match_Param_Replace_int ( "z_absol",	  	z_absol,   	0 )
			Match_Param_Replace_int ( "cut_seds",  		cut_seds, 	0 )
			Match_Param_Replace_int ( "cut_Block",  		cut_Block, 	0 )
			Match_Param_Replace_int ( "cut_Blocks",  	cut_Blocks, 	0 )
			Match_Param_Replace_int ( "cut_all",  		cut_all, 	0 )
			Match_Param_Replace_int ( "topoest",		switch_topoest, 	0 )
			Match_Param_Replace_flt ( "fill_up_to", 	fill_up_to, 	0 )
			/*Old versions:*/
			Match_Param_Replace_int ( "cut_unit",  		cut_Block, 	0 )
			Match_Param_Replace_int ( "cut_units",  	cut_Blocks, 	0 )
			Match_Param_Replace_flt ( "erodability",	erodibility_aux, 	1 )
				}
				if (strcmp(str1, "thickness_distribution")==0) break;
		}
		rewind(file); 
	}
	if (fill_up_to == NO_DATA) 
		readinterplin(file, h_last_unit, cfg->Nx, cfg->x0, cfg->xf) ;
	else {
		for (int i=0; i<cfg->Nx; i++) h_last_unit[i] = MAX_2(0, fill_up_to-ctx->topo[i]);
	}
	fclose(file);

	nloads++;

	/*Check incompatibilities between unit file signals*/
	if (switch_gradual && switch_move) {
		PRINT_WARNING("Gradual+moving units not implemented. It won't be gradual.");
		switch_gradual = false;
	}

	vel *= 1e3/Matosec;
	time_stop *= Matosec;

	/*ACT ACCORDING TO THE SIGNALS*/
	if (fault) {
		switch_move = true;
	}
	/*Creates a Block of infill if switch_topoest; it will be filled during Deflection*/
	if (switch_topoest) {
		insert_new_Block(cfg, ctx, i_first_Block_load);
		Blocks[i_first_Block_load].type = 'I'; /*stands for Infill*/
		Blocks[i_first_Block_load].density = cfg->densinfill;
		if (cfg->densinfill<2550) Blocks[i_first_Block_load].erodibility = erodibility_sed;
		i_first_Block_load++; 	i_Block_insert++;
	}
	if (insert) {
		i_Block_insert = 0;
	}
	if (top) {
		for (int k=ctx->numBlocks-1; k>=0; k--) {
			if (Blocks[k].density != cfg->denssedim) {
				i_Block_insert = k+1;
				break;
			}
		}
	}
	if (cut_all) {
		cut_Blocks = true;
	}

	if (fault && !top && !cut_all) i_Block_insert = 0;

	if (fault && !cut_Block) {
		int numBlocks0=ctx->numBlocks;
		/*Make a copy of all Blocks*/
		PRINT_DEBUG("Cutting Blocks: numBlocks= %d", ctx->numBlocks);
		for (int k=0; k<numBlocks0; k++) {
			float *thick_aux;
			insert_new_Block(cfg, ctx, ctx->numBlocks);
			thick_aux = Blocks[ctx->numBlocks-1].thick;
			Blocks[ctx->numBlocks-1] = Blocks[k];
			Blocks[ctx->numBlocks-1].thick = thick_aux;
			Blocks[ctx->numBlocks-1].vel = vel;
			Blocks[ctx->numBlocks-1].last_vel_time = ctx->Time-ctx->dt;/*!!*/
			Blocks[ctx->numBlocks-1].last_shift = 0;
			Blocks[ctx->numBlocks-1].time_stop = time_stop;
			if (Blocks[ctx->numBlocks-1].type == 'S') {
				Blocks[ctx->numBlocks-1].detr_ratio  = calloc(cfg->Nx, sizeof(float));
				Blocks[ctx->numBlocks-1].detr_grsize = calloc(cfg->Nx, sizeof(float));
			}
			if (density		 != NO_DATA && Blocks[ctx->numBlocks-1].type != 'S') Blocks[ctx->numBlocks-1].density	 = density;
			if (erodibility_aux != NO_DATA && Blocks[ctx->numBlocks-1].type != 'S') Blocks[ctx->numBlocks-1].erodibility = erodibility_aux;
		}
		/*Modify Blocks_base and cut above the fault*/
		for (int i=0; i<cfg->Nx; i++) {
			float z_fault = -h_last_unit[i], base_of_Block = Blocks_base[i];
			if (z_absol) base_of_Block -= w[i];
			h_last_unit[i] = MAX_2(0, Blocks_base[i] - z_fault);
			Blocks_base[i]  = MIN_2(Blocks_base[i], z_fault);
			if (cut_Blocks) {
				for (int k=0; k<numBlocks0; k++) {
						float top_of_Block=base_of_Block+Blocks[k].thick[i];
					if (Blocks[k].density == cfg->denssedim && !cut_all) {
						break;
					}
					if (z_fault <= base_of_Block) {
						Blocks[k+numBlocks0].thick[i] += Blocks[k].thick[i];
						Blocks[k].thick[i]			= 0;
					}
					else {
						Blocks[k+numBlocks0].thick[i] += MAX_2 (0, top_of_Block-z_fault);
						Blocks[k].thick[i]		   -= MAX_2 (0, top_of_Block-z_fault);
					}
						base_of_Block = top_of_Block;
				}
			}
		}
	}
	if (density		==NO_DATA) density		 = cfg->denscrust;
	if (erodibility_aux==NO_DATA) erodibility_aux = erodibility;

	/*Create a new Block for this file*/
	insert_new_Block(cfg, ctx, i_Block_insert);

	/*Add the thickness in file to the new Block; Thin Blocks and basement if the thickness is negative*/
	if (!switch_gradual && !hidden && !cut_Block) {
		for (int i=0; i<cfg->Nx; i++) {
		if (h_last_unit[i]>0) 
			Blocks[i_Block_insert].thick[i] = h_last_unit[i];
		else {
			float   h_load_aux, h_load_aux2;
			int	 j;
			h_load_aux = fabs(h_last_unit[i]);
			for (j=i_Block_insert-1; h_load_aux>0 && j>=0; j--) {
					h_load_aux2 = MIN_2(Blocks[j].thick[i], h_load_aux);
					h_load_aux -= h_load_aux2;
					Blocks[j].thick[i] -= h_load_aux2;
			}
			/*j is the deepest eroded Block in i*/
			if (j==-1) {
					Blocks_base[i] -= h_load_aux;
			}
		}
		}
	}
	if (hidden) Blocks[i_Block_insert].type = 'H';
	Blocks[i_Block_insert].density = density;
	Blocks[i_Block_insert].erodibility = erodibility_aux;
	Blocks[i_Block_insert].vel = vel;
	Blocks[i_Block_insert].time_stop = time_stop;

	if (ride) {
		PRINT_DEBUG("ride %d %d", i_Block_insert, ctx->numBlocks);
		for (int i_Block=i_Block_insert+1; i_Block<ctx->numBlocks; i_Block++) {
			Blocks[i_Block].vel		   = Blocks[i_Block_insert].vel; 
			Blocks[i_Block].last_shift	= 0; 
			Blocks[i_Block].last_vel_time = ctx->Time-ctx->dt;/*!!*/ 
			Blocks[i_Block].time_stop	 = Blocks[i_Block_insert].time_stop; 
		}
	}

	if (cut_Block) {
		int numBlocks0=ctx->numBlocks;
		float *thick_aux;
		if (cut_Block<0) {
			/*Search the biggest Block with density close to -cut_Block*/
			float vol, volmax=0, Blockvolmax=-1;
			for (int k=0; k<ctx->numBlocks; k++) {
				for (int i=vol=0; i<cfg->Nx; i++)  vol += Blocks[k].thick[i];
				if (vol>volmax && fabs(Blocks[k].density-fabs(cut_Block))<.5) {volmax=vol; Blockvolmax=k;}
			}
			cut_Block = Blockvolmax+1;
		}
		PRINT_INFO("Cutting Block %d", cut_Block);
		insert_new_Block(cfg, ctx, cut_Block);
		thick_aux = Blocks[cut_Block].thick;
		Blocks[cut_Block] = Blocks[cut_Block-1];
		Blocks[cut_Block].thick = thick_aux;
		Blocks[cut_Block].vel = vel;
		Blocks[cut_Block].density = density;
		Blocks[cut_Block].last_vel_time = ctx->Time-ctx->dt;/*!!*/
		Blocks[cut_Block].last_shift = 0;
		Blocks[cut_Block].time_stop = time_stop;
		/*Cut the Block*/
		for (int i=0; i<cfg->Nx; i++) {
			float z_fault = -h_last_unit[i], base_of_Block = Blocks_base[i], top_of_Block;
			if (z_absol) base_of_Block -= w[i];
			for (int k=0; k<cut_Block-1; k++) base_of_Block += Blocks[k].thick[i];
			top_of_Block = base_of_Block + Blocks[cut_Block-1].thick[i];
			if (z_fault <= base_of_Block) {
				Blocks[cut_Block].thick[i]   += Blocks[cut_Block-1].thick[i];
				Blocks[cut_Block-1].thick[i]  = 0;
			}
			else {
				Blocks[cut_Block].thick[i]   += MAX_2 (0, top_of_Block-z_fault);
				Blocks[cut_Block-1].thick[i] -= MAX_2 (0, top_of_Block-z_fault);
			}
		}
	}

	/*Don't Repare_Blocks() in case of: 
		Gradual load, because then h_last_unit[] will be empty until tectload()
		Topoest load, because the infill Block will be filled upon deflection.
	*/
	if (!switch_gradual && !switch_topoest) Repare_Blocks(cfg, ctx);

	/*Increment the isostatic load for this time interval*/
	if (!switch_gradual && !fault && !cut_Block) /*!cut_Block condition needed for Alice's setup (redefining densities without 'fault' option*/
		for (int i=0; i<cfg->Nx; i++) Dq[i] += (density-cfg->densenv)*g*h_last_unit[i];

	PRINT_INFO("Unit read from '%s'. ", filename);
	PRINT_DEBUG("%d params; dens=%.0f kg/m3; erodibility=%.1e; ", nparams, density, erodibility_aux);
	if (switch_gradual) PRINT_INFO("Will be gradually loaded until %.2fMy.", time_stop/Matosec);
	if (switch_move) PRINT_INFO("Vel= %.2f km/My till T=%.1f My", vel*Matosec/1000, time_stop/Matosec);

	fprintf(stdout, " l") ;
	if (fault)  		fprintf(stdout, "F") ;
	if (hidden)		fprintf(stdout, "H") ;
	if (insert) 		fprintf(stdout, "I") ;
	if (top) 		fprintf(stdout, "P") ;
	if (switch_gradual)	fprintf(stdout, "G") ;
	if (switch_move) 	fprintf(stdout, "M") ;
	if (ride) 		fprintf(stdout, "R") ;
	if (switch_topoest)	fprintf(stdout, "T") ;

	return(1);
}




int syntax()
{
	/*
		Displays the hardcoded command line syntax of the program
	*/
	fprintf(stderr, "\nSyntax:\n");
	fprintf(stderr, "  tao  project  -A[1|2] -B<bound_type> -D[x0/xf] -d<dx> -F[file] -f[2] \n");
	fprintf(stderr, "        -h[i|u|p|c] -M<lih_type>[t] -m<app_mom> -N<Nx> -o -P[c[geom]|p]\n");
	fprintf(stderr, "        -p<tec_force> -q<param=value> -Q<file>  -r<a|c|i|m|a><density>\n");
	fprintf(stderr, "        -S<b>/<n> -s<app_force> -T<eet> -t<i|f|d|v|r><time> -V[<level>] \n");
	fprintf(stderr, "        -v[<num>/<vel>]\n\n");
	fprintf(stderr, "  Options:\n");
	fprintf(stderr, "    'project'\tRoot name for all files (e.g. 'test' looks for 'test.PRM')\n");
	fprintf(stderr, "    -F\t\tResume a previous run from a .all file\n");
	fprintf(stderr, "    -P\t\tProduce graphic output (-Pp for Python, -Pc for GMT)\n");
	fprintf(stderr, "    -q\t\tOverride any parameter in the PRM file (-qparam=value)\n");
	fprintf(stderr, "    -Q\t\tDirect elastic deflection mode for a given load file\n");
	fprintf(stderr, "    -h\t\tShow this help (or -hc for clean PRM, -hu for UNIT example)\n");
	fprintf(stderr, "    --help\tOpen the full documentation file\n\n");
	fprintf(stderr, "For full documentation, please read doc/tAo_Documentation.md\n");
	
	return (1);
}





int The_End(ModelConfig *cfg, ModelContext *ctx)
{
	int	i, j;
	char	command[MAXLENLINE];	
	float	total_load=0, total_restitutive_force=0, total_Blocks_mass=0, total_sed_mass=0, total_sed_grain_mass=0, 
		x, xleft, xright, Krest;

	Krest =  (switch_topoest) ?  (densasthen-densinfill)*g : (densasthen-densenv)*g ;

	if (grav_anom_type) gravity_anomaly(cfg, ctx);
	if (isost_model>=3 && !switch_YSE_file) write_file_Temperature(cfg, ctx);

	if (verbose_level>=1) {
		fprintf(stdout, "\n\nFinal statistics:");
		for (i=0; i<cfg->Nx; i++)  {total_load += q[i]; total_restitutive_force += (Krest*w[i]);}
		total_load *= (cfg->dx*cfg->Nx/(cfg->Nx+1));
		total_restitutive_force *= (cfg->dx*cfg->Nx/(cfg->Nx+1));

		for (i=0; i<ctx->numBlocks; i++)  for (j=0; j<cfg->Nx; j++) {
			total_Blocks_mass += Blocks[i].thick[j] * Blocks[i].density;
			if (Blocks[i].density == denssedim) {
				total_sed_mass += Blocks[i].thick[j] * Blocks[i].density;
				total_sed_grain_mass += THICK2SEDMASS_1D(cfg, Blocks[i].thick[j]);
			}
		}
		total_Blocks_mass *= (cfg->dx*cfg->Nx/(cfg->Nx+1));
		total_sed_mass  *= (cfg->dx*cfg->Nx/(cfg->Nx+1));
		total_sed_grain_mass  *= ((float) cfg->Nx/ (float) (cfg->Nx+1));
		fprintf(stdout, "\n\tTotal_load_weight = %10.3e N/m", total_load);
		fprintf(stdout, "\n\tTotal_rest._force = %10.3e N/m", total_restitutive_force);
		fprintf(stdout, "\n\tTotal_Blocks_mass = %10.3e kg/m", total_Blocks_mass);
		fprintf(stdout, "\n\tTotal_sedim_mass  = %10.3e kg/m (%.3g kg/m of grain)", total_sed_mass, total_sed_grain_mass);
	}

	fprintf(stdout, "\n\n%d Blocks:", ctx->numBlocks);
	fprintf(stdout, "\nNo.\tDensity\tAge \tVolume \tVel   \tShift \tFrom x\tTo x \tStop\ttype");
	if (erosed_model) fprintf(stdout, "\tErodabil.");
	fprintf(stdout, "\n   \t(kg/m3)\t(My)\t(km2) \t(km/My)\t(km)  \t(km)  \t(km) \t (My)\t(-)");
	if (erosed_model) fprintf(stdout, "\t		 ");
	for (i=ctx->numBlocks-1; i>=0; i--) {
		float volume;
		for (xleft=cfg->x0,j=0; j<cfg->Nx; j++) {
			x=cfg->x0+j*cfg->dx; 
			if (Blocks[i].thick[j]>1) break;
			xleft=x;
		}
		for (xright=cfg->xf,j=cfg->Nx-1; j>=0; j--) {
			x=cfg->x0+j*cfg->dx; 
			if (Blocks[i].thick[j]>1) break;
			xright=x;
		}
		for (volume=j=0; j<cfg->Nx; j++) {
			volume += Blocks[i].thick[j];
		}
		volume *= cfg->dx*cfg->Nx/(cfg->Nx+1);
		fprintf(stdout, "\n%2d:\t%.0f \t%.1f\t%.1f\t%.2f\t%.1f\t%.1f\t%.1f\t%.1f\t%c ", 
			i, Blocks[i].density, Blocks[i].age/Matosec, volume/1e6, 
			Blocks[i].vel/1e3*Matosec, Blocks[i].shift/1e3, xleft/1e3, xright/1e3, Blocks[i].time_stop/Matosec, Blocks[i].type); 
		if (erosed_model) fprintf(stdout, "\t%.2e", Blocks[i].erodibility);
	}
	fprintf(stdout, "\n -:\t%.0f\t%.1f\t-\t0\t-\t-\t-\t-\tbasement\n", denscrust, Timeini/Matosec);

	if (switch_write_file_Blocks) write_file_resume(cfg, ctx);
	if (switch_ps<2) Write_Ouput(cfg, ctx);

	sprintf(command, "rm -f %s*.tao.tmp", projectname);
	system(command);
	{
		char filename[MAXLENFILE];
		FILE *file;
		sprintf(filename, "%s/.taodefaults", getenv ("HOME"));
		if ((file = fopen(filename, "rt")) == NULL) {
			sprintf(command, "echo First use of %s by `whoami` at `hostname` > %s; date >> %s", version, filename, filename);
			if (verbose_level>=3) fprintf(stdout, "\n%s", command);
			system(command);
			sprintf(command, "mail danielgc@ictja.csic.es < %s", filename);
			if (verbose_level>=3) fprintf(stdout, "\n%s", command);
			system(command);
			fclose (file);
		}
	}
	if (verbose_level>=3) AUTHORSHIP;
	fprintf(stdout, "\n\n");

	if (run_type==10 || run_type==2) exit(0);
	return(1);
}




int Viscous_Relaxation(ModelConfig *cfg, ModelContext *ctx)
{
	int 	i, NDs=3, NDi=3;
	double	**A, *b;
	float	*dwdt, *moment;

	/*
	CALL SUBROUTINES TO CALCULATE VISCOELASTIC RELAXATION FROM PRESENT 
	DEFLECTION AND TOTAL LOAD.   ADD NEW DEFLECTION INCREMENT TO THE 
	PREVIOUS DEFLECTION
	*/

	if (isost_model!=2 || !Te_default)  return(0);

	A = alloc_matrix_dbl (cfg->Nx, NDi+1+NDs);
	b = (double *) calloc (cfg->Nx , sizeof(double));
	moment = (float *) calloc (cfg->Nx , sizeof(float));
	dwdt   = (float *) calloc (cfg->Nx , sizeof(float));

	fprintf(stdout, " v");
	LES_matrix(cfg, ctx, A, b, D, q, Dq, w, true) ;
	solveLES(A, b, cfg->Nx, NDs, NDi, dwdt) ;
	for (i=0; i<cfg->Nx; i++) {
		w[i] += dwdt[i]*ctx->dt;	Dw[i] = dwdt[i]*ctx->dt;
	}
	for (i=1; i<cfg->Nx-1; i++) {
		moment[i] = (- D[i]*(1-nu*nu)/(1-.25)  * tau * 
				(dwdt[i+1] - 2*dwdt[i] + dwdt[i-1]) / 
				pow(cfg->dx, 2) + moment[i]*tau/ctx->dt ) / 
				(1 + tau/ctx->dt) ;
		/*If the deflection was made constant then:*/
		/*  	/= exp(dt/tau) ;   */
	}
	moment[0]=moment[1];
	moment[cfg->Nx-1]=moment[cfg->Nx-2];

	if (switch_topoest) {
		/*Defines the thickness of last infill Block*/
		for (i=0; i<cfg->Nx; i++) 
			Blocks[i_first_Block_load-1].thick[i] +=  MAX_2(Dw[i], 0) ;
	}

	flexural_stats(cfg, ctx, moment);

	free(moment);
	free_matrix_dbl(A, cfg->Nx);
	free(b);
	free(dwdt);
	return(1);
}



int Write_Ouput(ModelConfig *cfg, ModelContext *ctx)
{
	if (isost_model) write_file_time(cfg, ctx, w, ctx->topo);

	if (switch_write_file_Blocks) write_file_Blocks(cfg, ctx);
	if (cfg->erosed_model) write_file_erosed(cfg, ctx, total_erosion);
	if (isost_model>=3) write_file_Te(cfg, ctx);
	if (isost_model>=3) write_file_stress(cfg, ctx);
	if (isost_model>=3) write_file_maxmompoint(cfg, ctx);

	/*Make GMT Postscript*/
	if (switch_ps) {
		char 	command[300];
		if (switch_ps<=2) {
			sprintf(command, "tao.gmt.job %s %.2f %.2f %.2f %.2f %d %d", 
				projectname, xmin/1000, xmax/1000, zmin, zmax, water_load, (isost_model<3)? 0:1);
			if (verbose_level>=3) 
				fprintf(stdout, "\nPostscript file '%s.ps' is going to be produced with command:\n%s\n", projectname, command);
			system(command);
		}
		else {
			sprintf(command, "tao.plot.py %s; mv -f %s.png %s%03d.png", projectname, projectname, projectname, n_image);
			if (verbose_level>=3) 
				fprintf(stdout, "\nPlot files '%s.xvg' and %s%03d.png to be produced with command:\n%s\n", projectname, projectname, n_image, command) ;
			system(command);
			n_image++;
		}
		if (switch_ps==2) {
			/*crop by default to the border*/
			if (strlen(gif_geom)<2) sprintf(gif_geom, "-trim -background Khaki -label 'tAo software: %s' -gravity South -append", projectname);
			sprintf(command, "magick -density 200 %s.ps %s -interlace NONE  %s%03d.jpg", /*-fill \"#ffffff\" -draw \"rectangle 70,10 130,25\" -fill \"#000000\" -font helvetica -draw \"text 74,22 t_%+3.2f_My \" */
				projectname, gif_geom, projectname, n_image);
			if (verbose_level>=3)
				fprintf(stdout, "\n%s\n", command) ;
			system(command);
			n_image++;
		}
	}

	return (1);
}
