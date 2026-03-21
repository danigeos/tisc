/*******************************************************************************
*****                              TISC                                    *****
********************************************************************************
	For compilation and installation check the file tisc/README.md
*******************************************************************************/


#include "tisc.h"
#include "tisclib.h"
#include "tiscio.h"
#include "libreria.h"
#include <string.h>

int main(int argc, char **argv)
{
	float **topo_ant;
	ModelConfig cfg = {0};
	ModelContext ctx = {0};

	/*get input parameters and files*/
	inputs(argc, argv) ;

	fprintf(stdout, "\nT= %.4f My", Time/Matosec);

	topo_ant = alloc_matrix(Ny, Nx);

	/* Populate configuration context */
	cfg.Nx = Nx; cfg.Ny = Ny;
	cfg.dx = dx; cfg.dy = dy; cfg.dxy = dxy;
	cfg.xmin = xmin; cfg.xmax = xmax; cfg.ymin = ymin; cfg.ymax = ymax;
	cfg.hydro_model = hydro_model;
	cfg.erosed_model = erosed_model;
	cfg.verbose_level = verbose_level;
	cfg.densenv = densenv; cfg.denssedim = denssedim;
	cfg.denscrust = denscrust; cfg.denswater = denswater; cfg.densice = densice;
	cfg.Px = Px; cfg.Py = Py; cfg.Pxy = Pxy;

	ctx.Timeini = Timeini; ctx.Timefinal = Timefinal;
	ctx.dt_eros = dt_eros;
	ctx.topo = topo;
	ctx.numBlocks = numBlocks;
	ctx.Time = Time;
	ctx.dt = dt;
	ctx.sea_level = sea_level;
	ctx.nlakes = nlakes;
	ctx.total_accum_gypsum = total_accum_gypsum;
	ctx.total_accum_halite = total_accum_halite;
	ctx.total_precip_gypsum_rate = total_precip_gypsum_rate;
	ctx.total_precip_halite_rate = total_precip_halite_rate;

	float **compaction_prev = alloc_matrix(cfg.Ny, cfg.Nx);
	calculate_compaction(&cfg, &ctx, compaction_prev);

	if (plotting_mode>=2) {calculate_topo(&cfg, &ctx, topo); Write_Ouput(&cfg, &ctx);}

	/*MAIN LOOP: In this loop time increments from Timeini to Timefinal*/
	do {
		/*Remember topography before tectonics and flexure*/
		calculate_topo(&cfg, &ctx, topo_ant);

		/*Calculate tectonic deformation and tectonic load*/
		tectload(&cfg, &ctx);
		
		/*Sea level variations*/
		ctx.sea_level = calculate_sea_level(ctx.Time);

		/*Calculates water column load*/
		calculate_water_load(&cfg, &ctx);

		/* Adjust Dq for compaction water loss */
		float **compaction_now = alloc_matrix(cfg.Ny, cfg.Nx);
		calculate_compaction(&cfg, &ctx, compaction_now);
		for (int i=0; i<cfg.Ny; i++) {
			for (int j=0; j<cfg.Nx; j++) {
				Dq[i][j] -= (compaction_now[i][j] - compaction_prev[i][j]) * cfg.denswater * g;
				compaction_prev[i][j] = compaction_now[i][j];
			}
		}
		free_matrix(compaction_now, cfg.Ny);

		/*Define & solve elastic flexure equation*/
		Elastic_Deflection(&cfg, &ctx);

		/*Define & solve viscoelastic flexure equation*/
		Viscous_Relaxation(&cfg, &ctx);

		/* Update dynamic context */
		ctx.Time = Time;
		ctx.dt = dt;
		ctx.sea_level = sea_level;
		ctx.nlakes = nlakes;
		ctx.numBlocks = numBlocks;

		/*Calculates surface water flow and sediment-erosion. Distributes all previous effects on topo over ninters substeps of dt_st.*/
		surface_processes(topo_ant, &cfg, &ctx);

		/* Synchronize context with changes made by the physics engine */
		ctx.nlakes = nlakes;
		ctx.numBlocks = numBlocks;

		Time += dt;
		ctx.Time = Time;
		if (verbose_level>=2) fprintf(stdout, "\n--------------------------------------------------------------------------");
		fprintf(stdout, "\nT= %.4f My", Time/Matosec);

		if (plotting_mode>=2) Write_Ouput(&cfg, &ctx);
	} while (Time < Timefinal-dt/10);

	The_End(&cfg, &ctx);

	free_matrix (topo_ant, Ny);
	free_matrix (compaction_prev, Ny);

	return(1);
}




/**************************************************/
/****      SUBROUTINES  IN  RANDOM  ORDER     *****/
/**************************************************/


int inputs (int argc, char **argv)
{
	int	i, j, iarg, reformat=0 ;
	char 	command[MAXLENLINE], 
		resume_filename[MAXLENFILE],
		load_file_name[MAXLENLINE];
	bool	success_def_prm=false, switch_initial_geom=false;

	run_type=0;
	solver_type = 'l';
	setbuf(stdout, NULL);

	putenv("tisc_dir=" TISCDIR); /*printf("\nEnvir. Read: %s", getenv ("tisc_dir")); system ("echo test: $tisc_dir"); system ("printenv | grep tisc");*/

	/*Version of TISC is matched against the parameters file *.PRM*/
	/*¡¡ UPDATE template.PRM !!*/
	//defined in tisc/config.mk
	strncpy(version, VERSION, sizeof(version) - 1);
	version[sizeof(version) - 1] = '\0';

	/*Default parameter values are read from ./tisc/doc/template.PRM:*/
	snprintf(projectname, sizeof(projectname), "%s/tisc/doc/template", TISCDIR);
	success_def_prm = read_file_parameters(0, 0);
	projectname[0] = '\0';

	for (iarg=1; iarg<argc; iarg++) {
		if (argv[iarg][0] == '-') {
			int ilet;
			float 	value;
			char 	prm[MAXLENLINE];
			for (ilet=2; ilet < strlen(argv[iarg])+2; ilet++) 
				prm[ilet-2] = argv[iarg][ilet];
			value  = atof(prm);
			switch (argv[iarg][1]) {
				case 'f':
					reformat=1;
					if (argv[iarg][2]) reformat = (int)value;
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
							fprintf(stderr, "\nFile ./tisc/doc/template.PRM (sample parameters file) follows in stdout:\n") ;
							snprintf(command, sizeof(command), "cat %s/tisc/doc/template.PRM", TISCDIR);
							system(command) ;
							break;
						case 'c':
							fprintf(stderr, "\nFile ./tisc/doc/template.PRM (sample parameters file) follows in stdout:\n") ;
							// The cleanPRM script might not exist or be needed anymore. If it's just for display, cat is enough.
							// If it's a real script, ensure it's in PATH. For now, assuming it's a script that processes the PRM.
							snprintf(command, sizeof(command), "cat %s/tisc/doc/template.PRM | %s/tisc/script/cleanPRM", TISCDIR, TISCDIR);
							system(command) ;
							break;
						case 'u':
							fprintf(stderr, "\nFile ./tisc/doc/template.UNIT (sample unit file) follows in stdout:\n") ;
							snprintf(command, sizeof(command), "cat %s/tisc/doc/template.UNIT", TISCDIR);
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
						fprintf(stderr, "\nFile ./tisc/doc/tisc.info.txt follows:\n") ;
						snprintf(command, sizeof(command), "more %s/tisc/doc/tisc.info.txt", TISCDIR);
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
					if (argv[iarg][2]) verbose_level = (int)value;
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
	  fprintf(stdout, "\n* TISC: TECTONICS, ISOSTASY, SURFACE PROCESSES, AND CLIMATE PLANFORM MODELING *");
	  fprintf(stdout, "\nVersion: %s", version) ;
	  fflush(stdout);
	}

	if (!run_type) {
		syntax();
		fprintf(stdout, "\nType %s -h for further information.\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	nloads=0; n_image=0; nlakes=0;
	numBlocks=0; i_first_Block_load=0; i_Block_insert=0;
	nwrotenfiles=0; switch_topoest=false;

	switch (run_type)
	{
		case 0:
			fprintf(stdout, "\n\n\t*** END of run ***\n\n\n\n");
			exit(0);
		case 1:
			interpr_command_line_opts(argc, argv);
			Direct_mode(load_file_name);
			exit(0);
		case 2:
			read_file_resume(resume_filename);
			interpr_command_line_opts(argc, argv);
			if (verbose_level>=1) fprintf(stdout, "\nResuming project '%s'. Timefinal=%.1f My", projectname, Timefinal/Matosec);
			return(1);
		case 10:
			if (!read_file_parameters(verbose_level>=1, 0)) {
				syntax();
				fprintf(stderr, "\n\tAvailable parameter files in this directory:\n");
				system("ls *.PRM");
				if (!success_def_prm) {
					PRINT_ERROR("\aDefault parameters file './tisc/doc/template.PRM' could not be read.\n"); 
				}
				exit(EXIT_FAILURE);
			}
			interpr_command_line_opts(argc, argv);
			if (reformat) {
				snprintf(projectname, sizeof(projectname), "%s/tisc/doc/template", TISCDIR);
				read_file_parameters(0, reformat);
				exit(0);
			}
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
				PRINT_INFO("standard output redirected to %s.\n", filename);
			stdout=file;
			}
		}
		else {
		remove(filename);
		}
	}

	Allocate_Memory();	 

	dx = (xmax-xmin) / (Nx-1);
	dy = (ymax-ymin) / (Ny-1);
	dxy = sqrt(dx*dx+dy*dy);
	/*Change units*/
	dt *= Matosec;
	dt_eros *= Matosec;
	dt_ice *= Matosec; 
	tau *= Matosec;
	dt_record *= Matosec;
	Timefinal *= Matosec;
	Timeini *= Matosec;
	Time = Timeini;
	last_time_file_time = Timeini - 9999*dt_record;	/*very old*/
	Kerosdif *= 1e6/Matosec;
	Keroseol /= Matosec;
	Ksedim   /= Matosec;
	A_ice_rheo /= secsperyr;
	A_ice_slide /= secsperyr;
	rain *= 1e6/Matosec/1e3;
	if (hydro_model==1) Krain *= 1e6/Matosec/1e3/1e3;
	rain_per *= Matosec;
	evaporation_ct *= 1e6/Matosec/1e3;
	lost_rate *= 1e-2 * 1e-3;
	temp_sea_level += TEMP_FREEZE_WATER; /*converts from C to K*/
	switch_write_file_Blocks = true; 
	if (strlen(boundary_conds)<4)  boundary_conds[1]= boundary_conds[2]= boundary_conds[3]= boundary_conds[0];
	if (strlen(eros_bound_cond)<4) eros_bound_cond[1]=eros_bound_cond[2]=eros_bound_cond[3]=eros_bound_cond[0];

	PRINT_INFO("Nodes: %dx%d; Cell: %.2fx%.2fkm; Domain: [%.1f,%.1f]x[%.1f,%.1f]; Surface: %.0f km2", Nx, Ny, dx/1e3, dy/1e3, xmin/1e3, xmax/1e3, ymin/1e3, ymax/1e3, Nx*dx*Ny*dy/1e6);
	PRINT_INFO("Densities: asthenosphere=%f; crust=%.2f", densasthen, denscrust);
	PRINT_INFO("Boundary conditions: %s (flexure); %s (transport).", boundary_conds, eros_bound_cond);
	PRINT_INFO("Timing: from %.2f to %.2f each %.2f My\n", Timeini/Matosec, Timefinal/Matosec, dt/Matosec);

	if (verbose_level>=1) {time_t ltime; time(&ltime); fprintf(stdout, "\nTime start: %s", ctime(&ltime));}

	/*Test of incompatibilities between parameters*/
	if (switch_topoest && water_load)	{ water_load=0; 	PRINT_WARNING("Sea not possible when topoest<>0. Sea load switch turned off.") ; }
	if (densenv && water_load)			{ water_load=0; 	PRINT_WARNING("Sea not possible when densenv<>0. Sea load switch turned off.") ; }
	if (evaporation_ct && evaporation_ct<rain)	{ 			PRINT_WARNING("Evaporation should exceed the rain to produce endorheic lakes.") ; }
	if (K_ice_eros && !hydro_model)		{ K_ice_eros=0; 	PRINT_WARNING("No ice flow if hydro_model==0; K_ice_eros turned off.") ; }
	if (switch_topoest && !densinfill) 	{ 					PRINT_WARNING("Infill density has 0 value while topography has been selected to remain at initial level.") ; }
	if (plotting_mode && !switch_write_file){	 				PRINT_WARNING("switch_write_file should be turned on. Postscript may not be produced.") ; }
	if (tau<=0 && isost_model==2)  		{ isost_model=1; 	PRINT_WARNING("Negative tau will be ignored: elastic plate is assumed.");}
	if (Nx<6 || Ny<6)					{					PRINT_WARNING("Too coarse %dx%d gridding.", Nx, Ny); }

	snprintf(command, sizeof(command), "rm -f .*.tisc.tmp %s.mtrz", projectname);
	system(command);

	read_file_sea_level(); sea_level = calculate_sea_level(Time);
	read_file_insolation();
	read_file_horiz_record_time();
	read_file_Te();
	read_file_rain(precipitation_file);

	switch_initial_geom = read_file_initial_deflection(w) + read_file_initial_topo(topo) ;
	for (i=0; i<Ny; i++)  for (j=0; j<Nx; j++)  
		topo[i][j] += random_topo * ((((float) rand()) / ((float) RAND_MAX)) -.5);
	read_file_initial_rivers();
	for (i=0; i<Ny; i++)  for (j=0; j<Nx; j++)  {
		topo[i][j] += zini;
		Blocks_base[i][j] = topo[i][j];
		topo[i][j] -= w[i][j];
		h_water[i][j] = MAX_2(sea_level-topo[i][j], 0);
	}
	if (switch_initial_geom) {
		PRINT_GRID_INFO (topo, "topo_ini ", "m");
	}

	return(1);
}



int interpr_command_line_opts(int argc, char **argv) 
{
	/*Interpretates the command line options*/

	PRINT_INFO("Starting command line interpretation.");
	for (int iarg=1; iarg<argc; iarg++) {
		if (argv[iarg][0] == '-') {
			float 	value, value2;
			char 	prm[MAXLENLINE], prm2[MAXLENLINE], *ptr;
			for (int ilet=2; ilet<strlen(argv[iarg])+2; ilet++) 
				prm[ilet-2] = argv[iarg][ilet];
			for (int ilet=3; ilet < strlen(argv[iarg])+2; ilet++) 
				prm2[ilet-3] = argv[iarg][ilet];
			value  = atof(prm);
			value2 = atof(prm2);
			PRINT_DEBUG("\aArgument: %s", argv[iarg]);
			switch (argv[iarg][1]) {
				case 'B':
					strncpy(boundary_conds, prm, sizeof(boundary_conds) - 1);
					boundary_conds[sizeof(boundary_conds) - 1] = '\0';
					break;
				case 'D':
					if (run_type!=2) {
						xmin = atof(strtok(prm, "/"));
						xmax = atof(strtok(NULL, "/"));
						ymin = atof(strtok(NULL, "/"));
						ymax = atof(strtok(NULL, "/"));
					}
					else PRINT_WARNING("Impossible to change the domain when resuming a model.");
					break;
				case 'd':
					if (run_type!=2) {
						Nx = (int) (xmax-xmin)/atof(strtok(prm, "/")) + 1;
						ptr=strtok(NULL, "/");
						if (ptr != NULL) Ny = (int) (ymax-ymin)/atof(ptr) + 1;
						else Ny = Nx;
					}
					else PRINT_WARNING("Impossible to change gridding values when resuming a model.");
					break;
				case 'e':
					K_river_cap = atof(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) Kerosdif = atof(ptr);
					else Kerosdif = 0;
					break;
				case 'l':
					lake_instant_fill=1;
					break;
				case 'M':
					isost_model = value;
					break;
				case 'N':
					if (run_type!=2) {
						Nx = atoi(strtok(prm, "/"));
						ptr=strtok(NULL, "/");
						if (ptr != NULL) Ny = atoi(ptr);
						else Ny = Nx;
					}
					else PRINT_WARNING("Impossible to change grid when resuming a model.");
					break;
				case 'o':
					switch_file_out=true;
					break;
				case 'P':
					plotting_mode=(strlen(prm)>0)?(int)value:1; /*default is 1 to keep old command line syntax*/
					switch_write_file_Blocks=true;
					if (argv[iarg][2] == 'c') {
						plotting_mode = 2;
						gif_geom[0] = '\0';
						if (strlen(prm2)>0) {
							strncpy(gif_geom, prm2, sizeof(gif_geom) - 1);
							gif_geom[sizeof(gif_geom) - 1] = '\0';
						}
					} else if (argv[iarg][2] == 'p' || argv[iarg][2] == '3') {
						plotting_mode = 3;
						gif_geom[0] = '\0';
						if (strlen(prm2)>0) {
							strncpy(gif_geom, prm2, sizeof(gif_geom) - 1);
							gif_geom[sizeof(gif_geom) - 1] = '\0';
						}
					} else if (plotting_mode >= 2) {
						gif_geom[0] = '\0'; // Clear gif_geom
						if (strlen(prm2)>0) { // Check if prm2 is not empty
							strncpy(gif_geom, prm2, sizeof(gif_geom) - 1); // Use strncpy for safety
							gif_geom[sizeof(gif_geom) - 1] = '\0'; // Ensure null termination
						}
					}
					break;
				case 'p':
					rain = atof(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) Krain = atof(ptr);
					ptr=strtok(NULL, "/");
					if (ptr != NULL) evaporation_ct = atof(ptr);
					break;
				case 'q':
					ptr = strtok(prm, "=");
					match_parameter(ptr, strtok(NULL, "/"), 1, 0, prm);
					break;
				case 'R':
					random_topo = atof(strtok(prm, "/"));
					PRINT_INFO("random_topo: %f", random_topo);
					ptr=strtok(NULL, "/");
					if (ptr != NULL) {
						int seed_int;
						seed_int = atoi(ptr); 
						srand(seed_int);
						PRINT_INFO("random seed: %.1d", seed_int);
					}
					break;
				case 'r':
					switch (argv[iarg][2]) {
						case 'e':	densenv    = value2; 	break;
						case 'c':	denscrust  = value2; 	break;
						case 'i':	densinfill = value2; 	break;
						case 'm':	densmantle = value2; 	break;
						case 'a':	densasthen = value2; 	break;
					}
					break;
				case 'S':
					{
						int iblock, nblocks;
						struct BLOCK Block_aux;
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
					solver_type = argv[iarg][2];
					break;
				case 'T':
					Te_default = value;
					break;
				case 't':
					switch (argv[iarg][2]) {
						case 'i':	Timeini    = value2; 	if (run_type==2) Timeini *= Matosec;	break;
						case 'f':	Timefinal  = value2; 	if (run_type==2) Timefinal *= Matosec;	break;
						case 'd':	dt         = value2; 	if (run_type==2) dt *= Matosec;		break;
						case 'e':	dt_eros    = value2; 	if (run_type==2) dt_eros *= Matosec;	break;
						case 'v':	tau        = value2; 	if (run_type==2) tau *= Matosec; 	break;
						case 'r':	dt_record = value2; 	if (run_type==2) dt_record *= Matosec;	break;
					}
					break;
				case 'V':
					verbose_level = 2;
					if (argv[iarg][2]) verbose_level = (int)value;
					break;
				case 'v':
					{
						double density, vel_x, vel_y;
						density = atof(strtok(prm, "/"));
						vel_x = atof(strtok(NULL, "/"));
						vel_y = atof(strtok(NULL, "/"));
						for (int iu=0; iu<numBlocks; iu++) {
							if (Blocks[iu].density==-density || iu==density) {
								Blocks[iu].vel_x[0][0]=vel_x*1e3/Matosec;
								Blocks[iu].vel_y[0][0]=vel_y*1e3/Matosec;
								Blocks[iu].last_vel_time=Time-dt;/*!!*/
								Blocks[iu].last_shift_x=0;
								Blocks[iu].last_shift_y=0;
							}
						}
					}
					break;
			}
		}
	}
	return(1);
}



int Direct_mode(char *load_file_name)
{
	int 	i, j;
	FILE 	*file;
	ModelConfig cfg = {0}; 
	ModelContext ctx = {0};

	/*Solves flexure problem in direct mode: taking a single external load file and 
	writing deflection in standar output. There are no other input files neither 
	output files*/

	PRINT_INFO("Entering direct mode. xmin/xmax/ymin/ymax=%.1f/%.1f/%.1f/%.1f Nx/Ny=%d/%d\n", xmin,xmax,ymin,ymax, Nx,Ny);
	dx = (xmax-xmin) / (Nx-1) ;
	dy = (ymax-ymin) / (Ny-1) ;
	cfg.Nx = Nx; cfg.Ny = Ny; cfg.dx = dx; cfg.dy = dy; cfg.Px = Px; cfg.Py = Py; cfg.Pxy = Pxy;
	cfg.verbose_level = verbose_level;
	ctx.Time = Time; ctx.Timeini = Timeini;
	Allocate_Memory();
	if (strcmp(load_file_name, "")) {
		if ((file = fopen(load_file_name, "rt")) == NULL) {
			PRINT_INFO("Load file '%s' not found.\n", load_file_name);
			exit(0);
		}
		readinterp2D(file, h_last_unit, mode_interp, 0, xmin, xmax, ymin, ymax, Nx, Ny) ;
		fclose(file);
	}
	for (i=0; i<Ny; i++) for (j=0; j<Nx; j++) {
		D[i][j] = ET2RIG(Te_default); 
		Dq[i][j] = h_last_unit[i][j];
	}
	Elastic_Deflection(&cfg, &ctx);
	fprintf(stdout, "\n\nx[km]\t\ty[km]\t\tw[m]\t\tpressure[Pa]\n"); 
	for (i=0; i<Ny; i++) for (j=0; j<Nx; j++) 
		fprintf(stdout, "%8.1f\t%8.1f\t%8.1f\t%8.1f\n", 
			(xmin+j*dx)/1e3, (ymax-i*dy)/1e3, w[i][j], h_last_unit[i][j]);
	fprintf(stdout, "\n"); 
	return (1);
}




int tectload(ModelConfig *cfg, ModelContext *ctx)
{
	/*
	CALCULATES NEW LOAD INCREMENT FROM UNIT FILES, Returns 1 if elastic
	flexure must be done (i.e, if changes in load  occurred), 0 otherwise.
	*/
	
	PRINT_GRID_INFO (ctx->topo, "topogr.  ", "m");

	/*Moves Blocks*/
	move_Blocks(cfg, ctx);

	/*Reads external load from file*/
	while (read_file_unit(cfg, ctx));

	/*Distributes the emplacement of a unit along time*/
	gradual_Block(cfg, ctx);

	Repare_Blocks(cfg, ctx);

	return (1);
}



int Elastic_Deflection(ModelConfig *cfg, ModelContext *ctx)
{
	int 	i, j, NDi=2*cfg->Ny, NDs=2*cfg->Ny, Neqs=cfg->Nx*cfg->Ny;
	double	**A, *b;
	bool	load_changes=false;

	for (i=0; i<cfg->Ny; i++) for (j=0; j<cfg->Nx; j++) if (Dq[i][j]) load_changes = true;
	if (isost_model>0 && (load_changes || (ctx->Time==ctx->Timeini && (cfg->Px || cfg->Py || cfg->Pxy)))) {
    	    if (!Te_default) {
				/*LOCAL ISOSTASY*/
				float Krest;
				for (i=0; i<cfg->Ny; i++) for (j=0; j<cfg->Nx; j++)  {
				    GET_KREST(Krest, q, i,j)
				    Dw[i][j] = Dq[i][j] / Krest;
				}
    	    }
    	    else {
    	    	/*REGIONAL ISOSTASY*/
    	    	switch (solver_type) {
    		case 'l':
    		    /*Requires 4*Nx*Ny*Ny cells*/
    			b = (double *) calloc (Neqs, sizeof(double));
				A = alloc_matrix_dbl (Neqs, NDi+1+NDs);
    			defineLESalmostdiagonalmatrix(cfg, A, b, q, Dq, w, 0);
    			solveLESalmostdiagonal(cfg, A, b, Dw);
    			free_matrix_dbl (A, Neqs);
    			free(b);
    			break;
    		case 'm':
			{
    		    int nonzeroes=13*cfg->Nx*cfg->Ny, ESP, PATH, FLAG, NSP, *R, *C, *IC, *IA, *JA, *ISP;
    		    float *B, *Z, *mathlib_matrix, *RSP;
    		    NSP = 40*(6*Neqs+2+nonzeroes);
    		    PATH = 1;
		    B = (float *) calloc (Neqs, sizeof(float));
		    Z = (float *) calloc (Neqs, sizeof(float));
    		    mathlib_matrix = (float *) calloc (nonzeroes, sizeof(float));
    		    R =     (int *) calloc (Neqs, sizeof(int));
    		    C =     (int *) calloc (Neqs, sizeof(int));
    		    IC =    (int *) calloc (Neqs, sizeof(int));
    		    IA =    (int *) calloc (Neqs+1, sizeof(int));
    		    JA =    (int *) calloc (nonzeroes, sizeof(int));
    		    ISP =   (int *) calloc (1*NSP, sizeof(int));
    		    RSP =   (float *) calloc (NSP, sizeof(float));
    		    defineLESmatrix_for_mathlib(cfg, mathlib_matrix, IA, JA, B, Dq, w, &nonzeroes, 0);
    		    for (i=0; i<Neqs; i++)  R[i]=IC[i]=C[i]=i+1;
#ifdef MATHLIB_SOLVER
    		    printf("\nNeqs=%d NSP=%d nonzeroes=%d", Neqs, NSP, nonzeroes);
    		    cdrv(&Neqs, R, C, IC, IA, JA, mathlib_matrix, B, Z, &NSP, ISP, RSP, &ESP, &PATH, &FLAG);
    		    if (FLAG!=0) {PRINT_ERROR("\aTDRV exit #%d. Memory excess=%d\n", FLAG, ESP); }
    		    fprintf(stdout, "\tMemory excess=%d\n", ESP);
#endif
    		    for (i=0; i<cfg->Ny; i++) for (j=0; j<cfg->Nx; j++)	Dw[i][j] = Z[j*cfg->Ny+i];

    		    free(Z); free(B); free(mathlib_matrix);
    		    free(R); free(C);free(IC); free(IA); free(JA); free(ISP); free(RSP);
    		    break;
			}
    	      }
	    }

    	    for (i=0; i<cfg->Ny; i++) for (j=0; j<cfg->Nx; j++)  w[i][j] += Dw[i][j];
    	    if (switch_topoest) {
   		    /*Defines the thickness of last infill Block*/
    		    for (i=0; i<cfg->Ny; i++)  for (j=0; j<cfg->Nx; j++)  
		    	Blocks[i_first_Block_load-1].thick[i][j] +=  MAX_2(Dw[i][j], 0) ;
    	    }


    	    //calculate_topo(topo);

    	    /*Statistics on load*/
    	    {
    		    int i, j, imax, jmax, imin, jmin;
    		    float maxloadinc=-1e25, minloadinc=1e25, total_load=0;
    		    for (i=0;i<cfg->Ny;i++) for (j=0;j<cfg->Nx;j++) {
    			q[i][j] += Dq[i][j];
    			total_load += Dq[i][j]*cfg->dx*cfg->dy;
    			if (maxloadinc<Dq[i][j]) {maxloadinc=Dq[i][j]; imax=i; jmax=j;}
    			if (minloadinc>Dq[i][j]) {minloadinc=Dq[i][j]; imin=i; jmin=j;}
    		    }
    		    if (load_changes) 
		    	PRINT_SUMLINE("load  now:  max = %+8.2e N/m2  min = %+8.2e N/m2   Total: %+8.2e N\tMax at %.0f,%.0f km", 
		    		maxloadinc, minloadinc, total_load, (cfg->xmin+cfg->dx*jmax)/1e3, (cfg->ymax-cfg->dy*imax)/1e3); 
    	    }
    	    /*Statistics on deflection*/
    	    {
    		    float total_load=0, total_restitutive_force=0, Krest;
    		    for (i=0; i<cfg->Ny; i++) for (j=0; j<cfg->Nx; j++)	{
    			    total_load += Dq[i][j];
    			    GET_KREST(Krest, q, i,j);
    			    total_restitutive_force += (Krest*w[i][j]);
    		    }
    		    if (cfg->verbose_level>=1) {
    			    PRINT_SUMLINE("load: %+8.2e N   restit_force: %+8.2e N", total_restitutive_force*cfg->dx*cfg->dy, total_load*cfg->dx*cfg->dy);
			    PRINT_GRID_INFO (w, "deflect. ", "m");
    		    }
    	    }
	}

	/*Resets deflection and load grids*/
	for (i=0; i<cfg->Ny; i++)  for (j=0; j<cfg->Nx; j++)  Dq[i][j]=/*Dw[i][j]=*/0;

	return(1);
}





int move_Blocks(ModelConfig *cfg, ModelContext *ctx)
{
	int	*nshift_x, *nshift_y;
	float	**new_thick;
#ifdef THIN_SHEET
	char 	tmpTSBCfilename[84]="";
#endif

	/*
	  Moves the Blocks and calculates the isostatic load and thickness change at each cell.
	  Deforms the sediment above.
	*/

	new_thick = alloc_matrix(cfg->Ny, cfg->Nx);
	nshift_x = calloc(ctx->numBlocks, sizeof(int));
	nshift_y = calloc(ctx->numBlocks, sizeof(int));

	for (int iu=0; iu<ctx->numBlocks; iu++) {
		if (Blocks[iu].density == cfg->denssedim) {
/*!!*/
//Blocks[iu].vel_y[0][0] = 2.5e3/Matosec;
//Blocks[iu].last_vel_time = Blocks[iu].age;
//Blocks[iu].last_shift_x  = 0;
//Blocks[iu].last_shift_y  = 0;
//Blocks[iu].time_stop     = 1e19;
		/*DEFORM SEDIMENT Blocks*/
		for (int i=0; i<cfg->Ny; i++) for (int j=0; j<cfg->Nx; j++) 
			new_thick[i][j] = Blocks[iu].thick[i][j];
		for (int i=0; i<cfg->Ny; i++) for (int j=0; j<cfg->Nx; j++){
			float sedthick;
			sedthick=Blocks[iu].thick[i][j];
			/*Find the uppermost moving Block below this point in this sedim. Block*/
			for (int ju=iu-1; ju>=0; ju--) {
			/*Calculate the thickness of sediments between the top of this sed. Block and the moving Block*/
			if (Blocks[ju].density == cfg->denssedim) {
				sedthick += Blocks[ju].thick[i][j];
			}
			else {
   			  /*Amount of cells to propagate the deformation: ~20 deg assumed.*/
			  int nprop_x = SIGN(nshift_x[ju]) * (int) ceil(sedthick*3/cfg->dx);
			  int j_unprop = j-nprop_x;
			  int nprop_y = SIGN(nshift_y[ju]) * (int) ceil(sedthick*3/cfg->dy);
			  int i_unprop = i+nprop_y;
			  DOMAIN_LIMIT(i_unprop, j_unprop);
			  if (Blocks[ju].thick[i_unprop][j_unprop]>.1) {
				if (!nshift_y[ju] && !nshift_x[ju]) 
					break;
				else {
 					int i_shift, j_shift, i_unshift, j_unshift;
					i_shift =   i-nshift_y[ju];	j_shift =   j+nshift_x[ju];
					i_unshift = i+nshift_y[ju];	j_unshift = j-nshift_x[ju];
					/*If block ju is moving below [i][j] then shift seds.*/
					if (deform_sed && IN_DOMAIN(i_shift,j_shift)) 
					new_thick[i_shift][j_shift] += Blocks[iu].thick[i][j];
					if (deform_sed && IN_DOMAIN(i,j) && IN_DOMAIN(i_unshift,j_unshift)) 
					new_thick[i][j] -= Blocks[iu].thick[i][j];
					break;
				}
			  }
			}
			}
		}
		for (int i=0; i<cfg->Ny; i++) for (int j=0; j<cfg->Nx; j++) {
			Dq[i][j] += g * (new_thick[i][j] - Blocks[iu].thick[i][j]) * Blocks[iu].density;
			if (new_thick[i][j]<-1) PRINT_ERROR("negative sediment thickness: %.1f m", new_thick[i][j]);
			Blocks[iu].thick[i][j] = new_thick[i][j];
		}
		}
		else //!!
		{
		  if (Blocks[iu].type == 'V' && ctx->Time < Blocks[iu].time_stop) {
#ifdef THIN_SHEET
			int i,j, n, m, nn, nincogn, nbanda, 
			thicken_BC=1, /*0 means No temporal thickening variations on the boundaries, thickness=thickness_old; IBC_thicken!=0 -> No lateral variations of the thickening on the boundaries, d(thickness)/dx=d(thickness)/dy=0*/
			nitermax=0;  /*0 means viscosity is independent from strain rate => no iterations; normal=50*/
		double elapsed_time, dt_aux, Lx, Ly, tallmax=.015, alfa=.5, g=9.81, 
			vissup=1e25, visinf=1e21, /*1e25 & 1e22*/
			roalfa=3.5e-5, Zcompens=130e3;
		double *viscTer, *viscosity, 
			*average_pressure, 
			*vel_x_array, *vel_y_array, *vert_strain_rate, 
			*layer_thickness;
		/*THIN SHEET Block DEFORMATION*/
		elapsed_time = ctx->Time-Blocks[iu].last_vel_time;
		dt_aux = ctx->dt;
		Lx = cfg->xmax-cfg->xmin; Ly = cfg->ymax-cfg->ymin;
		n=cfg->Nx-1; m=cfg->Ny-1; nn=cfg->Nx*cfg->Ny; nincogn=2*nn; nbanda=4*(cfg->Nx)+7;
		viscTer        =	alloc_array_dbl(cfg->Nx*cfg->Ny);
		viscosity      =	alloc_array_dbl(cfg->Nx*cfg->Ny);
		average_pressure=	alloc_array_dbl(cfg->Nx*cfg->Ny);
		vel_x_array    =	alloc_array_dbl(cfg->Nx*cfg->Ny);
		vel_y_array    =	alloc_array_dbl(cfg->Nx*cfg->Ny);
		vert_strain_rate =	alloc_array_dbl(cfg->Nx*cfg->Ny);
		layer_thickness =	alloc_array_dbl(cfg->Nx*cfg->Ny);
		sprintf(tmpTSBCfilename, "%s"".TSBC.tmp", projectname);
		reformat_file_thin_sheet_BC(cfg, ctx, tmpTSBCfilename);
		PRINT_INFO("Deforming thin_sheet Block %d", iu);
		for (i=0; i<cfg->Ny; i++) for (j=0; j<cfg->Nx; j++)  {
			int k=(cfg->Ny-1-i)*cfg->Nx+j, l, m;
			viscTer[k] = Blocks[iu].viscTer[i][j];  /*.1e7*/  /* viscTer[Pa] = 1/2 * strength[Pa*m] / layer.thickness[m] */
			/*AVERAGE PRESSURE*/
			/*Water term*/
			average_pressure[k] = g*cfg->denswater*h_water[i][j];
			/*Blocks terms from the top to the thin_sheet*/
			for (l=ctx->numBlocks-1; l>iu; l--) {
				average_pressure[k] += g*Blocks[l].density*Blocks[l].thick[i][j];
			}
			/*Thin sheet term*/
			average_pressure[k] += g/2*Blocks[iu].density*Blocks[iu].thick[i][j];
			/*Asthenosphere restitutive push*/
			average_pressure[k] += (densasthen-cfg->densenv)*g*w[i][j];
			/*change sign to <0, meaning compresion*/
			average_pressure[k] *= -1; 

			vert_strain_rate[k] = 0;
			layer_thickness[k] = Blocks[iu].thick[i][j];
			if (elapsed_time > ctx->dt) {
				vel_x_array[k] = Blocks[iu].vel_x[i][j];
				vel_y_array[k] = Blocks[iu].vel_y[i][j];
				viscosity[k]   = Blocks[iu].visc[i][j];
			}
			else {
				vel_x_array[k] = 0;
				vel_y_array[k] = 0;
				viscosity[k] =  (double) Blocks[iu].viscTer[i][j] / 3.17e-17; /*viscosity[Pa*s] = viscTer/str.rate  3.17e-16 s^-1==1%/My*/
			}
		}
		{
			float apmin=1e19, apmax=-1e19;
			float vimin=1e19, vimax=-1e19;
			for (i=0; i<cfg->Ny; i++) for (j=0; j<cfg->Nx; j++)  {
				int k=(cfg->Ny-1-i)*cfg->Nx+j;
				apmin = MIN_2(apmin, average_pressure[k]);
				apmax = MAX_2(apmax, average_pressure[k]);
 				vimin = MIN_2(vimin, viscosity[k]);
				vimax = MAX_2(vimax, viscosity[k]);
 				/*printf("\t%d %d av_pres=%.2e", i,j, average_pressure[k]);*/
				/*if (i==1 && j==1)   printf("\n%d,%d: ap=%.4e", i,j, average_pressure[k]);
				if (i==19 && j==19) printf("\n%d,%d: ap=%.4e", i,j, average_pressure[k]);*/
			}
			PRINT_INFO("Viscosity  min = %.4e ;  max = %.4e ;  elapsed_time= %.2e s", vimin, vimax, elapsed_time);
			PRINT_INFO("Av. press. min = %.4e ;  max = %.4e", apmin, apmax);
		}
			velocity_field_(&elapsed_time, &dt_aux, &Lx, &Ly, 
			&n, &m, &nn, &vissup, &visinf, viscTer, viscosity, 
			&nincogn, &tallmax, &alfa, &nitermax, average_pressure, 
			vel_x_array, vel_y_array,
			tmpTSBCfilename, &nbanda
		);
		vertical_strain_rate_(&Lx, &Ly, &n, &m, &nn, 
			vel_x_array, vel_y_array, vert_strain_rate
		);
		thicken_(&dt_aux, &n, &m, &nn, &Lx, &Ly, 
			vel_x_array, vel_y_array, vert_strain_rate, layer_thickness, 
			&thicken_BC
		);
		for (i=0; i<cfg->Ny; i++)  for (j=0; j<cfg->Nx; j++)  {
			int k=(cfg->Ny-1-i)*cfg->Nx+j;
			Dq[i][j] += g * (layer_thickness[k] - Blocks[iu].thick[i][j]) * Blocks[iu].density;
			Blocks[iu].thick[i][j] = layer_thickness[k];
			Blocks[iu].vel_x[i][j] = vel_x_array[k];
			Blocks[iu].vel_y[i][j] = vel_y_array[k];
			Blocks[iu].visc[i][j]  = viscosity[k];
		}
		free(viscTer);
		free(viscosity);
		free(average_pressure);
		free(vel_x_array);
		free(vel_y_array);
		free(vert_strain_rate);
		free(layer_thickness);
		remove(tmpTSBCfilename); 
#endif
		  }
		  else {
			int i, j, i_unshifted, j_unshifted;
			float theor_shift_x, theor_shift_y;
			/*MOVE BLOCKS*/
			theor_shift_x = Blocks[iu].vel_x[0][0] * (ctx->Time-Blocks[iu].last_vel_time);
			theor_shift_y = Blocks[iu].vel_y[0][0] * (ctx->Time-Blocks[iu].last_vel_time);
			nshift_x[iu] = floor((theor_shift_x - Blocks[iu].last_shift_x) /cfg->dx +.5);
			nshift_y[iu] = floor((theor_shift_y - Blocks[iu].last_shift_y) /cfg->dy +.5);
			if (ctx->Time > Blocks[iu].time_stop + .1*ctx->dt) {nshift_x[iu]=0; nshift_y[iu]=0;}
			Blocks[iu].shift_x += nshift_x[iu]*cfg->dx;
			Blocks[iu].shift_y += nshift_y[iu]*cfg->dy;
			Blocks[iu].last_shift_x += nshift_x[iu]*cfg->dx;
			Blocks[iu].last_shift_y += nshift_y[iu]*cfg->dy;
			for (i=0; i<cfg->Ny; i++) for (j=0; j<cfg->Nx; j++){
				i_unshifted = i+nshift_y[iu];	j_unshifted = j-nshift_x[iu];	
				DOMAIN_LIMIT(i_unshifted, j_unshifted);	/*[i][j], unshifted*/
				new_thick[i][j] = Blocks[iu].thick[i_unshifted][j_unshifted];
			}
			for (i=0; i<cfg->Ny; i++) for (j=0; j<cfg->Nx; j++) {
				Dq[i][j] += g * (new_thick[i][j] - Blocks[iu].thick[i][j]) * Blocks[iu].density;
				Blocks[iu].thick[i][j] = new_thick[i][j];
			}
		  }
		}
	}

	free_matrix(new_thick, cfg->Ny);
	free(nshift_x);
	free(nshift_y);
	return(1);
}




int read_file_unit(ModelConfig *cfg, ModelContext *ctx)
{
	/*
	  READS UNIT FILE NAMED 'projectname[NUM].UNIT' WHERE 'NUM' IS 1 FOR THE 
	  FIRST UNIT, 2 FOR THE SECOND, ETC. Interpolates this input. Creates
	  new unit to store its properties and cuts sediment units when file
	  contains fault depth rather than a thickness itself.
	*/

	int 	nparams=0;
	float	time_stop=9999/*My*/, time_unit, 
		erodibility_aux=NO_DATA, fill_up_to=NO_DATA, 
		vel_x=0, vel_y=0, density=NO_DATA;
	bool 	insert, cut_Blocks, cut_all, top, fault, switch_move, 
		thin_sheet, ride, hidden, z_absol;
	FILE 	*file;
	char 	filename[MAXLENFILE];

	/*Read the next unit age*/
	sprintf(filename, "%s%d.UNIT", projectname, nloads+1);
	if ((file = fopen(filename, "rt")) == NULL) {
		PRINT_INFO("Cannot read unit file '%s'.", filename);
		return (0);
	}
	time_unit = Timeini/Matosec;
	{
		int nlines=0, nread, show, replace=0;
		char str1[MAXLENLINE], str2[MAXLENLINE], line[MAXLENLINE+200], *lineptr;
		show=(cfg->verbose_level>=4)? 1 : 0;
		rewind(file);
		while ((lineptr=fgets(line, MAXLENLINE+200-1, file)) != NULL && nlines<NMAXHEADERLINES) {
			nlines++; nread=sscanf(lineptr, "%s %s", str1, str2);
			if (nread == 2) {
				Match_Param_Replace_flt ( "time",  	time_unit, 0 )
				/*Old versions:*/
				Match_Param_Replace_flt ( "time_load",	time_unit, 1 )
			}
		}
	}
	time_unit *= Matosec;
	/*Return if it isn't time yet to read the new unit file*/
	if (time_unit>ctx->Time+.1*ctx->dt || time_unit<ctx->Timeini) return(0);

	PRINT_INFO("Reading '%s'", filename);
	switch_move = fault = switch_gradual = 
		insert = hidden = cut_Blocks = cut_all = 
		thin_sheet = top = ride = z_absol = false;
	i_Block_insert = ctx->numBlocks;

	/*READS AND INTERPOLATES UNIT/LOAD FILE*/
	{
		int nlines=0, nread, show, replace=0;
		char str1[MAXLENLINE], str2[MAXLENLINE], line[MAXLENLINE+200], *lineptr;
		show=(cfg->verbose_level>=3)? 1 : 0;
		rewind(file); 
		while ((lineptr=fgets(line, MAXLENLINE+200-1, file)) != NULL && nlines<NMAXHEADERLINES) {
				nlines++; nread=sscanf(lineptr, "%s %s", str1, str2);
				if (nread == 2) {
				Match_Param_Replace_flt ( "vel_x",		vel_x,   	0 )
				Match_Param_Replace_flt ( "vel_y",		vel_y,   	0 )
				Match_Param_Replace_flt ( "time_stop",  	time_stop,   	0 )
				Match_Param_Replace_flt ( "density", 	 	density,   	0 )
				Match_Param_Replace_flt ( "erodibility",	erodibility_aux,   	0 )
				Match_Param_Replace_int ( "gradual",		switch_gradual,   	0 )
				Match_Param_Replace_int ( "hidden",		hidden,   	0 )
				Match_Param_Replace_int ( "ride",		ride,   	0 )
				Match_Param_Replace_int ( "insert",		insert,   	0 )
				Match_Param_Replace_int ( "top",		top,   	0 )
				Match_Param_Replace_int ( "move",  		switch_move,   	0 )
				Match_Param_Replace_int ( "fault",	  	fault,   	0 )
				Match_Param_Replace_int ( "z_absol",	  	z_absol,   	0 )
				Match_Param_Replace_int ( "cut_Blocks",  	cut_Blocks,   	0 )
				Match_Param_Replace_int ( "cut_all",  		cut_all,   	0 )
				Match_Param_Replace_int ( "thin_sheet",		thin_sheet,   	0 )
				Match_Param_Replace_int ( "topoest",		switch_topoest,   	0 )
				Match_Param_Replace_int ( "densenv",		cfg->densenv,   	0 )
				Match_Param_Replace_flt ( "fill_up_to", 	fill_up_to,   	0 )
				/*Old versions:*/
				Match_Param_Replace_int ( "fault_load", 	fault,   	1 )
				Match_Param_Replace_int ( "interp_load",	switch_gradual ,   	1)
				Match_Param_Replace_int ( "hidden_load",	hidden,   	1 )
				Match_Param_Replace_flt ( "dens_load",  	density,   	1 )
				Match_Param_Replace_int ( "insert_load",	insert,   	1 )
				Match_Param_Replace_int ( "top_load",		top,   	1 )
				Match_Param_Replace_int ( "move_load",  	switch_move,   	1 )
				Match_Param_Replace_int ( "cut_loads",  	cut_Blocks,   	1 )
				Match_Param_Replace_flt ( "erodability",	erodibility_aux,   	1 )
				Match_Param_Replace_flt ( "l_fluv_eros",	erodibility_aux,   	1 )
				}
				if (strcmp(str1, "thickness_distribution")==0) break;
		}
		rewind(file); 
	}
	if (fill_up_to == NO_DATA) 
		readinterp2D(file, h_last_unit, mode_interp, 0, cfg->xmin, cfg->xmax, cfg->ymin, cfg->ymax, cfg->Nx, cfg->Ny);
	else {
		for (int i=0; i<cfg->Ny; i++) for (int j=0; j<cfg->Nx; j++) h_last_unit[i][j] = MAX_2(0, fill_up_to-ctx->topo[i][j]);
	}
	fclose(file);

	nloads++;

	vel_x *= 1e3/Matosec;
	vel_y *= 1e3/Matosec;
	time_stop *= Matosec;

	/*ACT ACCORDING TO THE SIGNALS*/
	if (thin_sheet) {
		switch_move = true; 
	}
	if (fault) {
		switch_move = true;
	}

	/*Check incompatibilities between unit file signals*/
	if (switch_gradual && switch_move) {
		PRINT_WARNING("Gradual & moving Blocks are not compatible. This one won't be gradual.");
		switch_gradual = false;
	}
	if (switch_gradual && hidden) {
		PRINT_WARNING("Gradual & hidden Blocks are not compatible. This one won't be hidden.");
		hidden = false;
	}

	/*Creates a Block of infill if switch_topoest; it will be filled later during Deflection*/
	if (switch_topoest) {
		insert_new_Block(i_first_Block_load);
		ctx->numBlocks = numBlocks;
		Blocks[i_first_Block_load].type = 'I'; 	/*stands for Infill*/
		Blocks[i_first_Block_load].density = densinfill;
		if (densinfill<2550) Blocks[i_first_Block_load].erodibility = erodibility_sed;
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

	if (fault && !top) i_Block_insert = 0;

	if (fault) {
		int numBlocks0=ctx->numBlocks;
		/*CUT BlockS*/
		/*Make deep copies of all Blocks*/ // Deep copy for faulting
		PRINT_DEBUG("Cutting Blocks: numBlocks= %d", ctx->numBlocks);
		for (int k=0; k<numBlocks0; k++) {
			// Create a new block at the end
			insert_new_Block(ctx->numBlocks);
			ctx->numBlocks = numBlocks;
			struct BLOCK *new_block = &Blocks[ctx->numBlocks-1];
			struct BLOCK *src_block = &Blocks[k];

			// Save pointers allocated by insert_new_Block
			float **new_thick = new_block->thick;
			float **new_vel_x = new_block->vel_x;
			float **new_vel_y = new_block->vel_y;
			float **new_visc = new_block->visc;
			float **new_viscTer = new_block->viscTer;

			// Handle 'V' type blocks
			if (src_block->type == 'V') {
				free_matrix(new_vel_x, 1);
				free_matrix(new_vel_y, 1);
				free_matrix(new_visc, 1);
				free_matrix(new_viscTer, 1);
				new_vel_x = alloc_matrix(cfg->Ny, cfg->Nx);
				new_vel_y = alloc_matrix(cfg->Ny, cfg->Nx);
				new_visc = alloc_matrix(cfg->Ny, cfg->Nx);
				new_viscTer = alloc_matrix(cfg->Ny, cfg->Nx);
			}

			// Allocate 'S' type arrays if needed
			float **new_detr_ratio = NULL;
			float **new_detr_grsize = NULL;
			float **new_thickgypsum = NULL;
			float **new_thickhalite = NULL;
			if (src_block->type == 'S') {
				new_detr_ratio = alloc_matrix(cfg->Ny, cfg->Nx);
				new_detr_grsize = alloc_matrix(cfg->Ny, cfg->Nx);
				new_thickgypsum = alloc_matrix(cfg->Ny, cfg->Nx);
				new_thickhalite = alloc_matrix(cfg->Ny, cfg->Nx);
			}

			// Copy scalar members
			*new_block = *src_block; // Shallow copy of scalars

			// Restore pointers
			new_block->thick = new_thick;
			new_block->vel_x = new_vel_x;
			new_block->vel_y = new_vel_y;
			new_block->visc = new_visc;
			new_block->viscTer = new_viscTer;
			if (src_block->type == 'S') {
				new_block->detr_ratio = new_detr_ratio;
				new_block->detr_grsize = new_detr_grsize;
				new_block->thickgypsum = new_thickgypsum;
				new_block->thickhalite = new_thickhalite;
			}

			// Deep copy contents of 2D arrays
			memcpy(new_block->thick[0], src_block->thick[0], cfg->Nx * cfg->Ny * sizeof(float));
			if (src_block->type == 'S') {
				memcpy(new_block->detr_ratio[0], src_block->detr_ratio[0], cfg->Nx * cfg->Ny * sizeof(float));
				memcpy(new_block->detr_grsize[0], src_block->detr_grsize[0], cfg->Nx * cfg->Ny * sizeof(float));
				memcpy(new_block->thickgypsum[0], src_block->thickgypsum[0], cfg->Nx * cfg->Ny * sizeof(float));
				memcpy(new_block->thickhalite[0], src_block->thickhalite[0], cfg->Nx * cfg->Ny * sizeof(float));
			}
			if (src_block->type == 'V') {
				memcpy(new_block->vel_x[0], src_block->vel_x[0], cfg->Nx * cfg->Ny * sizeof(float));
				memcpy(new_block->vel_y[0], src_block->vel_y[0], cfg->Nx * cfg->Ny * sizeof(float));
				memcpy(new_block->visc[0], src_block->visc[0], cfg->Nx * cfg->Ny * sizeof(float));
				memcpy(new_block->viscTer[0], src_block->viscTer[0], cfg->Nx * cfg->Ny * sizeof(float));
			} else {
				new_block->vel_x[0][0] = src_block->vel_x[0][0];
				new_block->vel_y[0][0] = src_block->vel_y[0][0];
				new_block->visc[0][0] = src_block->visc[0][0];
				new_block->viscTer[0][0] = src_block->viscTer[0][0];
			}

			new_block->vel_x[0][0] = vel_x; // Override velocity
			new_block->vel_y[0][0] = vel_y; // Override velocity
			new_block->last_vel_time = ctx->Time;
			new_block->last_shift_x = 0;
			new_block->last_shift_y = 0;
			new_block->time_stop = time_stop;
			if (density         != NO_DATA && Blocks[ctx->numBlocks-1].type != 'S') Blocks[ctx->numBlocks-1].density = density;
			if (erodibility_aux != NO_DATA && Blocks[ctx->numBlocks-1].type != 'S') Blocks[ctx->numBlocks-1].erodibility = erodibility_aux;
		}
		PRINT_DEBUG("Updating Blocks_base: numBlocks= %d", ctx->numBlocks);
		/*Modify Blocks_base and cut above the fault*/
		for (int i=0; i<cfg->Ny; i++) for (int j=0; j<cfg->Nx; j++) {
			float z_fault=-h_last_unit[i][j], base_of_Block=Blocks_base[i][j];
			if (z_absol) base_of_Block -= w[i][j];
			h_last_unit[i][j] = MAX_2(0, Blocks_base[i][j]-z_fault);	/*Block thickness below fault (to create the new Block, see below)*/
			Blocks_base[i][j]  = MIN_2(Blocks_base[i][j], z_fault);	/*new base of Blocks*/
			if (cut_Blocks) {
				for (int k=0; k<numBlocks0; k++) {
					float top_of_Block=base_of_Block+Blocks[k].thick[i][j];
					if (Blocks[k].density == cfg->denssedim && !cut_all) {
						break;
					}
					if (z_fault <= base_of_Block) {
						Blocks[k+numBlocks0].thick[i][j] += Blocks[k].thick[i][j];
						Blocks[k].thick[i][j]            = 0;
					}
					else {
						Blocks[k+numBlocks0].thick[i][j] += MAX_2 (0, top_of_Block-z_fault);
						Blocks[k].thick[i][j]           -= MAX_2 (0, top_of_Block-z_fault);
					}
					base_of_Block = top_of_Block;
				}
			}
		}
	}
	if (density        ==NO_DATA) density         = denscrust;
	if (erodibility_aux==NO_DATA) erodibility_aux = erodibility;

	PRINT_DEBUG("Creating Block for this file: i_Block_insert= %d", i_Block_insert);
	/*Create a new Block for the thickness in this file*/
	insert_new_Block(i_Block_insert);
	ctx->numBlocks = numBlocks;

	/*Add the thickness in file to the new Block; Shrink the Blocks and basement if the thickness is negative*/
	if (!switch_gradual && !hidden) {
		for (int i=0; i<cfg->Ny; i++)  for (int j=0; j<cfg->Nx; j++) {
			if (h_last_unit[i][j]>=0) {
				Blocks[i_Block_insert].thick[i][j] = h_last_unit[i][j];
			}
			else {
				float 	h_unit_aux, h_unit_aux2;
				int	k;
				/*Excavating rock for negative thickness in file, starting from base of water*/
				h_unit_aux = fabs((double) h_last_unit[i][j]);
				for (k=i_Block_insert-1; h_unit_aux>0 && k>=0; k--) {
					h_unit_aux2 = MIN_2(Blocks[k].thick[i][j], h_unit_aux);
					h_unit_aux -= h_unit_aux2;
					Blocks[k].thick[i][j] -= h_unit_aux2;
				}
				/*k is the deepest eroded Block*/
				if (k==-1) {
					Blocks_base[i][j] -= h_unit_aux;
				}
			}
		}
	}
	if (hidden) Blocks[i_Block_insert].type = 'H';
	if (switch_gradual) Blocks[ctx->numBlocks-1].type = 'G';
	if (thin_sheet) {
		float default_viscTerm = .5e7; /*.1e7*/
		char filename[MAXLENFILE];
		Blocks[i_Block_insert].type    = 'V';
		free_matrix(Blocks[i_Block_insert].vel_x, 1);
		free_matrix(Blocks[i_Block_insert].vel_y, 1);
		free_matrix(Blocks[i_Block_insert].visc, 1);
		free_matrix(Blocks[i_Block_insert].viscTer, 1);
		Blocks[i_Block_insert].vel_x   = alloc_matrix(cfg->Ny, cfg->Nx);
		Blocks[i_Block_insert].vel_y   = alloc_matrix(cfg->Ny, cfg->Nx);
		Blocks[i_Block_insert].visc    = alloc_matrix(cfg->Ny, cfg->Nx);
		Blocks[i_Block_insert].viscTer = alloc_matrix(cfg->Ny, cfg->Nx);
		sprintf(filename, "%s%d.VISC", projectname, nloads);
		if ((file = fopen(filename, "rt")) == NULL) {
			PRINT_WARNING("Cannot read thermal viscosity file '%s'.", filename);
			for (int i=0; i<cfg->Ny; i++)  for (int j=0; j<cfg->Nx; j++) Blocks[i_Block_insert].viscTer[i][j] = default_viscTerm;
		}
		else {
			PRINT_INFO("Reading viscosity for unit %d from file '%s'.", nloads, filename);
			readinterp2D(file, Blocks[i_Block_insert].viscTer, mode_interp, default_viscTerm, cfg->xmin, cfg->xmax, cfg->ymin, cfg->ymax, cfg->Nx, cfg->Ny);
		}
	}
	Blocks[i_Block_insert].density = density;
	Blocks[i_Block_insert].erodibility = erodibility_aux;
	Blocks[i_Block_insert].vel_x[0][0] = vel_x;
	Blocks[i_Block_insert].vel_y[0][0] = vel_y;
	Blocks[i_Block_insert].time_stop = time_stop;

	if (ride) {
		for (int i_Block=i_Block_insert+1; i_Block<ctx->numBlocks; i_Block++) {
			if (Blocks[i_Block].type == 'V' || Blocks[i_Block_insert].type == 'V') {
				PRINT_WARNING("Cannot apply ride to thin_sheet blocks or from thin_sheet blocks.");
			} else {
				Blocks[i_Block].vel_x[0][0]   = Blocks[i_Block_insert].vel_x[0][0]; 
				Blocks[i_Block].vel_y[0][0]   = Blocks[i_Block_insert].vel_y[0][0]; 
			}
			Blocks[i_Block].last_shift_x  = 0; 
			Blocks[i_Block].last_shift_y  = 0; 
			Blocks[i_Block].last_vel_time = ctx->Time; 
			Blocks[i_Block].time_stop     = Blocks[i_Block_insert].time_stop; 
		}
	}

	/*Don't Repare_Blocks() in case of: 
		Gradual load, because then h_last_unit[] will be empty until tectload()
		Topoest load, because the infill Block will be filled upon deflection.
	*/
	if (!switch_gradual && !switch_topoest) Repare_Blocks(cfg, ctx);

	/*Increment the isostatic load for this time interval*/
	if (!switch_gradual && !fault) 
		for (int i=0; i<cfg->Ny; i++) for (int j=0; j<cfg->Nx; j++) Dq[i][j] += (density-cfg->densenv)*g*h_last_unit[i][j];

	PRINT_INFO("Unit read from '%s'. ", filename);
	PRINT_DEBUG("%d params; dens=%.0f kg/m3; erodibility=%.1e; ", nparams, density, erodibility_aux);
	if (!fault) {
		float load=0;
		for (int i=0; i<cfg->Ny; i++) for (int j=0; j<cfg->Nx; j++) load += (density-cfg->densenv)*h_last_unit[i][j];
		fprintf(stdout, "%.2e kg. ", load*cfg->dx*cfg->dy);
	}
	if (switch_gradual) PRINT_INFO("Will be gradually loaded till %.2fMy.", time_stop/Matosec);
	if (switch_move) PRINT_INFO("Vel= %.2fE,%.2fN km/My till T=%.2f My", vel_x*Matosec/1000, vel_y*Matosec/1000, time_stop/Matosec);

	fprintf(stdout, " unit") ;
	if (fault)  		fprintf(stdout, "F") ;
	if (thin_sheet)  	fprintf(stdout, "V") ;
	if (hidden)			fprintf(stdout, "H") ;
	if (insert) 		fprintf(stdout, "I") ;
	if (top) 			fprintf(stdout, "P") ;
	if (switch_gradual)	fprintf(stdout, "G") ;
	if (switch_move) 	fprintf(stdout, "M") ;
	if (ride) 			fprintf(stdout, "R") ;
	if (switch_topoest)	fprintf(stdout, "T") ;

	return(1);
}




int syntax()
{
	/*
		Displays the hardcoded command line syntax of the program
	*/
	fprintf(stderr, "\nSyntax:\n");
	fprintf(stderr, "  tisc  project  -B<bound_type> -D[xmin/xmax/ymin/ymax] -d<dx>[/<dy>]\n");
	fprintf(stderr, "        -e<Kriv>[/<Kdif>] -f[2] -F[<file>] -h[i|u|p|c] -l -M<lith_type>\n");
	fprintf(stderr, "        -N<Nx>[/<Ny>] -o -P[0|1|2|3][geom] -p<rain>[/<Krain>][/<evap>] -Q<file>\n");
	fprintf(stderr, "        -q<param>=<value> -R<random_topo>[/<seed>] -r<e|c|i|m|a><dens>\n");
	fprintf(stderr, "        -s<solver> -T<eet> -t<i|f|d|e|v|r><time> -V[<level>] \n");
	fprintf(stderr, "        -v<num>/<vx>/<vy>\n\n");
	fprintf(stderr, "  Options:\n");
	fprintf(stderr, "    'project'\tRoot name for all files (e.g. 'test' looks for 'test.PRM')\n");
	fprintf(stderr, "    -F\t\tResume a previous run from a .all file\n");
	fprintf(stderr, "    -P\t\tProduce graphic output (-Pp for Python, -Pc for GMT)\n");
	fprintf(stderr, "    -q\t\tOverride any parameter in the PRM file (-qparam=value)\n");
	fprintf(stderr, "    -Q\t\tDirect elastic deflection mode for a given load file\n");
	fprintf(stderr, "    -h\t\tShow this help (or -hc for clean PRM, -hu for UNIT example)\n");
	fprintf(stderr, "    --help\tOpen the full documentation file\n\n");
	fprintf(stderr, "For full documentation, please read doc/TISC_Documentation.md\n");
	return (1);
}




int surface_processes (float **topo_ant, ModelConfig *cfg, ModelContext *ctx)
{
	/* 
	CALCULATES EROSION AND SEDIMENTATION:
	*/
	bool	switch_horiz_record=false;

	total_sed_mass=total_bedrock_eros_mass=0;
	ctx->total_sed_mass = 0;
	ctx->total_bedrock_eros_mass = 0;
	if (!cfg->erosed_model && !cfg->hydro_model) return (0) ;
	switch_topoest=false;

#ifdef SURFACE_TRANSPORT
	/*Creates a new sediment Block if necessary*/
	if (cfg->erosed_model) {
	    int i;
	    float TimelastBlock=-9999*Matosec;
  	    for (int i=0; i<cfg->Ny; i++) for (int j=0; j<cfg->Nx; j++) eros_now[i][j]=0;
  	    for (int i=0; i<numBlocks; i++)
  		  if (Blocks[i].age > TimelastBlock && Blocks[i].density==denssedim) TimelastBlock = Blocks[i].age;
  	    for (int i=0; i<n_record_times; i++)
  		  if (Time>horiz_record_time[i]-dt/2 && Time<=horiz_record_time[i]+dt/2)
  			  switch_horiz_record=true;
  	    if (Time == Timeini
  	      || ((Time-TimelastBlock)>(dt_record-.001*dt) && dt_record && !n_record_times)
  	      || switch_horiz_record) {
  		  insert_new_Block(numBlocks);
  		  Blocks[numBlocks-1].type = 'S' ;
  		  Blocks[numBlocks-1].density = denssedim ;
  		  Blocks[numBlocks-1].erodibility = erodibility_sed ;
		  Blocks[numBlocks-1].detr_ratio = alloc_matrix(cfg->Ny, cfg->Nx);
		  Blocks[numBlocks-1].detr_grsize = alloc_matrix(cfg->Ny, cfg->Nx);
		  Blocks[numBlocks-1].thickgypsum = alloc_matrix(cfg->Ny, cfg->Nx);
		  Blocks[numBlocks-1].thickhalite = alloc_matrix(cfg->Ny, cfg->Nx);
		  ctx->numBlocks = numBlocks; /* IMPORTANT: Sync Context Struct! */
  	    }
	}


	/*Fluvial Transport: adds to the topo and the next load Dq and removes material from Blocks*/
	Surface_Transport (cfg, ctx, topo_ant, lake_instant_fill);

	/*Diffusive Erosion: adds to the topo and the next load Dq and removes material from Blocks*/
	/*For grids of 100x100 needs dt of .01 My to converge*/
	Diffusive_Eros (cfg, ctx, Kerosdif);

	Landslide_Transport (cfg, ctx, critical_slope);

	/*Adds background erosion and sea sedimentation*/
	constant_rate_eros (cfg, ctx, Keroseol, Ksedim, water_load);
#endif




	{
#define MASS2SEDTHICK_LOCAL(cfg, mass)	((mass) /(cfg->denssedim-sed_porosity*cfg->denswater)/cfg->dx/cfg->dy)
	    int i, j;  float eros_meters, max=-1e19, min=1e19, vol=0;
	    for (i=0; i<cfg->Ny; i++) for (j=0; j<cfg->Nx; j++) {
	    	    eros_meters=MASS2SEDTHICK_LOCAL(cfg, eros_now[i][j])/ctx->dt;
	    	    vol += eros_meters*cfg->dx*cfg->dy;
	    	    if (max<eros_meters)  max=eros_meters;
	    	    if (min>eros_meters)  min=eros_meters;
	    }
	    PRINT_SUMLINE("eros-sed.:  max= %.2f mm/yr\tmin= %.2f mm/yr\tdifference=%.2e %s", max*1e3*secsperyr, min*1e3*secsperyr, vol*secsperyr, "m3/yr");
	}

	{
	    float volume=0, total_vol_seds=0;
    	    for (int i=ctx->numBlocks-1; i>=0; i--) {
    		    int j, k;
    		    for (volume=j=0; j<cfg->Ny; j++) for (k=0; k<cfg->Nx; k++) {
    			    volume += Blocks[i].thick[j][k];
    		    }
    		    volume *= (cfg->dx*cfg->dy);
    		    if (Blocks[i].density==cfg->denssedim) total_vol_seds += volume;
    	    }
    	    PRINT_SUMLINE("sediment_vol: %.2e km3\n", total_vol_seds/1e9);
	}

	return(1);
}




int The_End(ModelConfig *cfg, ModelContext *ctx)
{
	int 	i, j, k;
	char 	command[MAXLENLINE];
	float	volume, total_vol_seds=0, surface;

	fprintf(stdout, "\n\n%d Blocks:", ctx->numBlocks);
	fprintf(stdout, "\nNo. Density Age  Volume  Surf.     Vel_x  Vel_y   Shft_x Shft_y erosL");
	if (cfg->verbose_level>=2) fprintf(stdout, " AgeStop");
	fprintf(stdout, "\n     kg/m3  My   1e3km3  1e3km2    km/My  km/My   km     km     [m] ");
	if (cfg->verbose_level>=2) fprintf(stdout, "  My    ");

	for (i=ctx->numBlocks-1; i>=0; i--) { 
		float vel_x=0, vel_y=0;
		for (volume=surface=j=0; j<cfg->Ny; j++) for (k=0; k<cfg->Nx; k++) {
			volume += Blocks[i].thick[j][k];
			if (Blocks[i].thick[j][k] > 1) surface += cfg->dx*cfg->dy;
		}
		volume *= (cfg->dx*cfg->dy);
		if (Blocks[i].type=='V') {
			for (j=0; j<cfg->Ny; j++) for (k=0; k<cfg->Nx; k++) {
				vel_x += Blocks[i].vel_x[j][k];
				vel_y += Blocks[i].vel_y[j][k];
			}
			vel_x /= (cfg->Nx*cfg->Ny);
			vel_y /= (cfg->Nx*cfg->Ny);
		}
		else {
			vel_x = Blocks[i].vel_x[0][0];
			vel_y = Blocks[i].vel_y[0][0];
		}
		if (Blocks[i].density==cfg->denssedim) total_vol_seds += volume;
		fprintf(stdout, "\n%2d: %5.0f %6.2f %7.1f %6.1f %6.1f %6.1f %6.1f %6.1f   %6.1e", 
			i, Blocks[i].density, Blocks[i].age/Matosec, volume/1e12, surface/1e9, 
			vel_x/1e3*Matosec, vel_y/1e3*Matosec, 
			Blocks[i].shift_x/1e3, Blocks[i].shift_y/1e3, 
			Blocks[i].erodibility); 
		if (cfg->verbose_level>=2) fprintf(stdout, "%7.1f", 
			Blocks[i].time_stop/Matosec);
		fprintf(stdout, "  %c", Blocks[i].type);
		if (i==i_first_Block_load) 	fprintf(stdout, "  1st Block");
	}
	fprintf(stdout, "\n -: %5.0f %6.2f   -       -       0      0      0      0     %6.1e  ", cfg->denscrust, ctx->Timeini/Matosec, erodibility);
	if (cfg->verbose_level>=2) fprintf(stdout, " -   ");
	fprintf(stdout, "-  basement\n");
	fprintf(stdout, "\nFinal total sediment volume: %.2f 1e3 km3\n", total_vol_seds/1e12);

	if (plotting_mode<2) Write_Ouput(cfg, ctx);

	if (cfg->verbose_level>=1) {time_t ltime; time(&ltime); fprintf(stdout, "\nTime end: %s", ctime(&ltime));}

	sprintf(command, "rm -f %s*.tisc.tmp", projectname);
	system(command) ;
	{
		char filename[MAXLENFILE];
		FILE *file;
		/*printf("\nEnvir. Read: %s\n", getenv ("tisc_dir")); system ("echo test: $tisc_dir"); system ("printenv | grep tisc");*/
		sprintf(filename, "%s/.tiscdefaults", getenv ("HOME"));
		if ((file = fopen(filename, "rt")) == NULL) {
			sprintf(command, "echo First use of %s by `whoami` at `hostname` > %s; date >> %s", version, filename, filename);
			if (cfg->verbose_level>=3) fprintf(stdout, "\n%s", command);
			system(command);
			sprintf(command, "mail danielgc@ictja.csic.es < %s", filename);
			if (cfg->verbose_level>=3) fprintf(stdout, "\n%s", command);
			system(command);
			fclose (file);
		}
	}

	if (cfg->verbose_level>=3) AUTHORSHIP;
	fprintf(stdout, "\n");
	free(Lake);

	if (run_type==10 || run_type==2) exit(0);
	return(1);
}




int Viscous_Relaxation(ModelConfig *cfg, ModelContext *ctx)
{
	int 	i, j, NDi=2*cfg->Ny, NDs=2*cfg->Ny, Neqs=cfg->Nx*cfg->Ny;
	double	**A, *b;
	float 	**dwdt;

	if (isost_model!=2 || !Te_default) return(0);

	switch (solver_type) {
	  case 'l': 
		A = alloc_matrix_dbl (Neqs, NDi+1+NDs);
		b = (double *) calloc (Neqs, sizeof(double)); 
		dwdt=alloc_matrix (cfg->Ny, cfg->Nx);
		defineLESalmostdiagonalmatrix(cfg, A, b, q, Dq, w, 1);
		solveLESalmostdiagonal(cfg, A, b, dwdt);
		for (i=0; i<cfg->Ny; i++) for (j=0; j<cfg->Nx; j++)   Dw[i][j] = dwdt[i][j]*ctx->dt;
		free_matrix_dbl (A, Neqs);
		free_matrix (dwdt, cfg->Ny);
		free(b);
		break;
	  case 'm': 
		break;
	}

	for (i=0; i<cfg->Ny; i++) for (j=0; j<cfg->Nx; j++)  w[i][j] += Dw[i][j]; 
	if (switch_topoest) {
		/*Defines the thickness of last infill Block*/
		for (i=0; i<cfg->Ny; i++)  for (j=0; j<cfg->Nx; j++)  Blocks[i_first_Block_load-1].thick[i][j] +=  Dw[i][j] /*MAX(Dw[i][j], 0)*/ ;
	}

	//calculate_topo(topo);

	/*Prints DEFLECTION characteristics*/
	if (cfg->verbose_level>=1) {
		PRINT_GRID_INFO (ctx->topo, "topogr.V ", "m");
		PRINT_GRID_INFO (w, "deflect.V", "m");
	}

	return(1);
}



int Write_Ouput(ModelConfig *cfg, ModelContext *ctx)
{
	//write_file_time(cfg, ctx, w, h_water);
	write_file_deflection(cfg, ctx);
	write_file_Blocks(cfg, ctx);
	write_file_grainsize(cfg, ctx);
	write_file_thicksalt(cfg, ctx);
	write_file_cross_section(cfg, ctx);
	write_file_drainage(cfg, ctx);
	write_file_ice(cfg, ctx);
	write_file_surftransp(cfg, ctx);
	write_file_river_basins(cfg, ctx);
	write_file_lakes(cfg, ctx);
	write_file_velocity_field(cfg, ctx);
	write_file_resume(cfg, ctx);

	/*Make GMT Postscript*/
	if (plotting_mode) {
		char 	command[300];
		if (plotting_mode<=2) {
			sprintf(command, "tisc.gmt.job %s", projectname);
			if (cfg->verbose_level>=3) 
				fprintf(stdout, "\nPostscript file '%s.ps' is going to be produced with command:\n%s\n", projectname, command) ;
			system(command);
		}
		else {
			if (strlen(gif_geom)<2) sprintf(gif_geom, " --elev 30 --fov 140");
			sprintf(command, "tisc.plot.py %s %s; mv -f %s.jpg %s%03d.jpg", projectname, gif_geom, projectname, projectname, n_image);
			if (cfg->verbose_level>=3) 
				fprintf(stdout, "\nPlot files '%s.svg' and %s%03d.jpg to be produced with command:\n%s\n", projectname, projectname, n_image, command) ;
			system(command);
			n_image++;
		}
		if (plotting_mode==2) {
			/*crop by default to the border*/
			if (strlen(gif_geom)<2) sprintf(gif_geom, "-trim -background Khaki -label 'TISC software: %s' -gravity South -append", projectname);
			sprintf(command, "magick -density 300 %s.ps %s -interlace NONE  %s%03d.jpg", /*-fill \"#ffff99\" -draw \"rectangle 8,8 90,25\" -fill \"#000055\" -font helvetica -draw \"text 12,20 t_%+3.2f_My \" */
				projectname, gif_geom, projectname, n_image);
			if (cfg->verbose_level>=3)
				fprintf(stdout, "\n%s\n", command) ;
			system(command);
			n_image++;
		}
	}

	return (1);
}
