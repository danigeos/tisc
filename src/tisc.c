/*******************************************************************************
*****                              TISC                                    *****
********************************************************************************
	For compilation and installation check the file tisc/README
	Main author: Daniel Garcia-Castellanos, danielgc@ictja.csic.es 
	Copyright details and other information in tisc/doc/ 
********************************************************************************

Memory debugging with: 
/usr/local/Cellar/valgrind/3.13.0/bin/valgrind --dsymutil=yes --track-origins=yes --tool=memcheck --leak-check=full tisc linear_range -tf0

	COMMENTS (programmer's agenda)
	-Edit here.
	-Track sediment composition (carbonates, salt, and detrital grain size) in another class in structure Blocks. This to calculate grain size distribution in basin, and as a first step for sed-size dependent erosion, once transitory flow is implemented. 
	-La flexion no es estable con cambios bruscos de Te (Mayo 2001).
	-DONE by M. Berry. Implement sediment load effects on transport and erosion (Sklar). 
	-Implement transitory water flow. Interesting for acceleration of erosion during lake overtopping.
	-DONE. Filter part of the surface water to the lowest surrounding node at 2-cell distance (16 candidates), to simulate underground that accelerates capture. -Alternative: smooth out the discharge grid to simulate underground flow.
	-DONE. Solved bug in erosed_model 6.
	-DONE. Implement sediment compaction (in calculate_topo and when writting the hrz/pfl files).
	-DONE. Add the water layer to the pfl profile, as in tAo.
	-DONE. Find water divides: the maximum swimming distance with the neigbouring cells. In write_file_drainage.
	-DONE. It should take a mean erodibility when eroding.
	-DONE, works if hydro_model!=0. Water load of lakes does not work properly. 
*/


#include "tisc.h"
#include "tisclib.c"
#include "tiscio.c"
#include "../lib/libreria.h"

int main(int argc, char **argv)
{

	/*get input parameters and files*/
	inputs(argc, argv) ;

	fprintf(stdout, "\nT= %.4f My", Time/Matosec);

	if (switch_dt_output) {calculate_topo(topo); Write_Ouput();}

	/*MAIN LOOP: In this loop time increments from Timeini to Timefinal*/
	do {
		/*Remember topography before tectonics and flexure*/
		calculate_topo(topo);

		/*Calculate tectonic deformation and tectonic load*/
		tectload();
		
		/*Sea level variations*/
		calculate_sea_level();

		/*Calculates water column load*/
		calculate_water_load();

		/*Define & solve elastic flexure equation*/
		Elastic_Deflection();

		/*Define & solve viscoelastic flexure equation*/
		Viscous_Relaxation();

		/*Calculates surface water flow and sediment-erosion*/
		surface_processes(topo);

		Time += dt;
		fprintf(stdout, "\nT= %.4f My", Time/Matosec);

		if (switch_dt_output) Write_Ouput();
	} while (Time < Timefinal-dt/10);

	The_End();

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
	BOOL	success_def_prm=NO, switch_initial_geom=NO;

	run_type=0;
	solver_type = 'l';
	setbuf(stdout, NULL);

	putenv("tisc_dir=" TISCDIR); /*printf("\nEnvir. Read: %s", getenv ("tisc_dir")); system ("echo test: $tisc_dir"); system ("printenv | grep tisc");*/

	/*Version of TISC is matched against the parameters file *.PRM*/
	/*¡¡ UPDATE template.PRM !!*/
	strcpy(version, "TISC_2019-05-03");

	/*Default parameter values are read from ./tisc/doc/template.PRM:*/
	sprintf(projectname, "%s/doc/template", TISCDIR);
	success_def_prm = read_file_parameters(0, 0);
	sprintf(projectname, "");

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
					if (argv[iarg][2]) reformat = value;
					break;
				case 'F':
					run_type=2;
					if (strlen(prm)>0) strcpy (resume_filename, prm);
					else sprintf(resume_filename, "%s"".all", projectname);
					break;
				case 'h':
					switch (argv[iarg][2]) {
						case 'p':
							fprintf(stderr, "\nFile ./tisc/doc/template.PRM (sample parameters file) follows in stdout:\n") ;
							sprintf(command, "cat %s/doc/template.PRM", TISCDIR);
							system(command) ;
							break;
						case 'c':
							fprintf(stderr, "\nFile ./tisc/doc/template.PRM (sample parameters file) follows in stdout:\n") ;
							sprintf(command, "cat %s/doc/template.PRM | %s/script/cleanPRM", TISCDIR, TISCDIR);
							system(command) ;
							break;
						case 'u':
							fprintf(stderr, "\nFile ./tisc/doc/template.UNIT (sample unit file) follows in stdout:\n") ;
							sprintf(command, "cat %s/doc/template.UNIT", TISCDIR);
							system(command) ;
							break;
						default:
							fprintf(stderr, "\nFile ./tisc/doc/tisc.info.txt follows:\n") ;
							sprintf(command, "more %s/doc/tisc.info.txt", TISCDIR);
							AUTHORSHIP;
							system(command) ;
							break;
					}
					fprintf(stderr, "\n") ;
					exit(0);
				case 'Q':
					run_type=1;
					strcpy(load_file_name, prm);
					break;
				case 'V':
					verbose_level = 1;
					if (argv[iarg][2]) verbose_level = value;
					break;
			}
		}
		else {
			if (run_type != 2) run_type=10;
			if (strlen(projectname)<1) strcpy(projectname, argv[iarg]);
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
		exit(0);
	}

	nloads=0; n_image=0; nlakes=0;
	numBlocks=0; i_first_Block_load=0; i_Block_insert=0;
	nwrotenfiles=0; switch_topoest=NO;

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
			if (switch_dt_output) n_image--; /*Don't produce 2 jpg's of the same stage of restart*/
			return(1);
		case 10:
			if (!read_file_parameters(verbose_level>=1, 0)) {
				syntax();
				fprintf(stderr, "\n\tAvailable parameter files in this directory:\n");
				system("ls *.PRM");
				if (!success_def_prm) {
					PRINT_ERROR("\aDefault parameters file './tisc/doc/template.PRM' could not be read.\n"); 
				}
				exit(0);
			}
			interpr_command_line_opts(argc, argv);
			if (reformat) {
				sprintf(projectname, "%s/doc/template", TISCDIR);
				read_file_parameters(0, reformat);
				exit(0);
			}
			break;
	}

	{
		char filename[MAXLENLINE]; FILE *file;
		sprintf(filename, "%s.out", projectname);
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
	dt *= Matosec;
	dt_eros *= Matosec;
	tau *= Matosec;
	dt_record *= Matosec;
	Timefinal *= Matosec;
	Timeini *= Matosec;
	Time = Timeini;
	last_time_file_time = Timeini - 9999*dt_record;	/*very old*/
	Kerosdif *= 1e6/Matosec;
	Keroseol /= Matosec;
	Ksedim   /= Matosec;
	rain *= 1e6/Matosec/1e3;
	if (hydro_model==1) Krain *= 1e6/Matosec/1e3/1e3;
	evaporation_ct *= 1e6/Matosec/1e3;
	lost_rate *= 1e-2 * 1e-3;
	temp_sea_level += TEMP_FREEZE_WATER; /*converts from C to K*/
	switch_write_file_Blocks = YES; 
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
	if (switch_ps && !switch_write_file){	 				PRINT_WARNING("switch_write_file should be turned on. Postscript may not be produced.") ; }
	if (tau<=0 && isost_model==2)  		{ isost_model=1; 	PRINT_WARNING("Negative tau will be ignored: elastic plate is assumed.");}
	if (Nx<6 || Ny<6)					{					PRINT_WARNING("Too coarse %dx%d gridding.", Nx, Ny); }

	sprintf(command, "rm -f .*.tisc.tmp %s.mtrz", projectname);
	system(command);

	read_file_sea_level(); calculate_sea_level();
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

	PRINT_INFO("Enetering command line interpretation.");
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
					strcpy(boundary_conds, prm);
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
					switch_file_out=YES;
					break;
				case 'P':
					switch_ps=YES;
					switch_write_file_Blocks=YES;
					if (argv[iarg][2] == 'c') {
						switch_dt_output=YES;
						strcpy(gif_geom, "");
						if (strlen(prm2)>0) strcpy(gif_geom, prm2);
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
					if (argv[iarg][2]) verbose_level = value;
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

	/*Solves flexure problem in direct mode: taking a single external load file and 
	writing deflection in standar output. There are no other input files neither 
	output files*/

	PRINT_INFO("Entering direct mode. xmin/xmax/ymin/ymax=%.1f/%.1f/%.1f/%.1f Nx/Ny=%d/%d\n", xmin,xmax,ymin,ymax, Nx,Ny);
	dx = (xmax-xmin) / (Nx-1) ;
	dy = (ymax-ymin) / (Ny-1) ;
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
	Elastic_Deflection();
	fprintf(stdout, "\n\nx[km]\t\ty[km]\t\tw[m]\t\tpressure[Pa]\n"); 
	for (i=0; i<Ny; i++) for (j=0; j<Nx; j++) 
		fprintf(stdout, "%8.1f\t%8.1f\t%8.1f\t%8.1f\n", 
			(xmin+j*dx)/1e3, (ymax-i*dy)/1e3, w[i][j], h_last_unit[i][j]);
	fprintf(stdout, "\n"); 
	return (1);
}




int tectload()
{
	/*
	CALCULATES NEW LOAD INCREMENT FROM UNIT FILES, Returns 1 if elastic
	flexure must be done (i.e, if changes in load  occurred), 0 otherwise.
	*/
	
	PRINT_GRID_INFO (topo, "topogr.  ", "m");

	/*Moves Blocks*/
	move_Blocks();

	/*Reads external load from file*/
	while (read_file_unit());

	/*Distributes the emplacement of a unit along time*/
	gradual_Block();

	Repare_Blocks();

	return (1);
}



int Elastic_Deflection()
{
	int 	i, j, NDi=2*Ny, NDs=2*Ny, Neqs=Nx*Ny, 
		nonzeroes=13*Nx*Ny, ESP, PATH, FLAG, NSP, 
		*R, *C, *IC, *IA, *JA, *ISP;
	double	**A, *b, *w_aux;
	float 	*B, *Z, *mathlib_matrix, *RSP;
	BOOL	load_changes=NO;

	for (i=0; i<Ny; i++) for (j=0; j<Nx; j++) if (Dq[i][j]) load_changes = YES;
	if (isost_model>0 && (load_changes || (Time==Timeini && (Px || Py || Pxy)))) {
    	    if (!Te_default) {
    		/*LOCAL ISOSTASY*/
    		float Krest;
    		for (i=0; i<Ny; i++) for (j=0; j<Nx; j++)  {
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
    		    defineLESalmostdiagonalmatrix(A, b, q, Dq, w, 0);
    		    solveLESalmostdiagonal(A, b, Dw);
    		    free_matrix_dbl (A, Neqs);
    		    free(b);
    		    break;
    		case 'm':
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
    		    defineLESmatrix_for_mathlib(mathlib_matrix, IA, JA, B, Dq, w, &nonzeroes, 0);
    		    for (i=0; i<Neqs; i++)  R[i]=IC[i]=C[i]=i+1;
#ifdef MATHLIB_SOLVER
    		    printf("\nNeqs=%d NSP=%d nonzeroes=%d", Neqs, NSP, nonzeroes);
    		    cdrv(&Neqs, R, C, IC, IA, JA, mathlib_matrix, B, Z, &NSP, ISP, RSP, &ESP, &PATH, &FLAG);
    		    if (FLAG!=0) {PRINT_ERROR("\aTDRV exit #%d. Memory excess=%d\n", FLAG, ESP); }
    		    fprintf(stdout, "\tMemory excess=%d\n", ESP);
#endif
    		    for (i=0; i<Ny; i++) for (j=0; j<Nx; j++)	Dw[i][j] = Z[j*Ny+i];

    		    free(Z); free(B); free(mathlib_matrix);
    		    free(R); free(C);free(IC); free(IA); free(JA); free(ISP); free(RSP);
    		    break;
    	      }
	    }

    	    for (i=0; i<Ny; i++) for (j=0; j<Nx; j++)  w[i][j] += Dw[i][j];
    	    if (switch_topoest) {
   		    /*Defines the thickness of last infill Block*/
    		    for (i=0; i<Ny; i++)  for (j=0; j<Nx; j++)  
		    	Blocks[i_first_Block_load-1].thick[i][j] +=  MAX_2(Dw[i][j], 0) ;
    	    }


    	    calculate_topo(topo);

    	    /*Statistics on load*/
    	    {
    		    int i, j, imax, jmax, imin, jmin;
    		    float maxloadinc=-1e25, minloadinc=1e25, total_load=0;
    		    for (i=0;i<Ny;i++) for (j=0;j<Nx;j++) {
    			q[i][j] += Dq[i][j];
    			total_load += Dq[i][j]*dx*dy;
    			if (maxloadinc<Dq[i][j]) {maxloadinc=Dq[i][j]; imax=i; jmax=j;}
    			if (minloadinc>Dq[i][j]) {minloadinc=Dq[i][j]; imin=i; jmin=j;}
    		    }
    		    if (load_changes) 
		    	PRINT_SUMLINE("load  now:  max = %+8.2e N/m2  min = %+8.2e N/m2   Total: %+8.2e N\tMax at %.0f,%.0f km", 
		    		maxloadinc, minloadinc, total_load, (xmin+dx*jmax)/1e3, (ymax-dy*imax)/1e3); 
    	    }
    	    /*Statistics on deflection*/
    	    {
    		    float total_load=0, total_restitutive_force=0, Krest;
    		    for (i=0; i<Ny; i++) for (j=0; j<Nx; j++)	{
    			    total_load += Dq[i][j];
    			    GET_KREST(Krest, q, i,j);
    			    total_restitutive_force += (Krest*w[i][j]);
    		    }
    		    if (verbose_level>=1) {
    			    PRINT_SUMLINE("load: %+8.2e N   restit_force: %+8.2e N", total_restitutive_force*dx*dy, total_load*dx*dy);
			    PRINT_GRID_INFO (w, "deflect. ", "m");
    		    }
    	    }
	}

	/*Resets deflection and load grids*/
	for (i=0; i<Ny; i++)  for (j=0; j<Nx; j++)  Dq[i][j]=Dw[i][j]=0;

	return(1);
}





int move_Blocks()
{
	int	*nshift_x, *nshift_y;
	float	**new_thick;
	char 	tmpTSBCfilename[84]="";

	/*
	  Moves the Blocks and calculates the isostatic load and thickness change at each cell.
	  Deforms the sediment above.
	*/

	new_thick = alloc_matrix(Ny, Nx);
	nshift_x = calloc(numBlocks, sizeof(int));
	nshift_y = calloc(numBlocks, sizeof(int));

	for (int iu=0; iu<numBlocks; iu++) {
		if (Blocks[iu].density == denssedim) {
/*!!*/
//Blocks[iu].vel_y[0][0] = 2.5e3/Matosec;
//Blocks[iu].last_vel_time = Blocks[iu].age;
//Blocks[iu].last_shift_x  = 0;
//Blocks[iu].last_shift_y  = 0;
//Blocks[iu].time_stop     = 1e19;
		/*DEFORM SEDIMENT Blocks*/
		for (int i=0; i<Ny; i++) for (int j=0; j<Nx; j++) 
			new_thick[i][j] = Blocks[iu].thick[i][j];
		for (int i=0; i<Ny; i++) for (int j=0; j<Nx; j++){
			float sedthick;
			sedthick=Blocks[iu].thick[i][j];
			/*Find the uppermost moving Block below this point in this sedim. Block*/
			for (int ju=iu-1; ju>=0; ju--) {
			/*Calculate the thickness of sediments between the top of this sed. Block and the moving Block*/
			if (Blocks[ju].density == denssedim) {
				sedthick += Blocks[ju].thick[i][j];
			}
			else {
   			  /*Amount of cells to propagate the deformation: ~20 deg assumed.*/
			  int nprop_x = SIGN(nshift_x[ju]) * (int) ceil(sedthick*3/dx);
			  int j_unprop = j-nprop_x;
			  int nprop_y = SIGN(nshift_y[ju]) * (int) ceil(sedthick*3/dy);
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
		for (int i=0; i<Ny; i++) for (int j=0; j<Nx; j++) {
			Dq[i][j] += g * (new_thick[i][j] - Blocks[iu].thick[i][j]) * Blocks[iu].density;
			if (new_thick[i][j]<-1) PRINT_ERROR("negative sediment thickness: %.1f m", new_thick[i][j]);
			Blocks[iu].thick[i][j] = new_thick[i][j];
		}
		}
		else //!!
		{
		  if (Blocks[iu].type == 'V' && Time < Blocks[iu].time_stop) {
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
		elapsed_time = Time-Blocks[iu].last_vel_time;
		dt_aux = dt;
		Lx = xmax-xmin; Ly = ymax-ymin;
		n=Nx-1; m=Ny-1; nn=Nx*Ny; nincogn=2*nn; nbanda=4*(Nx)+7;
		viscTer        =	alloc_array_dbl(Nx*Ny);
		viscosity      =	alloc_array_dbl(Nx*Ny);
		average_pressure=	alloc_array_dbl(Nx*Ny);
		vel_x_array    =	alloc_array_dbl(Nx*Ny);
		vel_y_array    =	alloc_array_dbl(Nx*Ny);
		vert_strain_rate =	alloc_array_dbl(Nx*Ny);
		layer_thickness =	alloc_array_dbl(Nx*Ny);
		sprintf(tmpTSBCfilename, "%s"".TSBC.tmp", projectname);
		reformat_file_thin_sheet_BC(tmpTSBCfilename);
		PRINT_INFO("Deforming thin_sheet Block %d", iu);
		for (i=0; i<Ny; i++) for (j=0; j<Nx; j++)  {
			int k=(Ny-1-i)*Nx+j, l, m;
			viscTer[k] = Blocks[iu].viscTer[i][j];  /*.1e7*/  /* viscTer[Pa] = 1/2 * strength[Pa*m] / layer.thickness[m] */
			/*AVERAGE PRESSURE*/
			/*Water term*/
			average_pressure[k] = g*denswater*h_water[i][j];
			/*Blocks terms from the top to the thin_sheet*/
			for (l=numBlocks-1; l>iu; l--) {
				average_pressure[k] += g*Blocks[l].density*Blocks[l].thick[i][j];
			}
			/*Thin sheet term*/
			average_pressure[k] += g/2*Blocks[iu].density*Blocks[iu].thick[i][j];
			/*Asthenosphere restitutive push*/
			average_pressure[k] += (densasthen-densenv)*g*w[i][j];
			/*change sign to <0, meaning compresion*/
			average_pressure[k] *= -1; 

			vert_strain_rate[k] = 0;
			layer_thickness[k] = Blocks[iu].thick[i][j];
			if (elapsed_time > dt) {
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
			for (i=0; i<Ny; i++) for (j=0; j<Nx; j++)  {
				int k=(Ny-1-i)*Nx+j;
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
		for (i=0; i<Ny; i++)  for (j=0; j<Nx; j++)  {
			int k=(Ny-1-i)*Nx+j;
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
			theor_shift_x = Blocks[iu].vel_x[0][0] * (Time-Blocks[iu].last_vel_time);
			theor_shift_y = Blocks[iu].vel_y[0][0] * (Time-Blocks[iu].last_vel_time);
			nshift_x[iu] = floor((theor_shift_x - Blocks[iu].last_shift_x) /dx +.5);
			nshift_y[iu] = floor((theor_shift_y - Blocks[iu].last_shift_y) /dy +.5);
			if (Time > Blocks[iu].time_stop + .1*dt) {nshift_x[iu]=0; nshift_y[iu]=0;}
			Blocks[iu].shift_x += nshift_x[iu]*dx;
			Blocks[iu].shift_y += nshift_y[iu]*dy;
			Blocks[iu].last_shift_x += nshift_x[iu]*dx;
			Blocks[iu].last_shift_y += nshift_y[iu]*dy;
			for (i=0; i<Ny; i++) for (j=0; j<Nx; j++){
				i_unshifted = i+nshift_y[iu];	j_unshifted = j-nshift_x[iu];	
				DOMAIN_LIMIT(i_unshifted, j_unshifted);	/*[i][j], unshifted*/
				new_thick[i][j] = Blocks[iu].thick[i_unshifted][j_unshifted];
			}
			for (i=0; i<Ny; i++) for (j=0; j<Nx; j++) {
				Dq[i][j] += g * (new_thick[i][j] - Blocks[iu].thick[i][j]) * Blocks[iu].density;
				Blocks[iu].thick[i][j] = new_thick[i][j];
			}
		  }
		}
	}

	free_matrix(new_thick, Ny);
	free(nshift_x);
	free(nshift_y);
	return(1);
}




int read_file_unit()
{
	/*
	  READS UNIT FILE NAMED 'projectnameNUM.UNIT' WHERE 'NUM' IS 1 FOR THE 
	  FIRST UNIT, 2 FOR THE SECOND, ETC. Interpolates this input. Creates
	  new unit to store its properties and cuts sediment units when file
	  contains fault depth rather than a thickness itself.
	*/

	int 	nparams=0;
	float	time_stop=9999/*My*/, time_unit, 
		erodibility_aux=NO_DATA, fill_up_to=NO_DATA, 
		vel_x=0, vel_y=0, density=NO_DATA;
	BOOL 	insert, cut_Blocks, cut_all, top, fault, switch_move, 
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
		show=(verbose_level>=3)? 1 : 0;
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
	if (time_unit>Time+.1*dt || time_unit<Timeini) return(0);

	PRINT_INFO("Reading '%s'", filename);
	switch_move = fault = switch_gradual = 
		insert = hidden = cut_Blocks = cut_all = 
		thin_sheet = top = ride = z_absol = NO;
	i_Block_insert = numBlocks;

	/*READS AND INTERPOLATES UNIT/LOAD FILE*/
	{
		int nlines=0, nread, show, replace=0;
		char str1[MAXLENLINE], str2[MAXLENLINE], line[MAXLENLINE+200], *lineptr;
		show=(verbose_level>=3)? 1 : 0;
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
				Match_Param_Replace_int ( "densenv",		densenv,   	0 )
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
		readinterp2D(file, h_last_unit, mode_interp, 0, xmin, xmax, ymin, ymax, Nx, Ny);
	else {
		for (int i=0; i<Ny; i++) for (int j=0; j<Nx; j++) h_last_unit[i][j] = MAX_2(0, fill_up_to-topo[i][j]);
	}
	fclose(file);

	nloads++;

	vel_x *= 1e3/Matosec;
	vel_y *= 1e3/Matosec;
	time_stop *= Matosec;

	/*ACT ACCORDING TO THE SIGNALS*/
	if (thin_sheet) {
		switch_move = YES; 
	}
	if (fault) {
		switch_move = YES;
	}

	/*Check incompatibilities between unit file signals*/
	if (switch_gradual && switch_move) {
		PRINT_WARNING("Gradual & moving Blocks are not compatible. This one won't be gradual.");
		switch_gradual = NO;
	}
	if (switch_gradual && hidden) {
		PRINT_WARNING("Gradual & hidden Blocks are not compatible. This one won't be hidden.");
		hidden = NO;
	}

	/*Creates a Block of infill if switch_topoest; it will be filled later during Deflection*/
	if (switch_topoest) {
		insert_new_Block(i_first_Block_load);
		Blocks[i_first_Block_load].type = 'I'; 	/*stands for Infill*/
		Blocks[i_first_Block_load].density = densinfill;
		if (densinfill<2550) Blocks[i_first_Block_load].erodibility = erodibility_sed;
		i_first_Block_load++; 	i_Block_insert++;
	}
	if (insert) {
		i_Block_insert = 0;
	}
	if (top) {
		for (int k=numBlocks-1; k>=0; k--) {
			if (Blocks[k].density != denssedim) {
				i_Block_insert = k+1;
				break;
			}
		}
	}
	if (cut_all) {
		cut_Blocks = YES;
	}

	if (fault && !top) i_Block_insert = 0;

	if (fault) {
		int numBlocks0=numBlocks;
		/*CUT BlockS*/
		/*Make copies of all Blocks*/
		PRINT_DEBUG("Cutting Blocks: numBlocks= %d", numBlocks);
		for (int k=0; k<numBlocks0; k++) {
			float **thick_aux, **vel_x_aux, **vel_y_aux;
			insert_new_Block(numBlocks);
			thick_aux = Blocks[numBlocks-1].thick;
			vel_x_aux = Blocks[numBlocks-1].vel_x;
			vel_y_aux = Blocks[numBlocks-1].vel_y;
			Blocks[numBlocks-1] = Blocks[k];
			Blocks[numBlocks-1].thick = thick_aux;
			Blocks[numBlocks-1].vel_x = vel_x_aux;
			Blocks[numBlocks-1].vel_y = vel_y_aux;
			Blocks[numBlocks-1].vel_x[0][0] = vel_x;
			Blocks[numBlocks-1].vel_y[0][0] = vel_y;
			Blocks[numBlocks-1].last_vel_time = Time;
			Blocks[numBlocks-1].last_shift_x = 0;
			Blocks[numBlocks-1].last_shift_y = 0;
			Blocks[numBlocks-1].time_stop = time_stop;
			if (Blocks[numBlocks-1].type == 'V') {
			    if (thin_sheet) {
				Blocks[numBlocks-1].vel_x  = alloc_matrix(Ny, Nx);
				Blocks[numBlocks-1].vel_y  = alloc_matrix(Ny, Nx);
				Blocks[numBlocks-1].visc   = alloc_matrix(Ny, Nx);
				Blocks[numBlocks-1].viscTer= alloc_matrix(Ny, Nx);
			    }
			    else {
				Blocks[numBlocks-1].type = '-';
			    }
			}
			if (Blocks[numBlocks-1].type == 'S') {
				Blocks[numBlocks-1].detr_ratio  = alloc_matrix(Ny, Nx);
				Blocks[numBlocks-1].detr_grsize = alloc_matrix(Ny, Nx);
			}
			if (density         != NO_DATA && Blocks[numBlocks-1].type != 'S') Blocks[numBlocks-1].density = density;
			if (erodibility_aux != NO_DATA && Blocks[numBlocks-1].type != 'S') Blocks[numBlocks-1].erodibility = erodibility_aux;
		}
		PRINT_DEBUG("Updating Blocks_base: numBlocks= %d", numBlocks);
		/*Modify Blocks_base and cut above the fault*/
		for (int i=0; i<Ny; i++) for (int j=0; j<Nx; j++) {
			float z_fault=-h_last_unit[i][j], base_of_Block=Blocks_base[i][j];
			if (z_absol) base_of_Block -= w[i][j];
			h_last_unit[i][j] = MAX_2(0, Blocks_base[i][j]-z_fault);	/*Block thickness below fault (to create the new Block, see below)*/
			Blocks_base[i][j]  = MIN_2(Blocks_base[i][j], z_fault);	/*new base of Blocks*/
			if (cut_Blocks) {
				for (int k=0; k<numBlocks0; k++) {
					float top_of_Block=base_of_Block+Blocks[k].thick[i][j];
					if (Blocks[k].density == denssedim && !cut_all) {
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

	/*Add the thickness in file to the new Block; Shrink the Blocks and basement if the thickness is negative*/
	if (!switch_gradual && !hidden) {
		for (int i=0; i<Ny; i++)  for (int j=0; j<Nx; j++) {
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
	if (switch_gradual) Blocks[numBlocks-1].type = 'G';
	if (thin_sheet) {
		float default_viscTerm = .5e7; /*.1e7*/
		char filename[MAXLENFILE];
		Blocks[i_Block_insert].type    = 'V';
		Blocks[i_Block_insert].vel_x   = alloc_matrix(Ny, Nx);
		Blocks[i_Block_insert].vel_y   = alloc_matrix(Ny, Nx);
		Blocks[i_Block_insert].visc    = alloc_matrix(Ny, Nx);
		Blocks[i_Block_insert].viscTer = alloc_matrix(Ny, Nx);
		sprintf(filename, "%s%d.VISC", projectname, nloads);
		if ((file = fopen(filename, "rt")) == NULL) {
			PRINT_WARNING("Cannot read thermal viscosity file '%s'.", filename);
			for (int i=0; i<Ny; i++)  for (int j=0; j<Nx; j++) Blocks[i_Block_insert].viscTer[i][j] = default_viscTerm;
		}
		else {
			PRINT_INFO("Reading viscosity for unit %d from file '%s'.", nloads, filename);
			readinterp2D(file, Blocks[i_Block_insert].viscTer, mode_interp, default_viscTerm, xmin, xmax, ymin, ymax, Nx, Ny);
		}
	}
	Blocks[i_Block_insert].density = density;
	Blocks[i_Block_insert].erodibility = erodibility_aux;
	Blocks[i_Block_insert].vel_x[0][0] = vel_x;
	Blocks[i_Block_insert].vel_y[0][0] = vel_y;
	Blocks[i_Block_insert].time_stop = time_stop;

	if (ride) {
		for (int i_Block=i_Block_insert+1; i_Block<numBlocks; i_Block++) {
			Blocks[i_Block].vel_x         = Blocks[i_Block_insert].vel_x; 
			Blocks[i_Block].vel_y         = Blocks[i_Block_insert].vel_y; 
			Blocks[i_Block].last_shift_x  = 0; 
			Blocks[i_Block].last_shift_y  = 0; 
			Blocks[i_Block].last_vel_time = Time; 
			Blocks[i_Block].time_stop     = Blocks[i_Block_insert].time_stop; 
		}
	}

	/*Don't Repare_Blocks() in case of: 
		Gradual load, because then h_last_unit[] will be empty until tectload()
		Topoest load, because the infill Block will be filled upon deflection.
	*/
	if (!switch_gradual && !switch_topoest) Repare_Blocks();

	/*Increment the isostatic load for this time interval*/
	if (!switch_gradual && !fault) 
		for (int i=0; i<Ny; i++) for (int j=0; j<Nx; j++) Dq[i][j] += (density-densenv)*g*h_last_unit[i][j];

	PRINT_INFO("Unit read from '%s'. ", filename);
	PRINT_DEBUG("%d params; dens=%.0f kg/m3; erodibility=%.1e; ", nparams, density, erodibility_aux);
	if (!fault) {
		float load=0;
		for (int i=0; i<Ny; i++) for (int j=0; j<Nx; j++) load += (density-densenv)*h_last_unit[i][j];
		fprintf(stdout, "%.2e kg. ", load*dx*dy);
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
		Displays the command line syntax of the program
	*/
	char 	filename[MAXLENFILE], line[MAXLENLINE], *lineptr;
	FILE 	*file;
	BOOL 	print=NO;

	sprintf(filename, "%s/doc/tisc.info.txt", TISCDIR);
	if ((file = fopen(filename, "rt")) == NULL) {
		PRINT_WARNING("Cannot find file %s. You shouldn't have moved TISC directory after compilation...", filename);
		return(0);
	}
	fprintf(stderr, "\nSyntax:\n");
	while (1) {
		fgets(line, MAXLENLINE-1, file);
		if ((lineptr=strstr(line, "2.- INTRO"))) break;;
		if (print==YES) fprintf(stderr, "%s", line);
		if (strstr(line, "SYNTAX")) print=YES;
	}
	fclose(file);
	return (1);
}




int surface_processes (float **topo_ant)
{
	/* 
	CALCULATES EROSION AND SEDIMENTATION:
	*/
	BOOL	switch_horiz_record=NO;
	float	TimelastBlock;
	int 	test;

	total_sed_mass=total_bedrock_eros_mass=0;
	if (!erosed_model && !hydro_model) return (0) ;
	switch_topoest=NO;

#ifdef SURFACE_TRANSPORT
	/*Creates a new sediment Block if necessary*/
	if (erosed_model) {
	    int i;
	    float TimelastBlock=-9999*Matosec;
  	    for (int i=0; i<Ny; i++) for (int j=0; j<Nx; j++) eros_now[i][j]=0;
  	    for (int i=0; i<numBlocks; i++)
  		  if (Blocks[i].age > TimelastBlock && Blocks[i].density==denssedim) TimelastBlock = Blocks[i].age;
  	    for (int i=0; i<n_record_times; i++)
  		  if (Time>horiz_record_time[i]-dt/2 && Time<=horiz_record_time[i]+dt/2)
  			  switch_horiz_record=YES;
  	    if (Time == Timeini
  	      || ((Time-TimelastBlock)>(dt_record-.001*dt) && dt_record && !n_record_times)
  	      || switch_horiz_record) {
  		  insert_new_Block(numBlocks);
  		  Blocks[numBlocks-1].type = 'S' ;
  		  Blocks[numBlocks-1].density = denssedim ;
  		  Blocks[numBlocks-1].erodibility = erodibility_sed ;
		  Blocks[numBlocks-1].detr_ratio = alloc_matrix(Ny, Nx);
		  Blocks[numBlocks-1].detr_grsize = alloc_matrix(Ny, Nx);
  	    }
	}


	/*Fluvial Transport: adds to the topo and the next load Dq and removes material from Blocks*/
	Surface_Transport (topo, topo_ant, dt, dt_eros, erosed_model, lake_instant_fill);

	/*Diffusive Erosion: adds to the topo and the next load Dq and removes material from Blocks*/
	/*For grids of 100x100 needs dt of .01 My to converge*/
	Diffusive_Eros (Kerosdif, dt, dt_eros/5);

	Landslide_Transport (critical_slope, dt, dt_eros);

	/*Adds background erosion and sea sedimentation*/
	constant_rate_eros (topo, Keroseol, Ksedim, sea_level, water_load, dt, Time);
#endif




	{
#define MASS2SEDTHICK(mass)	((mass) /(denssedim-sed_porosity*denswater)/dx/dy)	/*converts sediment mass into sediment thickness*/
	    int i, j;  float eros_meters, max=-1e19, min=1e19, vol=0;
	    for (i=0; i<Ny; i++) for (j=0; j<Nx; j++) {
	    	    eros_meters=MASS2SEDTHICK(eros_now[i][j])/dt;
	    	    vol += eros_meters*dx*dy;
	    	    if (max<eros_meters)  max=eros_meters;
	    	    if (min>eros_meters)  min=eros_meters;
	    }
	    PRINT_SUMLINE("eros-sed.:  max= %.2f mm/yr\tmin= %.2f mm/yr\tdifference=%.2e %s", max*1e3*secsperyr, min*1e3*secsperyr, vol*secsperyr, "m3/yr");
	}

	{
	    float volume=0, total_vol_seds=0;
    	    for (int i=numBlocks-1; i>=0; i--) {
    		    int j, k;
    		    for (volume=j=0; j<Ny; j++) for (k=0; k<Nx; k++) {
    			    volume += Blocks[i].thick[j][k];
    		    }
    		    volume *= (dx*dy);
    		    if (Blocks[i].density==denssedim) total_vol_seds += volume;
    	    }
    	    PRINT_SUMLINE("sediment_vol: %.2e km3\n", total_vol_seds/1e9);
	}

	return(1);
}




int The_End()
{
	int 	i, j, k, action=2, status;
	char 	command[MAXLENLINE];
	float	volume, total_vol_seds=0, surface;

	fprintf(stdout, "\n\n%d Blocks:", numBlocks);
	fprintf(stdout, "\nNo. Density Age  Volume  Surf.     Vel_x  Vel_y   Shft_x Shft_y erosL");
	if (verbose_level>=2) fprintf(stdout, " AgeStop");
	fprintf(stdout, "\n     kg/m3  My   1e3km3  1e3km2    km/My  km/My   km     km     [m] ");
	if (verbose_level>=2) fprintf(stdout, "  My    ");

	for (i=numBlocks-1; i>=0; i--) { 
		float vel_x=0, vel_y=0;
		for (volume=surface=j=0; j<Ny; j++) for (k=0; k<Nx; k++) {
			volume += Blocks[i].thick[j][k];
			if (Blocks[i].thick[j][k] > 1) surface += dx*dy;
		}
		volume *= (dx*dy);
		if (Blocks[i].type=='V') {
			for (j=0; j<Ny; j++) for (k=0; k<Nx; k++) {
				vel_x += Blocks[i].vel_x[j][k];
				vel_y += Blocks[i].vel_y[j][k];
			}
			vel_x /= (Nx*Ny);
			vel_y /= (Nx*Ny);
		}
		else {
			vel_x = Blocks[i].vel_x[0][0];
			vel_y = Blocks[i].vel_y[0][0];
		}
		if (Blocks[i].density==denssedim) total_vol_seds += volume;
		fprintf(stdout, "\n%2d: %5.0f %6.2f %7.1f %6.1f %6.1f %6.1f %6.1f %6.1f   %6.1e", 
			i, Blocks[i].density, Blocks[i].age/Matosec, volume/1e12, surface/1e9, 
			vel_x/1e3*Matosec, vel_y/1e3*Matosec, 
			Blocks[i].shift_x/1e3, Blocks[i].shift_y/1e3, 
			Blocks[i].erodibility); 
		if (verbose_level>=2) fprintf(stdout, "%7.1f", 
			Blocks[i].time_stop/Matosec);
		fprintf(stdout, "  %c", Blocks[i].type);
		if (i==i_first_Block_load) 	fprintf(stdout, "  1st Block");
	}
	fprintf(stdout, "\n -: %5.0f %6.2f   -       -       0      0      0      0     %6.1e  ", denscrust, Timeini/Matosec, erodibility);
	if (verbose_level>=2) fprintf(stdout, " -   ");
	fprintf(stdout, "-  basement\n");
	fprintf(stdout, "\nFinal total sediment volume: %.2f 1e3 km3\n", total_vol_seds/1e12);

	if (!switch_dt_output) Write_Ouput();

	if (verbose_level>=1) {time_t ltime; time(&ltime); fprintf(stdout, "\nTime end: %s", ctime(&ltime));}

	sprintf(command, "rm -f %s*.tisc.tmp", projectname);
	system(command) ;
	{
		char filename[MAXLENFILE];
		FILE *file;
		/*printf("\nEnvir. Read: %s\n", getenv ("tisc_dir")); system ("echo test: $tisc_dir"); system ("printenv | grep tisc");*/
		sprintf(filename, "%s/.tiscdefaults", getenv ("HOME"));
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
	fprintf(stdout, "\n");
	free(Lake);

	if (run_type==10 || run_type==2) exit(0);
	return(1);
}




int Viscous_Relaxation()
{
	int 	i, j, NDi=2*Ny, NDs=2*Ny, Neqs=Nx*Ny, 
		nonzeroes=13*Nx*Ny, nmax, ESP, FLAG, 
		*R, *IC, *IA, *JA, *ISP;
	double	**A, *b, *dwdt_aux, *mathlib_matrix, *RSP;
	float 	**dwdt;

	if (isost_model!=2 || !Te_default) return(0);

	switch (solver_type) {
	  case 'l': 
		A = alloc_matrix_dbl (Neqs, NDi+1+NDs);
		b = (double *) calloc (Neqs, sizeof(double)); 
		dwdt=alloc_matrix (Ny, Nx);
		defineLESalmostdiagonalmatrix(A, b, q, Dq, w, 1);
		solveLESalmostdiagonal(A, b, dwdt);
		for (i=0; i<Ny; i++) for (j=0; j<Nx; j++)   Dw[i][j] = dwdt[i][j]*dt;
		free_matrix_dbl (A, Neqs);
		free_matrix (dwdt, Ny);
		free(b);
		break;
	  case 'm': 
		break;
	}

	for (i=0; i<Ny; i++) for (j=0; j<Nx; j++)  w[i][j] += Dw[i][j]; 
	if (switch_topoest) {
		/*Defines the thickness of last infill Block*/
		for (i=0; i<Ny; i++)  for (j=0; j<Nx; j++)  Blocks[i_first_Block_load-1].thick[i][j] +=  Dw[i][j] /*MAX(Dw[i][j], 0)*/ ;
	}

	calculate_topo(topo);

	/*Prints DEFLECTION characteristics*/
	if (verbose_level>=1) {
		PRINT_GRID_INFO (topo, "topogr.V ", "m");
		PRINT_GRID_INFO (w, "deflect.V", "m");
	}

	return(1);
}



int Write_Ouput()
{
	//write_file_time(w, h_water);
	write_file_deflection();
	write_file_Blocks();
	write_file_cross_section();
	write_file_drainage();
	write_file_ice();
	write_file_surftransp();
	write_file_river_basins();
	write_file_lakes();
	write_file_velocity_field();
	write_file_resume();

	/*Make GMT Postscript*/
	if (switch_ps) {
		char 	command[300];
		sprintf(command, "tisc.gmt.job %s", projectname);
		if (verbose_level>=3) 
			fprintf(stdout, "\nPostscript file '%s.ps' is going to be produced with command:", projectname) ;
		if (verbose_level>=3) 
			fprintf(stdout, "\n%s\n", command) ;
		system(command);
		if (switch_dt_output) {
			/*crop by default to the border*/
			if (strlen(gif_geom)<2) sprintf(gif_geom, "-trim -background Khaki -label 'TISC software: %s' -gravity South -append", projectname);
			sprintf(command, "convert -density 200 %s.ps %s -interlace NONE  %s%03d.jpg", /*-fill \"#ffff99\" -draw \"rectangle 8,8 90,25\" -fill \"#000055\" -font helvetica -draw \"text 12,20 t_%+3.2f_My \" */
				projectname, gif_geom, projectname, n_image);
			if (verbose_level>=3)
				fprintf(stdout, "\n%s\n", command) ;
			system(command);
			n_image++;
		}
	}

	return (1);
}


