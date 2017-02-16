/*
	Used for external calls to TISC libraries
	(i.e., from uhurutisc.f)

	Daniel Garcia-Castellanos, 2001-2012
*/

#include "tisc.h"

#include "tisclib.c"
#include "tiscio.c"
#include "../lib/libreria.h"


/*
TO CALL TISC FROM FORTRAN TRY THE NAMES OF ROUTINES HERE EITHER 
IN UPPERCASE WITH NO FINAL UNDERSCORE, OR JUST THE OPPOSITE.
*/
int init_tisc_ (
	short int *ad_Nx, short int *ad_Ny,	/*Number of nodes in x and y directions: column, rows index*/
	float *ad_xmin, float *ad_xmax, 	/*[m] x minimum and maximum of th domain. (xmax-xmin) domain length in x*/
	float *ad_ymin, float *ad_ymax, 	/*[m] y minimum and maximum of th domain. (ymax-ymin) domain length in y*/
	float *ad_dt_eros, 			/*[Ma] transport time increment approximed inside Dt, internal to 'surf_proc'*/
	float *ad_Kerosdif, 			/*[m2/a] Diffusive transport erosion coefficient*/
	float *ad_K_river_cap,			/*[kg/m3] Constant of river transport capacity*/
	float *ad_erodibility, 			/*[m] Length scale of fluvial erosion*/
	float *ad_erodibility_sed, 		/*[m] Length scale of fluvial erosion (for sediments)*/
	float *ad_l_fluv_sedim,			/*[m] Length scale of fluvial sedimentation*/
	short int  *ad_hydro_model, 		/*1 for rain proportional to elevation; 2 for orographic precipitation (wind-dependent)*/
	float *ad_rain, 			/*[l/m2/a]=[mm/a] Background runoff (water going to the drainage system)*/
	float *ad_Krain, 			/*[l/m2/a/km] Proportionality of runoff with altitude*/
	float *ad_windazimut, 			/*0 is north, then clockwise*/
	float *ad_relative_humidity, 		/*0 to 1, for incomming wind*/
	float *ad_CXrain, 			/*[m] Proportionality of runoff with x (distance of duplication). 0 means no variations in x.*/
	float *ad_CYrain, 			/*[m] Proportionality of runoff with y (distance of duplication). 0 means no variations in y.*/
	float *ad_evaporation_ct,		/*[l/m2/a]=[mm/a] Evaporation rate at lakes*/
	float *ad_Te_default			/*[m] Elastic Thickness*/
	
	)
{

	int i, j;

	/*
	  Subrutina de definicio de variables. Cridar-la al principi de l'UHURU
	*/

	Nx=*ad_Nx, Ny=*ad_Ny,
	xmin=*ad_xmin,  xmax=*ad_xmax,  ymin=*ad_ymin,  ymax=*ad_ymax, 
	dt_eros=*ad_dt_eros, 
	Kerosdif=*ad_Kerosdif,  K_river_cap=*ad_K_river_cap, 
	erodibility=*ad_erodibility,  erodibility_sed=*ad_erodibility_sed, 
	l_fluv_sedim=*ad_l_fluv_sedim,
	rain=*ad_rain, Krain=*ad_Krain,  
	CXrain=*ad_CXrain, CYrain=*ad_CYrain,  evaporation_ct=*ad_evaporation_ct;
	Px=Py=Pxy=0; K_ice_eros=0;
	Te_default=*ad_Te_default;
	Timeini=0; Time=1;

	verbose_level = 1;
	switch_topoest=NO;
	switch_write_file=1;
	lost_rate=.1;

	if (verbose_level>=2) fprintf(stdout, "\nInitializing TISC: \n\t%dx%d ; x=[%.1f,%.1f] y=[%.1f,%.1f]",
		Nx, Ny, xmin, xmax, ymin, ymax);
	if (verbose_level>=3) fprintf(stdout, "\n\terodibility=%.2f m ; erodibility_sed=%.2f m ; l_fluv_sedim=%.2f m",
		erodibility, erodibility_sed, l_fluv_sedim);

	if (rain || Krain) hydro_model = *ad_hydro_model; else hydro_model = 0;
	windazimut=*ad_windazimut;
	relative_humidity=*ad_relative_humidity;
	if (rain || Krain) erosed_model = 6; else erosed_model = 1;
	switch_write_file = YES;	switch_write_file_Units = YES;
	strcpy(projectname, "res_tisc");

	dx = (xmax-xmin) / (Nx-1);
	dy = (ymax-ymin) / (Ny-1);
	dxy = sqrt(dx*dx+dy*dy);
	dt_eros *= Matosec;
	Kerosdif *= 1e6/Matosec ;
	lost_rate *= 1e-2 * 1e-3;
	rain *= 1e6/Matosec/1e3;
	if (hydro_model < 2) Krain *= 1e6/Matosec/1e3/1e3;
	if (verbose_level>=4) fprintf(stdout, "\thydro_model=%d", hydro_model);

	evaporation_ct *= 1e6/Matosec/1e3;
	temp_sea_level += TEMP_FREEZE_WATER;

	sea_level = 0;
	
	numUnits=0;
	densenv = 0;
	denssedim = 2200;
	denscrust = 2780;
	densmantle = 3300;
	densasthen = 3250;

	Allocate_Memory_for_external_use();
	insert_new_unit(numUnits);
	Units[numUnits-1].density=denssedim;
	Units[numUnits-1].erodibility=erodibility_sed;
	
	for (i=0; i<Ny; i++) for (j=0; j<Nx; j++) D[i][j] = ET2RIG(Te_default);
	strcpy(boundary_conds, "6566");
	strcpy(eros_bound_cond, "2111");

	if (verbose_level>=3) fprintf(stdout, "\nTISC initialisation done.");	
	if (verbose_level>=3) fprintf(stdout, "\n\thydro_model=%d", hydro_model);
	if (verbose_level>=3) fprintf(stdout, "\n\train=%.2f l/m2/a ; Krain=%.2f l/m2/a/km", rain/(1e6/Matosec/1e3), Krain);
	if (verbose_level>=3) fprintf(stdout, "\n\tCXrain, CYrain = %.2e, %.2e m", CXrain, CYrain);
	fflush(stdout);
}




int call_surf_proc_ (
	float *ad_dt, 			/*elapsed time [s]*/
	float *topo_array,		/*skyline array of topography [m]*/
	float *sed_thick_array,		/*skyline array of sediment thickness [m]*/
	short int *ad_write_files	/*different of 0 to write output files*/
	)
{
	/*
	  Sediment thickness and topography variation. 
	    topo_array, sed_thick_array: vector from (ymax,xmin) to (ymin, xmax), per rows.
	    Input: topo_array and sed_thick_array, topography and Sediment thickness at t.
	    Output: topo_array and sed_thick_array, topography and Sediment thickness at t+dt.
	*/

	int i, j;
        total_lost_sed_mass=total_sed_mass=total_bedrock_eros_mass=0;

	dt = *ad_dt;
	if (verbose_level>=4) fprintf(stdout, "\nTransporting during %.2f Ma: %dx%d ; write files = %d ",
		dt/Matosec, Nx, Ny, *ad_write_files);

	for (i=0; i<Ny; i++)  for (j=0; j<Nx; j++)  {
		eros_now[i][j]=0;
		w[i][j]=0; /*to avoid changes in **topo in calculate_topo()*/
		topo[i][j]=topo_array[i*Nx+j]; 
		Units[numUnits-1].thick[i][j]=sed_thick_array[i*Nx+j]; 
		Units_base[i][j] = topo[i][j] - Units[numUnits-1].thick[i][j];
	}


	/*Diffusive Erosion: adds to the topo and the next load Dq and removes material from units*/
	if (verbose_level>=4) fprintf(stdout, "\nCalling Diffusive_Eros: ");	fflush(stdout);
	Diffusive_Eros (Kerosdif, dt, dt_eros/5);

	/*Fluvial Transport: adds to the topo and the next load Dq and removes material from units*/
	if (verbose_level>=4) fprintf(stdout, "\nCalling Fluvial_Transport: ");	fflush(stdout);
	if (hydro_model) Surface_Transport (topo, topo, dt, dt_eros, erosed_model, lake_instant_fill);
	
	if (*ad_write_files) {
		//if (verbose_level>=4) 
		fprintf(stdout, "\nWritting drainage file.  rain=%.2f l/m/a", rain/1e6*Matosec*1e3); fflush(stdout);
		write_file_drainage ();
		write_file_lakes();
		/*if (verbose_level>=4) fprintf(stdout, "\nWritting basin file.");
		write_file_river_basins ();*/
		if (verbose_level>=4) fprintf(stdout, "\nWritting st file."); fflush(stdout);
		write_file_surftransp ();
	}

	for (i=0; i<Ny; i++) for (j=0; j<Nx; j++) {
		topo_array[i*Nx+j]=topo[i][j]; sed_thick_array[i*Nx+j]=Units[numUnits-1].thick[i][j]; 
	}

        //fprintf(stdout, "\n  noSed er.: %+8.2e N    sed.incr.: %+8.2e N outp.seds: %+8.2e N  ",  total_bedrock_eros_mass*g, total_sed_mass*g, total_lost_sed_mass*g);
	fflush (stdout);
}



int call_flexure_ (
	float *load_array,		/*input:  skyline array of load [m]*/
	float *deflection_array,	/*output: skyline array of deflection [m]*/
	short int *ad_write_files	/*different of 0 to write output files*/
	)
{
	/*
	  To call flexure from a FORTRAN program 
	*/

	int i, j, NDi=2*Ny, NDs=2*Ny, Neqs=Nx*Ny;
	double **A, *b;

	if (verbose_level>=4) fprintf(stdout, "\nCalculating Flexure: %dx%d ; write files = %d ",
		Nx, Ny, *ad_write_files);

	for (i=0; i<Ny; i++)  for (j=0; j<Nx; j++)  {Dq[i][j]=load_array[i*Nx+j]; w[i][j]=0;}
	/*Diffusive Erosion: adds to the topo and the next load Dq and removes material from units*/
	b = (double *) calloc (Neqs, sizeof(double));
	A = alloc_matrix_dbl (Neqs, NDi+1+NDs);
	if (verbose_level>=4) fprintf(stdout, "\nCalling LES_Matrix");
	defineLESalmostdiagonalmatrix(A, b, q, Dq, w, 0); 
	if (verbose_level>=4) fprintf(stdout, "\nCalling solveLEScasiD");
	solveLESalmostdiagonal(A, b, w);
	free_matrix_dbl (A, Neqs);
	free(b);

	if (*ad_write_files) {
		/*if (verbose_level>=4) fprintf(stdout, "\nWritting deflection file.");
		nwrotenfiles=0;
		write_file_time(w, topo);
		if (verbose_level>=4) fprintf(stdout, "\nWritting units file.");
		write_file_Units();*/
	}

	for (i=0; i<Ny; i++) for (j=0; j<Nx; j++) {
		deflection_array[i*Nx+j]=w[i][j]; 
	}
	fflush (stdout);
}



