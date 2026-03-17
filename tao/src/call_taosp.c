/*
	Used for external calls to the SURFACE PROCESSES library
	(i.e., from PARAVOZ or UHURU)

	Daniel Garcia-Castellanos, 2003

The first routine init_surf_proc is to be called at the begining of the external program, 
as initialization, before the time loop.: To this routine you pass the 
parameter values you want to use. 

The second routine is the one to call at every time step of the loop.

*/

#include "tao.h"

#include "taolib.c"
#include "taoio.c"
#include "libreria.h"



int init_surf_proc_ (
	short int *ad_Nx, 			/*Number of nodes in x and y directions: column, rows index*/
	float *ad_xmin, float *ad_xmax, 	/*[m] x minimum and maximum of th domain. (xmax-xmin) domain length in x*/
	float *ad_dt_eros, 			/*[Ma] transport time increment approximed inside Dt, internal to 'surf_proc'*/
	float *ad_Kerosdif, 			/*[m2/a] Diffusive transport erosion coefficient*/
	float *ad_K_river_cap,			/*[kg/m3] Constant of river transport capacity*/
	float *ad_erodability, 			/*[m] Length scale of fluvial erosion*/
	float *ad_erodability_sed, 		/*[m] Length scale of fluvial erosion (for sediments)*/
	float *ad_l_fluv_sedim,			/*[m] Length scale of fluvial sedimentation*/
	float *ad_rain, 			/*[l/m2/a]=[mm/a] Background runoff (water going to the drainage system)*/
	float *ad_Krain, 			/*[l/m2/a/km] Proportionality of runoff with altitude*/
	float *ad_CXrain, 			/*[m] Proportionality of runoff with x (distance of duplication). 0 means no variations in x.*/
	float *ad_evaporation			/*[l/m2/a]=[mm/a] Evaporation rate at lakes*/
	)
{

	/*
	  Subroutine of definition and initialization of variables. 
	  Call this at the begining of the calling code.
	*/

	Nx=*ad_Nx, 
	x0=xmin=*ad_xmin,  xf=xmax=*ad_xmax,  
	dt_eros=*ad_dt_eros, 
	Kerosdif=*ad_Kerosdif,  K_river_cap=*ad_K_river_cap, 
	erodability=*ad_erodability,  erodability_sed=*ad_erodability_sed, 
	l_fluv_sedim=*ad_l_fluv_sedim,
	rain=*ad_rain,  Krain=*ad_Krain,  
	CXrain=*ad_CXrain, evaporation=*ad_evaporation;

	switch_verbose = 1;
	switch_debug = false;

	if (switch_verbose) fprintf(stdout, "\nRunning tAo model: \n\t%d nodes ; x=[%.1f,%.1f]",
		Nx, xmin, xmax);
	if (switch_debug) fprintf(stdout, "\n\train=%.2f l/m2/a ; Krain=%.2f l/m2/a/km",
		rain, Krain);
	if (switch_debug) fprintf(stdout, "\n\terodability=%.2f m ; erodability_sed=%.2f m ; l_fluv_sedim=%.2f m",
		erodability, erodability_sed, l_fluv_sedim);

/*	if (rain || Krain) switch_hydro = true; else switch_hydro = false;*/
	erosed_type = 1;
	switch_write_file = true;
	strcpy(projectname, "drainage");
	
	dx = (xmax-xmin) / (Nx-1);
	dt_eros *= Matosec;
	Kerosdif *= 1e6/Matosec ;
	lost_rate *= 1e-2 * 1e-3;
	rain *= 1e6/Matosec/1e3;
	Krain *= 1e6/Matosec/1e3/1e3;
	evaporation *= 1e6/Matosec/1e3;
	
	sea_level = 0;
	
	commstdout=stdout; 
	numBlocks=0;
	denssedim = 2200;
	denscrust = 2780;

	Allocate_Memory_for_external_use ();
	
	insert_new_Block(numBlocks);
	Blocks[numBlocks-1].density=denssedim;
	
	if (switch_verbose) fprintf(stdout, "\ntAo initialisation done.");
/*	if (switch_debug)   fprintf(stdout, "\tswitch_hydro=%d", switch_hydro);*/
	return 1;
}




int call_surf_proc_ (
	float *ad_dt, 			/*elapsed time [s]*/
	float *topo,			/*array of topography [m]*/
	float *sed_thick,		/*array of sediment thickness [m]*/
	short int *ad_write_files	/*<> 0 means write output files*/
	)
{
	/*
	  Sediment thickness and topography variation. 
	    topo, sed_thick: from xmin to xmax.
	    Input: topo and sed_thick, topography and Sediment thickness at t.
	    Output: topo and sed_thick, topography and Sediment thickness at t+dt.
	*/

	int i, j;
        total_lost_sed_mass=total_sed_mass=total_bedrock_eros_mass=0;

	dt = *ad_dt;
	if (switch_debug) fprintf(stdout, "\nTransporting during %.2f Ma: %d nodes ; write files = %d ",
		dt/Matosec, Nx, *ad_write_files);

	ModelConfig cfg = {0}; ModelContext ctx = {0};
	cfg.Nx = Nx; cfg.dx = dx; cfg.riverbasinwidth = riverbasinwidth; 
	cfg.denssedim = denssedim; cfg.denscrust = denscrust; cfg.erosed_model = erosed_type;
	ctx.dt = dt; ctx.dt_eros = dt_eros; ctx.topo = topo; ctx.numBlocks = numBlocks;

	for (j=0; j<Nx; j++)  {
		Blocks[numBlocks-1].thick[j]=sed_thick[j]; 
		Blocks_base[j] = topo[j] - Blocks[numBlocks-1].thick[j];
	}
	
	/*Diffusive Erosion: adds to the topo and the next load Dq and removes material from Blocks*/
	/*if (switch_debug) fprintf(stdout, "\nCalling Diffussive_Eros: ");
	Diffusive_Eros_2D (topo, Kerosdif, dt, dt_eros/5);
	*/

	/*Adds background erosion and sea sedimentation*/
	/*if (switch_debug) fprintf(stdout, "\nCalling constant_rate_eros: ");
	constant_rate_eros (topo, Keroseol, Ksedim, maxprofsed, sea_level, switch_sea, dt, dt_eros, Time, n_eros_level_input_points, var_eros_level, &eros_level);
	*/


	/*Fluvial Transport: adds to the topo and the next load Dq and removes material from Blocks*/
	if (switch_debug) fprintf(stdout, "\nCalling Fluvial_Transport: ");
	if (erosed_type>=2) Fluvial_Transport (&cfg, &ctx, topo, dt);
	numBlocks = ctx.numBlocks;

	if (*ad_write_files) {
		if (switch_debug) fprintf(stdout, "\nWritting ST & drainage file.  rain=%.2f l/m/a total_rain=%.2e m3/s", rain/1e6*Matosec*1e3, total_rain);
		write_file_erosed(&cfg, &ctx, total_erosion);
	}

	for (j=0; j<Nx; j++) {
		topo[j]=topo[j]; sed_thick[j]=Blocks[numBlocks-1].thick[j]; 
	}
	fflush (commstdout);

        fprintf(commstdout, "\n  noSed er.: %+8.2e N    sed.incr.: %+8.2e N outp.seds: %+8.2e N  ",  total_bedrock_eros_mass*g, total_sed_mass*g, total_lost_sed_mass*g);
	return 1;
}



int Allocate_Memory_for_external_use()
{
	/* Allocates dynamic memory for the arrays and initializes them to zero*/
	int i, j;
		
	if (switch_verbose) fprintf(stdout, "\ntAo memo initialisation.");

	topo		= calloc(Nx, sizeof(float));
	Blocks_base	= calloc(Nx, sizeof(float));
	Dq		= calloc(Nx, sizeof(float));
	w		= calloc(Nx, sizeof(float));

	Blocks =    	(struct BLOCK_1D *) calloc(NmaxBlocks, sizeof(struct BLOCK_1D));

	sortcell =    	calloc(Nx, sizeof(int));
	for (j=0; j<Nx; j++)  {sortcell[j]=j;}

	if (erosed_type>=2) {
		drainage =    	(struct DRAINAGE_1D *) calloc(Nx, sizeof(struct DRAINAGE_1D));
	}
	if (erosed_type) {
		eros_now	= calloc(Nx, sizeof(float));
		total_erosion	= calloc(Nx, sizeof(float));
	}

	ModelConfig cfg = {0}; ModelContext ctx = {0};
	cfg.Nx = Nx; cfg.dx = dx;

	/*Allocation for Lake #0, which is not used.*/
	Lake = (struct LAKE_INFO_1D *) calloc(1, sizeof(struct LAKE_INFO_1D));
	fflush (commstdout);

	return(1);
}
