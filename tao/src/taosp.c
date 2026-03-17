/*
LIBRARY  FOR  SURFACE PROCESSES of tao.c

Daniel Garcia-Castellanos, 1994-2003
*/

#include "tao.h"

/*Declaration of functions at libreria.c*/
float 	**alloc_matrix  (int num_fil, int num_col);
int 	free_matrix 	(float **matrix, int num_fil);


int Surface_Transport (ModelConfig *cfg, ModelContext *ctx, float *topo) 
{
	/*
	  THIS ROUTINE COMPUTES CHANGES IN TOPOGRAPHY AND ISOSTATIC LOAD DUE
	  TO FLUVIAL EROSION/TRANSPORT/SEDIMENTATION.
	*/

	int 	i, k, iter, n_iters;	/*Number of substeps to subdivide fluvial processes*/
	float	dt_fv, 
		total_lost_sed_mass=0;

	if (!cfg->hydro_model) return(0);

	PRINT_DEBUG("Calculating surface transport");

	calculate_topo(cfg, ctx, topo);

	if (cfg->verbose_level>=3) fprintf(stdout, " ft"); fflush(stdout);

	n_iters = MAX_2(floor(ctx->dt/ctx->dt_eros+.5), 1);
	dt_fv = ctx->dt/n_iters;
	PRINT_INFO("n_iters=%3d", n_iters);

	/*Distributes fluvial erosion in 'n_iters' substeps:*/
	for (iter=0; iter<n_iters; iter++) {
	    if (cfg->verbose_level>=3) fprintf(stdout, "\b\b\b%3d", n_iters-iter); fflush(stdout);

	    ctx->total_rain=total_lost_water=total_evap_water=0;

	    for (i=0; i<cfg->Nx; i++) {
	    	drainage[i].masstr = 0;
	    	drainage[i].discharge = 0;
	    }

	    /*Resorts the matrix of topography.*/
	    ReSort_Array (topo, sortcell, cfg->Nx);

	    Define_Drainage_Net(sortcell, cfg, ctx);

	    Calculate_Precipitation_Evaporation(cfg, ctx);

	    /*read_file_node_defs(dt_st);*/

	    Calculate_Discharge(sortcell, cfg, ctx);

	    Fluvial_Transport (cfg, ctx, topo, dt_fv);
	}
	if (cfg->verbose_level>=3) fprintf(stdout, "\b\b\b"); fflush(stdout);

	/*Print relevant statistics*/
	if (cfg->verbose_level>=1) {
		float 	error;
		PRINT_ARRAY_INFO_1D(secsperyr*precipitation, "precipit.", "m/yr", "m2/yr");
		PRINT_ARRAY_INFO_1D(secsperyr*evaporation,   "evaporat.", "m/yr", "m2/yr");
		PRINT_SUMLINE("rain_now : %+8.2e m3/s  evap_wat: %+8.2e m3/s outp_water: %+8.2e m3/s", ctx->total_rain, total_evap_water, total_lost_water); 
		if (ctx->total_rain) error=-(ctx->total_rain-total_evap_water-total_lost_water)/ctx->total_rain*100; else error = -9999;
			if (fabs(error)>=1)
				PRINT_WARNING("water_balance: %.1f%% (>0 => disch>rain)", error);
	}
	if (cfg->erosed_model || cfg->verbose_level>=3) {
		float error_seds=0;
		PRINT_SUMLINE("erosion.nosd   = %+10.2e N/m  sedim.inc= %+10.2e N/m  outp.seds=  %+8.2e N/m  ",  ctx->total_bedrock_eros_mass*g/cfg->riverbasinwidth, ctx->total_sed_mass*g/cfg->riverbasinwidth, total_lost_sed_mass*g/cfg->riverbasinwidth);
		if (ctx->total_bedrock_eros_mass) error_seds= -(ctx->total_bedrock_eros_mass-ctx->total_sed_mass-total_lost_sed_mass)/ctx->total_bedrock_eros_mass*100;
		if (fabs(error_seds)>=1 && (cfg->erosed_model>1 || cfg->verbose_level>=3))
			PRINT_WARNING("seds._balance: %.1f%% (>0 => sed>eros)", error_seds);
	}
	return (1);
}



int Calculate_Discharge (int *sortcell, ModelConfig *cfg, ModelContext *ctx)
{
	/*
	  CALCULATES THE WATER DISCHARGE IN THE NETWORK.
	  Modifies the drainage information to account for evaporation.
	*/

	int 	i, k, l, il, iter, isort, 
		ind, dind ;
	float	rain_here_now, 
		dd=0;

	/*
	  This bucle starts from the top point and
	  descends transferring water and eroded mass
	*/
	for (isort=0; isort<cfg->Nx; isort++) {
	    ind  = sortcell[isort];
	    dind = drainage[ind].dr;
	    il   = drainage[ind].lake;

	    /*Check if this is an underburden lake node*/
	    if (il) {
	    	float lake_evaporation=0;
		for (i=0; i<Lake[il].n; i++) lake_evaporation += evaporation[Lake[il].cell[i]] * cfg->riverbasinwidth * cfg->dx;
	    	if (Lake_Input_Discharge(cfg, il) < lake_evaporation && Lake[il].n > 1) {
			Attempt_Delete_Node_From_Lake (cfg, ctx, ind);
			dind = drainage[ind].dr;
    			il = drainage[ind].lake;
	    	}
	    }

	    /*Calculate the distance to the output node*/
	    switch (drainage[ind].type) {
	      case 'R':
	      case 'E':
	    	    if (IN_DOMAIN_1D(dind))
		    	dd = fabs(cfg->dx*(dind-ind));
		    else 
		    	dd = 0;
	    	    break;
	      case 'L':
	    	    dd = 0;
	    	    break;
	      default:
	    	    PRINT_ERROR("node <%d> has no defined type.", ind);
	    }

	    /*Adds the rainfall water (m3/s/m) to the water transported to this cell: */
	    rain_here_now = precipitation[ind] * cfg->riverbasinwidth * cfg->dx;
	    /*Put the rain of lakes and sea in their saddles*/
	    if (drainage[ind].type == 'L') {
	    	if (Lake[il].n_sd) rain_here_now = 0;
	    }
	    if (drainage[ind].type == 'E') {
		int i;
		/*Add the rain from lake nodes draining to this saddle.*/
		for (i=0; i<Lake[il].n; i++) {
		    if (drainage[Lake[il].cell[i]].dr == ind) 
			rain_here_now += precipitation[Lake[il].cell[i]] * cfg->riverbasinwidth * cfg->dx;
		}
	    }
	    
	    ctx->total_rain += rain_here_now;
	    drainage[ind].discharge += rain_here_now;
	    if (drainage[ind].type == 'L') 
	    	if (!Lake[il].n_sd) {
		    drainage[ind].discharge -= rain_here_now;
		    total_evap_water += rain_here_now;
		}

	    /*Remove evaporated water from the lake discharge in saddles*/
	    if (drainage[ind].type == 'E') {
	    	float lake_evaporation=0, factor;
		for (i=0; i<Lake[il].n; i++) lake_evaporation += evaporation[Lake[il].cell[i]] * cfg->riverbasinwidth * cfg->dx;
	    	if (Lake_Input_Discharge(cfg, il)) factor = MIN_2(1, lake_evaporation/Lake_Input_Discharge(cfg, il));
		else factor = 0;
	    	total_evap_water        += drainage[ind].discharge * factor;
	    	drainage[ind].discharge -= drainage[ind].discharge * factor;
	    }
	    /*Transfers water.*/
	    if (IN_DOMAIN_1D(dind)) {
		/*Remove evapotranspirated water from the rivers*/
		total_evap_water            += drainage[ind].discharge * MIN_2(lost_rate*dd, 1); /* lost_rate is global */
		drainage[ind].discharge     -= drainage[ind].discharge * MIN_2(lost_rate*dd, 1); /* lost_rate is global */
		switch (drainage[dind].type) {
		    case 'L':
			/*Check: this shouldn't happen (a node transferring to a lake with a higher level).*/
			if (Lake[drainage[dind].lake].n_sd) if (IN_DOMAIN_1D(drainage[dind].dr)) if (ctx->topo[drainage[dind].dr] > ctx->topo[ind])
				if (cfg->verbose_level>=1) PRINT_ERROR("fluvial transport: node %d transferring to lake in %d is < than the saddle %d:  %.2f<%.2f.", ind, dind, drainage[dind].dr, ctx->topo[ind], ctx->topo[drainage[dind].dr]);
			/*Drain to the lake node*/
			drainage[dind].discharge += drainage[ind].discharge;
			/*Drain also to its saddle if it has*/
			if (IN_DOMAIN_1D(drainage[dind].dr)) 
				drainage[drainage[dind].dr].discharge += drainage[ind].discharge;
			else {
				/*Do nothing, border saddle transfer will be calculated below*/
			}
			break;
		    case 'R':
		    	if (dind==ind) {
				PRINT_ERROR("\aNo deberia pasar por aqui!. ind=%d, discharge=%f", ind, drainage[ind].discharge);
				total_evap_water += drainage[ind].discharge; 
			}
			else 	
				drainage[dind].discharge += drainage[ind].discharge;
			break;
		    case 'E':
		    	/*Lake internal drainage is done above*/
			if (il != drainage[dind].lake) 
				drainage[dind].discharge += drainage[ind].discharge;
			break;
		    default:
		    	PRINT_ERROR("Fluvial Transport: node <%d> draining to <%d> is type missing.", ind, dind);
			break;
		}
	    }
	    else {
		    if (AT_BORDER_1D(ind)) {
			/*Transfers out of model.*/
			total_lost_water += drainage[ind].discharge;
		    }
		    else {
		    	/*Evaporates water from underburden lake nodes*/
			total_evap_water += drainage[ind].discharge;
		    }
	    }
	}
	
	/*Calculate lake altitudes and volumes*/
	for (il=1; il<=ctx->nlakes; il++) {
		if (Lake[il].n_sd) {
			if (ctx->topo[Lake[il].sd[0]] < ctx->sea_level && AT_BORDER_1D(Lake[il].sd[0])) {
				Lake[il].alt = ctx->sea_level;
			}
			else
				Lake[il].alt = ctx->topo[Lake[il].sd[0]];
			}
		else {
			Lake[il].alt = ctx->topo[Lake[il].cell[Lake[il].n-1]];
		}
		/*Check: altitude should be the same as the one of the last node (except in the sea)*/
		if (Lake[il].alt != ctx->topo[Lake[il].cell[Lake[il].n-1]]) {
		    if (Lake[il].n_sd) if (ctx->topo[Lake[il].sd[0]] >= ctx->sea_level || !AT_BORDER_1D(Lake[il].sd[0])) 
			PRINT_WARNING("lake %d (overburden but not sea) should have altitude %.1f m instead of %.1f m.", il, Lake[il].alt, ctx->topo[Lake[il].cell[Lake[il].n-1]]);
		}
		Lake[il].vol = 0;
		for (i=0; i<Lake[il].n; i++)  Lake[il].vol +=  (Lake[il].alt - ctx->topo[Lake[il].cell[i]]);
		Lake[il].vol *= cfg->dx;
	}

	return (1);
}



int constant_rate_eros (ModelConfig *cfg, ModelContext *ctx, float Keroseol, float Ksedim, int water_load, int n_eros_level_input_points, float **var_eros_level, float *eros_level) /*returns the eros_level*/
{
	int i;
	float Dh;

	/*
	  Adds background erosion and sea sedimentation
	*/

	if (!cfg->erosed_model) return (0);

	/*Finds erosion level*/
	if (n_eros_level_input_points) {
 	    for (i=0; i<n_eros_level_input_points; i++) {
 		if (var_eros_level[i][0]>=ctx->Time) break;
 	    }
 	    if (i!=0 && i!=n_eros_level_input_points) {
		*eros_level =	 ( (ctx->Time-var_eros_level[i-1][0])*var_eros_level[i][1] +
					  (var_eros_level[i][0]-ctx->Time)*var_eros_level[i-1][1] )
				     / (var_eros_level[i][0]-var_eros_level[i-1][0]);
	    }
	    else {
		if (i==0)			  *eros_level = var_eros_level[0][1];
		if (i==n_eros_level_input_points) *eros_level = var_eros_level[n_eros_level_input_points-1][1];
	    }
	}
	else *eros_level = 0;

	/*Calculate eros/sed*/
	for (i=0; i<cfg->Nx; i++) {
	    /*EROSION proportional to elevation*/
	    if (ctx->topo[i] > ctx->sea_level+*eros_level) {
		    Dh = Keroseol * ctx->dt * (ctx->topo[i]-ctx->sea_level-*eros_level);
		    Erode (cfg, ctx, THICK2SEDMASS_1D(cfg, Dh), i);
	    }
	    /*SEDIMENTATION*/
	    else if (water_load) {
		    Dh = MIN_2(Ksedim*ctx->dt, ctx->sea_level+*eros_level-ctx->topo[i]);
	    	    Sediment (cfg, ctx, THICK2SEDMASS_1D(cfg, Dh), i);
	    }
	}
	return 1;

}



int Define_Drainage_Net (int *sortcell, ModelConfig *cfg, ModelContext *ctx) 
{
	/*
	Here the drainage net is defined by classifying the nodes in the domain 
	in 3 types:
	'R' river;  'L' lake/sea;  'E' lake/sea saddle/exit;

	This routine will define the lakes (these include topographic minima,
	flats,  i.e., lakes themself, sea and plains) at every time
	substep  of the fluvial erosion. The topography is asumed to be
	perfectly  sorted in 'sortcell'. 
	The seas can be distinguished from proper lakes only because they are 
	the only lakes having saddles in the boundary and below sea level.
	Lakes occupy topographic minima and (before evaporation is calculated)
	always have a saddle where water exits. Sea is defined as one or more
	lakes including at least one node below sea level in the boundary. 
	This connection with the boundary is an exit/outlet. These outlets 
	below sea level, allow to differenciate sea-lakes from the rest of 
	normal lakes. Note that if a lake has no connection with the boundary 
	then it will either be a closed lake with no outlet or it will fill 
	above sea-level, looking for an outlet. 

	The next algorithm is followed for every node, starting in the lowest
	and ascending node by node. Initially nothing is known about the
	nodes,  therefore in the algorithm nothing is known about the nodes
	above the  current one. 
	I call a lake is 'defined' when a saddle for the overflowing has 
	already been found and so, the limits of the lake are well known.  I
	call an 'adjacent' to any of the 8 possible surrounding nodes of the 
	current node.
	I call a 'saddle' to any point of the lake that transfers out of the 
	lake. This can be either a saddle itself of the topography or a node 
	in the border, or the limit of a plain that has a lower adjacent.

	BEFORE (for the nodes with topo <= sea_level): these nodes are treated 
		as normal lake nodes. This nodes can form different lakes 
		which, in case they are closed, use a saddle.
		
	LOOP (for every node with a topo > sea_level):
	1. if has no undefined-lake adjacents
	    1.1 if < than all the adjacents and it is not in the border, 
	    	then mark it as a new undefined lake.
	    1.2 if = than at least one of the adjacents
	    	1.2.1 if > than one adjacent, then mark the node as a new 
	    		undefined lake and mark it as the saddle.
		1.2.2 if <= than all adjacents, then mark the node as a new 
			undefined lake.
	2. if has undefined lake adjacents 
	    2.1. if <= than all the adjacents that are not undefined 
		lakes, then mark it as the adjacent undefined lake.
	    2.2. if > than at least one adjacent that is not undefined 
	    	lake, then mark it as the adjacent undefined lake, mark it as 
	    	a saddle and mark the minimum adjacent that is not undefined 
	    	lake as the node to transfer. 
	3. if is a lake labelled different than an undefined lake 
	  	adjacent, then unify both lakes adding the saddles. If one of 
	  	the lakes is defined then the resulting lake is too.
	4. if is a lake in the border and it has not been defined as a saddle, 
		then mark it as a saddle.
	5. if the next node in the sortcell array is > than the present 
		node, then mark the lakes which saddles are at the present 
		height as defined.

	AFTER (checks): each lake has one or more saddles (all of them with 
		the same height) with a node next to them to which to transfer 
		water and seds. Lakes cannot be formed only with saddles.

	EXAMPLE in 1D assuming no evaporation:
         22  1111  333 444444555 666666  66   <--lake_number
	 ELRRSSSSRRELLRELLLLLLLLRLLLLLE  LERR <--node_type
      z	|                       #      |     |
	|              #        ##### #| ##  |
	|          #  ###   #  ########| ### |
	|  ##     ### #### ############| ####|
	|####    ######################| ####|
    z=0	|----#-------------------------| ----|
	|     ###                      |     |
	|       #                      |     |
         12120000--2212222122122-122212  22--  \
         2121      2111211111211 211111  11    | <-- applied rule
         2 1                 2   2             /

	The information about the lakes is stored in two structures: 'drainage' 
	(for every node in the grid) and 'Lake' (for every lake).
	*/


	int	i, isort, k, l, m, n, 
		total_lake_nodes;

	/*Delete all lakes*/
	for (i=ctx->nlakes; i>0; i--) Deallocate_Lake (ctx, i);
	for (i=0; i<cfg->Nx; i++) 
		{drainage[i].lake = 0;   drainage[i].type = '-';}

	/*Create lake #0 wich will contain nothing*/
	Lake = calloc(1, sizeof(struct LAKE_INFO));
	
	
	/*Define lakes ascending in the topo grid.*/
	for (isort=cfg->Nx-1; isort>=0; isort--) {
		int undef_lake_adj=0, n_nonposders=0, n_zeroders=0, n_negders=0, n_posders=0;
		int imaxderneg=SIGNAL, imaxderneg_noundef=SIGNAL, n_negders_nonundef_lake=0; 
		float deriv=SIGNAL,   maxderneg=0, maxderneg_noundef=0;
		bool switch_change_in_next_height=false, switch_saddle=false;

		i = sortcell[isort]; 

		/*Calculate derivates and look for undefined lakes in all directions: */
		if (i>0) {
			deriv = (ctx->topo[i-1]-ctx->topo[i])/cfg->dx;	/*left side*/
			if (deriv<0)   {
				n_negders++; 
				if (deriv<maxderneg) {maxderneg=deriv; imaxderneg=i-1;}
				if (drainage[i-1].lake >= 0) n_negders_nonundef_lake++;
			}
			if (drainage[i-1].lake<0) /*thus, if undefined*/
				undef_lake_adj = drainage[i-1].lake;
			else
				if (deriv<maxderneg_noundef) {maxderneg_noundef=deriv; imaxderneg_noundef=i-1;}
			if (deriv>0)  	n_posders++; 
			if (deriv<=0) 	n_nonposders++;
			if (deriv==0) 	n_zeroders++;
		}
		if (i<cfg->Nx-1) {
			deriv = (ctx->topo[i+1]-ctx->topo[i])/cfg->dx;	/*right side*/
			if (deriv<0)  {
				n_negders++; 
				if (deriv<maxderneg) {maxderneg=deriv; imaxderneg=i+1;}
				if (drainage[i+1].lake >= 0) n_negders_nonundef_lake++;
			}
			if (drainage[i+1].lake<0) /*thus, if undefined*/
				undef_lake_adj = drainage[i+1].lake;
			else
				if (deriv<maxderneg_noundef) {maxderneg_noundef=deriv; imaxderneg_noundef=i+1;}
			if (deriv>0)  	n_posders++; 
			if (deriv<=0) 	n_nonposders++;
			if (deriv==0) 	n_zeroders++;
		}


		/*START IDENTIFYING THE NEW NODE.*/
		/*Sea lakes:*/
		if (ctx->topo[i]<=ctx->sea_level) {
			if (undef_lake_adj) {
				Add_Node_To_Lake(i, undef_lake_adj);
			}
			else {
				int i_lake;
				Add_Node_To_Lake(i, i_lake=New_Lake(cfg, ctx));
			}
		}
		/*Normal (non-sea) lakes:*/
		else {
			if (undef_lake_adj) {
				/*Add saddle to the undefined lake.*/
				if (n_negders_nonundef_lake) {
					Add_Node_To_Lake(i, undef_lake_adj);
					Add_Saddle_To_Lake(i, imaxderneg_noundef, undef_lake_adj);
					switch_saddle=true;
				}
				/*Add node to the undefined lake.*/
				else {
					Add_Node_To_Lake(i, undef_lake_adj);
				}
			}
			else {
				/*Create a new lake if there is no way down/flat/out.*/
				if (!n_nonposders && NOT_AT_BORDER_1D(i)) {
					Add_Node_To_Lake(i, New_Lake(cfg, ctx));
				}
				/*Create a new lake if there is a flat.*/
				if (n_zeroders) {
					if (n_negders) {
						int i_lake;
						Add_Node_To_Lake(i, i_lake=New_Lake(cfg, ctx));
						Add_Saddle_To_Lake(i, imaxderneg, i_lake);
						switch_saddle=true;
					}
					else {
						Add_Node_To_Lake(i, New_Lake(cfg, ctx));
					}
				}
			}
		}
		/*Unify connected undefined lakes*/
		if (drainage[i].lake<0 && undef_lake_adj) {
			if (i>0)    if (drainage[i-1].lake<0 && drainage[i-1].lake != drainage[i].lake) {
				Unify_Lakes(cfg, ctx, drainage[i].lake, drainage[i-1].lake);
			}
			if (i<cfg->Nx-1) if (drainage[i+1].lake<0 && drainage[i+1].lake != drainage[i].lake) {
				Unify_Lakes(cfg, ctx, drainage[i].lake, drainage[i+1].lake);
			}
		}
		/*If it's a lake node in the border, then is a saddle*/
		if (drainage[i].lake && AT_BORDER_1D(i) && !switch_saddle) {
			if (eros_bound_cond[0] != 'c') 
				Add_Saddle_To_Lake(i, SIGNAL, drainage[i].lake);
		}
		if (ctx->topo[i]>ctx->sea_level) {
			/*
			  If the height is going to change in the next node or
			  this is the last node, then mark as defined all the
			  overburden lakes that have the present height.
			*/
			if (isort>0)  if (ctx->topo[i] != ctx->topo[sortcell[isort-1]]) 
				switch_change_in_next_height=true;
			if (switch_change_in_next_height || isort==0) {
				for (l=1; l<=ctx->nlakes; l++) {
					if (Lake[l].n_sd)
						if (ctx->topo[Lake[l].sd[0]] == ctx->topo[i]) {
							Define_Lake(l);
						}
				}
			}
		}
		else {
			/*
			  If this isort is the last under-sea-level node, or simply the
			  last (upper most) node, then define all 'sea-like' lakes.
			*/
			if (isort>0) {
			    if (ctx->topo[sortcell[isort-1]] > ctx->sea_level) {
				for (l=1; l<=ctx->nlakes; l++) {
					if (Lake[l].n_sd) {
						Define_Lake(l);
					}
				}
			    }
			}
			else {
				for (l=1; l<=ctx->nlakes; l++) {
					if (Lake[l].n_sd) {
						Define_Lake(l);
					}
				}
			}
		}
		/*
		  Determines drainage of the non-lake (river) nodes.
		  Add transferring and other information to 'drainage'.
		*/
		if (!drainage[i].lake) {
			drainage[i].type = 'R';
			if (imaxderneg != SIGNAL) {
			    	/*Drain to the lowest neighbour*/
				drainage[i].dr = imaxderneg;
			}
			else {
				drainage[i].dr = SIGNAL;
			}
		}
	}

	/*
	  Delete all lakes which nodes are all of them saddles. 
	  The saddle drainage is not mantained because in the borders does not work.
	*/
	for (l=1; l<=ctx->nlakes; l++) {
	    if (Lake[l].n == Lake[l].n_sd) {
		if (Lake[l].n > 2) PRINT_ERROR("'Lake' %d should consist of maximum 2 saddles.", l, Lake[l].n);
		/*
		  Determines drainage of the new non-lake (river) nodes.
		  Add transferring and other information to 'drainage'.
		*/
		drainage[Lake[l].cell[0]].lake = 0;
		drainage[Lake[l].cell[0]].type = 'R';
		if (Lake[l].n == 2) {
			drainage[Lake[l].cell[1]].lake = 0;
			drainage[Lake[l].cell[1]].type = 'R';
		}
	    	Deallocate_Lake(ctx, l);
	    	l--;
	    }
	}

	/*
	  Determines in which case are the lake nodes: 
	  	'L' lake;  'E' saddle/exit of lake;
	  Add transferring and other information to 'drainage'.
	*/
	for (l=1; l<=ctx->nlakes; l++) {
		float 	dist, mindist;
		int 	imindist=-1, ls;
		/*Saddles*/
		for (m=0; m<Lake[l].n_sd; m++) {
			ls = Lake[l].sd[m];
			drainage[ls].type = 'E';
			/*All saddles should be already 'defined' as lake node in 'drainage'*/
			if (drainage[ls].lake != l) 
				PRINT_ERROR("'Lake' %d and 'drainage' %d do not match in saddle [%d].", l, drainage[ls].lake, ls);
			if (drainage[ls].lake < 0) 
				PRINT_ERROR("saddle [%d] was supposed to be a defined lake node.", ls);
		}
		/*Lakes*/
		for (m=0; m<Lake[l].n; m++) {
			ls = Lake[l].cell[m];
			/*All lakes should be 'defined' in 'drainage'*/
			if (drainage[ls].lake < 0) {
				if (Lake[l].n_sd) PRINT_ERROR("Lake %d (%d nodes) needed to be defined late.", l, Lake[l].n);
			}
			if (drainage[ls].lake != l) {
				PRINT_ERROR("''Lake' %d (%dth of %d) and 'drainage' %d do not match in node [%d].", l, m, Lake[l].n, drainage[ls].lake, ls);
				for (i=0; i<cfg->Nx; i++) fprintf (stderr, "\n>>>%d\t%f", i, ctx->topo[i]);
				//if (drainage[ls].lake<0) drainage[ls].lake = -drainage[ls].lake;
			}
			/*Lake nodes drain to the nearest saddle:*/
			for (n=0,mindist=1e24; n<Lake[l].n_sd; n++) {
				dist = fabs(cfg->dx*(Lake[l].sd[n]-ls));
				if (dist<mindist) {mindist=dist; imindist=n;}
			}
			if (Lake[l].n_sd && drainage[ls].type != 'E') {
				drainage[ls].dr = Lake[l].sd[imindist];
				drainage[ls].type = 'L';
				if (!drainage[ls].lake) PRINT_ERROR("node [%d] should have an assigned lake.", ls);
			}
			if (!Lake[l].n_sd) 	PRINT_ERROR("Lake without exit before evaporation");
		}
	}

	/*For the same altitude, put the saddles first in sortcell*/
	for (isort=0; isort<cfg->Nx-1; isort++) {
	    int i = sortcell[isort];
	    if (drainage[i].type != 'E') {
		float topoisort; int j;
		topoisort = ctx->topo[i];
		for (j=isort+1; j<cfg->Nx; j++) {
		    float topoj;
		    topoj = ctx->topo[sortcell[j]];
		    if (topoj == topoisort) {
			if (drainage[sortcell[j]].type == 'E') {
		    	    int aux;
			    aux=sortcell[j];
			    sortcell[j]=i;
			    sortcell[isort]=aux;
			    break;
		    	}
		    }
		    else break;
		}
	    }
	}

	/*CHECK RESULTS*/
	/*Lakes 'defined' in 'drainage' should be as many as the total lake nodes.*/
	for (l=1,total_lake_nodes=0; l<=ctx->nlakes; l++) total_lake_nodes += Lake[l].n;
	for (i=k=0; i<cfg->Nx; i++) if (drainage[i].lake) k++;
	if (k != total_lake_nodes) 
		PRINT_ERROR("%d lake nodes were expected rather than %d.", total_lake_nodes, k);

	/*All saddles of a lake have the same height*/
	for (l=1; l<=ctx->nlakes; l++) {
	    float  sd_height;
	    for (m=1; m<Lake[l].n_sd; m++) {
		sd_height = ctx->topo[Lake[l].sd[m-1]];
		if (sd_height>ctx->sea_level && sd_height!=ctx->topo[Lake[l].sd[m]])
		    PRINT_ERROR("lake %d has saddles %d and %d with different height.",
		    	    l, m-1, m
		    );
	    }
	}

	/*Drainage direction and lake transferring info should match.*/
	/*Drainage type and lake info should match.*/
	for (l=1; l<=ctx->nlakes; l++) {
		for (m=0; m<Lake[l].n_sd; m++) {
			int i=Lake[l].sd[m];
			if (drainage[i].type != 'E') 
				PRINT_ERROR("node [%d] should have type 'saddle'.", i);
		}
	}

	/*All nodes transfer either to SIGNAL,SIGNAL or to a real node.*/
	for (i=0; i<cfg->Nx; i++) {
		if (drainage[i].dr < 0 && drainage[i].dr != SIGNAL)
			PRINT_ERROR("defining drainage network:  %d -->>  %d \t'%c'", i, drainage[i].dr, drainage[i].type);
	}

	
/*	for (isort=0; isort<Nx-1; isort++) {
		fprintf(stderr, "\n##Node [%d] %.1fm  %c", 
			sortcell[isort], topo[sortcell[isort]], drainage[sortcell[isort]].type);
	}
	for (i=0; i<Nx; i++) 
		fprintf(stderr, "\n**Node [%d]->[%d] %.1fm %.1fm3/s Lake %d: n=%d;  type '%c' ; inputw=%f", 
			i, drainage[i].dr, drainage[i].discharge, l=drainage[i].lake, Lake[l].n, drainage[i].type, Lake_Input_Discharge(l)); 
*/
	return(1);
}



int Diffusive_Eros_1D (ModelConfig *cfg, ModelContext *ctx, float Kerosdif) 
{
	int i, n_iters, conteros;
	float *Dheros;

	if (!cfg->erosed_model || !Kerosdif) return (0);
	Dheros = (float *) calloc(cfg->Nx, sizeof(float));

	/*Divides erosion into substeps to transmit the effect more far away*/
	n_iters = MAX_2(floor(ctx->dt/ctx->dt_eros+.5), 1);
	PRINT_INFO("n_iters=%3d", n_iters);

	for (conteros=0; conteros<n_iters; conteros ++) {
	    float   topo2deriv, mean_height, Dheros_dif_max, Dheros_dif, dteros=ctx->dt/n_iters ;
	    for (i=1; i<cfg->Nx-1; i++) {
		    /*Second derivate of topographic heigth*/
		    topo2deriv = ( ctx->topo[i-1]+Dheros[i-1] + ctx->topo[i+1]+Dheros[i+1] - 2*(ctx->topo[i]+Dheros[i]) ) / (cfg->dx*cfg->dx) ;
		    /*Mean height between points*/
		    mean_height = ( ctx->topo[i-1]+Dheros[i-1] + ctx->topo[i+1]+Dheros[i+1] ) / 2 ;
		    /*Maximum eros/sed is that which makes flat the topography*/
		    Dheros_dif_max = mean_height - (ctx->topo[i]+Dheros[i]) ;
		    /*Calculates diffusive erosion proportional to second derivate*/
		    Dheros_dif = Kerosdif * topo2deriv * dteros ;
		    if (Dheros_dif == 0) Dheros_dif_max =0;
		    if (Dheros_dif > 0) {
			    Dheros_dif = MIN_2(Dheros_dif, Dheros_dif_max);
		    }
		    else {
			    Dheros_dif = MAX_2(Dheros_dif, Dheros_dif_max);
		    }
		    /*Adds it to the total erosion during this time interval*/
		    Dheros[i] += Dheros_dif;
	    }
	}
	/*Adds results to the height and the next load Dq and removes material from the Blocks*/
	for (i=0; i<cfg->Nx; i++) {
		ctx->topo[i] += Dheros[i]; 
		if (Dheros[i]>0) {
			/*Sediment*/
			Sediment (cfg, ctx, THICK2SEDMASS_1D(cfg, Dheros[i]), i);
		}
		else {
			Erode (cfg, ctx, -THICK2SEDMASS_1D(cfg, Dheros[i]), i);
		}
	}
	free(Dheros);
	return 1;
}



int Fluvial_Transport (ModelConfig *cfg, ModelContext *ctx, float *topo, float dt_fv) 
{
	/*
	  THIS ROUTINE COMPUTES THE TOPOGRAPHY AND LOAD VARIATIONS DUE
	  TO FLUVIAL EROSION/TRANSPORT/SEDIMENTATION.
	  The model assumes Hack's Law and that the large scale 
	  erosion/sedimentation in a point depends on the water volume 
	  crossing the cell and the slope in that point.
	  >Input water comes only from the rain falling at each
	  cell. Water is transferred to the lowest sorrounding cell.
	  Water losts occur when streams reach the model boundary, and
	  when they reach the sea.
	  >Eroded material is transported with water. Mass is preserved.
	  Mass losts only occur through the domain boundaries.
	  >The criteria used from erosion/sedimentation is as follows:
	  Continent:
	  	  Sediments capacity proportional to slope and discharge.
	  	  eroded/sedimented mass proportional to
	  		  (capacity-tranported_sediment)*distance
	  	  Lake: Sedimentation of all the available material
	  		  until lake overflow.
	  Sea:
	  	  Sedimentation as in continent with discharge=0
	  	  Bathymetric minima as in continental lakes.
	*/
	int 	i, isort, k, ind, dind, n_limited_sed=0, n_limited_eros=0;
	float	total_lost_sed_mass=0, 
		minsorr_trib, minsorr, maxsorr=SIGNAL, 
		main_tribut_slope;

	    /*
	      This bucle starts from the top point and 
	      descends transferring the eroded mass
	    */
	    if (cfg->erosed_model>=2) for (isort=0; isort<cfg->Nx; isort++) {
		float d_mass; 	/*Increment of suspended mass in this cell (positive means erosion).*/
	    	ind = sortcell[isort];
		dind = drainage[ind].dr;
		d_mass = 0;
		/*Calculate max and min height in the 2 sorrounding nodes:*/
		maxsorr=minsorr=topo[ind];
		main_tribut_slope=0;
		if (IN_DOMAIN_1D(ind-1)) {
			minsorr = MIN_2(minsorr, topo[ind-1]);
			maxsorr = MAX_2(maxsorr, topo[ind-1]);
			if (topo[ind-1]>topo[ind] && drainage[ind-1].dr==ind) main_tribut_slope = MAX_2(main_tribut_slope, (topo[ind-1]-topo[ind])/cfg->dx);
		}
		if (IN_DOMAIN_1D(ind+1)) {
			minsorr = MIN_2(minsorr, topo[ind+1]);
			maxsorr = MAX_2(maxsorr, topo[ind+1]);
			if (topo[ind+1]>topo[ind] && drainage[ind+1].dr==ind) main_tribut_slope = MAX_2(main_tribut_slope, (topo[ind+1]-topo[ind])/cfg->dx);
		}

		minsorr_trib = 99999; /*very high*/
		if (IN_DOMAIN_1D(ind-1)) if (topo[ind-1]<minsorr_trib && drainage[ind-1].dr==ind) minsorr_trib = topo[ind-1];
		if (IN_DOMAIN_1D(ind+1)) if (topo[ind+1]<minsorr_trib && drainage[ind+1].dr==ind) minsorr_trib = topo[ind+1];

		/*
		  Calculates the potential d_mass increment (kg, >0 for erosion) 
		  in this cell's solid-earth column due to 
		  EROSION or SEDIMENTATION in this cell: 
		*/
		switch (drainage[ind].type) {
		    float dist, slope, transp_capacity_eq;
		    case 'L':
			/*Sedimentation will be done in Lake_Fill().*/
			d_mass  = 0;
			break;
		    case 'R':
		    case 'E':
			/*Calculate the slope to the output node*/
			if (IN_DOMAIN_1D(dind)) {
				int ild = drainage[dind].lake;
			  	dist = fabs(cfg->dx*(dind-ind));
				if (ild) slope = - (Lake[ild].alt - topo[ind]) / dist;
				else 	 slope = - (topo[dind]    - topo[ind]) / dist;
			}
			else {
				dist = cfg->dx;
				slope = main_tribut_slope;
			}

			switch (cfg->erosed_model) {
#define ERODED_ERODIBILITY   /*Takes a mean erodibility*/ float depth2average=10., dh, weight, totalweight=0, basedepth=0, erodibility_aux=0;\
				for (i=ctx->numBlocks-1; i>=0; i--) {\
					basedepth+=Blocks[i].thick[ind];\
					basedepth=MIN_2(basedepth,depth2average+.1);\
					weight=Blocks[i].thick[ind]/(basedepth+1.); totalweight+=weight;\
					erodibility_aux+=weight*Blocks[i].erodibility; \
					/*PRINT_ERROR("xxx %.2e %.2e %.2e %.2e", Blocks[i].erodibility, erodibility_aux, basedepth, weight);*/\
					if (basedepth>=depth2average) break;\
				};\
				if (basedepth<=depth2average) {\
					weight=(depth2average-basedepth)/(depth2average+1); totalweight+=weight;\
					erodibility_aux+=weight*erodibility;\
					/*PRINT_ERROR("XXX %.2e %.2e %.2e %.2e", erodibility, erodibility_aux, basedepth, weight);*/\
				}\
				if (totalweight) erodibility_aux/=totalweight; if (!erodibility_aux) erodibility_aux=erodibility; 
#define TRANSPORT_BOUNDARY_CONDITIONS \
			    if (AT_BORDER_1D(ind)) {\
				switch (eros_bound_cond[0]) {\
				    case '1':\
				    case '2':\
				    case 'c':					break;\
				    case '0':	transp_capacity_eq = 0;		break;\
				    case '3':	transp_capacity_eq *= .995;	break;\
			    }}
			  case 2:
			    /*Beaumont et al. (1992) stream power:*/
			    /*Transport capacity in equilibrium  [kg/s/m].*/
			    transp_capacity_eq = K_river_cap * drainage[ind].discharge * slope;	    /*Eq. 16 of Tucker&Slingerland, 1996*/
			    TRANSPORT_BOUNDARY_CONDITIONS;
			    /*EROSION*/
			    if (transp_capacity_eq >= drainage[ind].masstr) {
			    	ERODED_ERODIBILITY;
				d_mass  =  dist /  erodibility_aux * (transp_capacity_eq - drainage[ind].masstr) * dt_fv;
			    }
			    /*SEDIMENTATION*/
			    else {
				d_mass  =  dist / l_fluv_sedim   *  (transp_capacity_eq - drainage[ind].masstr) * dt_fv;
			    }
			    break;
			  case 3:
			    /*Tucker & Slingerland (1996) stream power:*/
			    transp_capacity_eq = K_river_cap * drainage[ind].discharge * slope;	    /*Eq. 16 of Tucker&Slingerland, 1996*/
			    TRANSPORT_BOUNDARY_CONDITIONS;
			    if (transp_capacity_eq > drainage[ind].masstr) {
				spl_m = 1/3;
				spl_n = 2/3;
				ERODED_ERODIBILITY;
	    			/*bedrock channel incision*/
				dh = erodibility_aux		/*Eq. 11 of T&S*/
					* pow((double)drainage[ind].discharge, (double)spl_m)
					* pow((double)slope,                   (double)spl_n)
					* dt_fv;
				d_mass = dh * cfg->denscrust * cfg->dx;
			    }
			    else{
				/*alluvial channel aggradation: sediment the excess*/
				d_mass = 				/*Eqs. 18 & 10 of T&S*/
					(transp_capacity_eq - drainage[ind].masstr)
					* dt_fv;
			    }
			    break;
			  case 4:
			    /*Critical slope criteria:*/
			    /*if (slope>criticslope) {d_mass  =   Keroswater * dischg[ind] * dist * (slope-criticslope);}
			      if (slope<criticslope) {d_mass  = - MIN_2(Ksedimwater * dist * (criticder-slope)/criticslope, 1) * masstr[ind];}
			    */
			    break;
			  case 6:
			    /*Garcia-Castellanos & Villasenor (2011, Nature) approach:*/
			    transp_capacity_eq = K_river_cap * drainage[ind].discharge * slope;	    /*Eq. 16 of Tucker&Slingerland, 1996*/
			    TRANSPORT_BOUNDARY_CONDITIONS;
			    if (transp_capacity_eq >= drainage[ind].masstr) {
			  	float a=1.5, Kw=1.1, aw=0.5;
				spl_m = 3*a*(1-aw)/5;
				spl_n = 7*a/10;
			  	ERODED_ERODIBILITY;
			  	/*bedrock channel incision*/
	    		    	dh = erodibility_aux/secsperyr * pow(1020*g, a)
	    			    * pow((double).05/Kw, (double)3*a/5)
				    * pow((double)drainage[ind].discharge, (double)spl_m) 
				    * pow((double)slope,(double)spl_n)
			  	    * dt_fv;
//PRINT_ERROR("??[%d] x=%.1f topo=%.1f Qw=%.2e Qeq=%.2e Qs=%.2e dh=%.2e S=%.2e erod=%.2e", ind, (x0+ind*dx)/1e3, topo[ind], drainage[ind].discharge, transp_capacity_eq, drainage[ind].masstr, dh, slope, erodibility_aux);
				if (transp_capacity_eq) dh *= (transp_capacity_eq-drainage[ind].masstr)/transp_capacity_eq;
			  	d_mass = THICK2SEDMASS_1D(cfg, dh);
//PRINT_ERROR("##[%d] x=%.1f topo=%.1f Qw=%.2e Qeq=%.2e Qs=%.2e dh=%.2e dmass=%.2e", ind, (x0+ind*dx)/1e3, topo[ind], drainage[ind].discharge, transp_capacity_eq, drainage[ind].masstr, dh, d_mass);
//PRINT_ERROR("### %.1f %.2e %c %.2fm %.2e %.2e", ind*dx+x0, d_mass, drainage[ind].type, dh, transp_capacity_eq, drainage[ind].masstr);
			    }
			    else{
			  	/*alluvial channel aggradation: sediment the excess*/
			  	d_mass =
			  		dist / l_fluv_sedim * (transp_capacity_eq - drainage[ind].masstr)
			  		* dt_fv;
//PRINT_ERROR(">>[%d] x=%.1f topo=%.1f Qw=%.2e Qeq=%.2e Qs=%.2e dmass=%.2e", ind, (x0+ind*dx)/1e3, topo[ind], drainage[ind].discharge, transp_capacity_eq, drainage[ind].masstr, d_mass);
//PRINT_ERROR("))) %.1f %.2e %c %.6e %.6e  %.2fm", ind*dx+x0, d_mass, drainage[ind].type, transp_capacity_eq, drainage[ind].masstr, dist);
			    }
			    break;
			  case 7:
			    /*Ferrier et al. (2013, Nature) unit stream power approach:*/
			    transp_capacity_eq = K_river_cap * drainage[ind].discharge * slope;	    /*Eq. 16 of Tucker&Slingerland, 1996*/
			    TRANSPORT_BOUNDARY_CONDITIONS;
			    if (transp_capacity_eq >= drainage[ind].masstr) {
			  	float a=1, Kw=1.1, aw=0.5;
				spl_m = (1-aw);
				spl_n = a;
			  	ERODED_ERODIBILITY;
			  	/*bedrock channel incision*/
			  	dh = erodibility_aux/secsperyr * pow(1020*g, a) / Kw 
			  		* pow((double)drainage[ind].discharge, (double)spl_m) 
					* pow((double)slope,  		       (double)spl_n)
			  		* dt_fv;
				if (transp_capacity_eq) dh *= (transp_capacity_eq-drainage[ind].masstr)/transp_capacity_eq;
			  	d_mass = THICK2SEDMASS_1D(cfg, dh);
			    }
			    else{
			  	/*alluvial channel aggradation: sediment the excess*/
			  	d_mass =
			  		dist / l_fluv_sedim * (transp_capacity_eq - drainage[ind].masstr)
			  		* dt_fv;
			    }
			    break;
			}
			break;
		    default:
		  	PRINT_ERROR("Fluvial Transport: node <%d> has no defined type.", ind);
		}
//PRINT_ERROR(">>> %.1f %.2e %c", ind*dx+x0, d_mass,drainage[ind].type );

		if (AT_BORDER_1D(ind)) if (eros_bound_cond[0]=='2') d_mass = 0;

		/*Limit d_mass with the sorrounding topo*/
		switch (drainage[ind].type) {
		  float Dhsed;
		  case 'L':
			/*Sedimentation will be done in Lake_Fill().*/
			d_mass = 0;
			break;
		  case 'R':
			/*RIVER EROSION/SEDIMENTATION*/
			Dhsed = -MASS2SEDTHICK_1D(cfg, d_mass);
			/*Sedimentation*/
			if (Dhsed>0) {
				float Dhsedmax=MAX_2(minsorr_trib-topo[ind]-1e-2, 0);
				if (Dhsed>Dhsedmax) {
					Dhsed = Dhsedmax;
					n_limited_sed++;
				}
			}
			/*Erosion*/
			if (Dhsed<0) {
				float Dhsedmin=MIN_2(minsorr-topo[ind]+1e-2, 0);
				if (Dhsed<Dhsedmin) {
					Dhsed = Dhsedmin;
					n_limited_eros++;
				}
			}
			/*limit*/
			d_mass = -THICK2SEDMASS_1D(cfg, Dhsed);
			break;
		  case 'E':
			/*SADDLE RIVER-LIKE EROSION*/
			if (topo[ind]>ctx->sea_level || !AT_BORDER_1D(ind)) {
				Dhsed = -MASS2SEDTHICK_1D(cfg, d_mass);
				if (Dhsed>0) {
					PRINT_WARNING("Fluvial Transport: node <%d> (topo=%.1f) is a saddle and shouldn't deposit %.1e m.", ind, topo[ind], Dhsed); 
					Dhsed=0;
				}
				else
					Dhsed = MAX_2(Dhsed, MIN_2(minsorr-topo[ind]+1, 0));
				d_mass = -THICK2SEDMASS_1D(cfg, Dhsed);
			}
			/*Sea outlets (border nodes below sea_level) are already sedimented in Lake_Fill*/
			else 
				d_mass = 0;
			break;
		  default:
		  	PRINT_ERROR("Fluvial Transport: node [%d] has no defined type.", ind);
		}

		/*Adds results to the topo and the next load Dq and removes/adds material to the Blocks*/
		if (d_mass<0) {
			/*SEDIMENTATION, limit d_mass with the masstr in this node*/
			d_mass = MAX_2(d_mass, -drainage[ind].masstr*dt_fv);
			Sediment (cfg, ctx, -d_mass, ind);
		}
		if (d_mass>0) {
			/*EROSION.*/
			Erode    (cfg, ctx, d_mass, ind);
		}

		/*Adds the mass increment to the transferring mass contained in this cell: */
		drainage[ind].masstr += d_mass/dt_fv;
		
		/*Transfers suspended solid mass.*/
		if (drainage[ind].type!='L') {
		  if (IN_DOMAIN_1D(dind)) {
		    int ild = drainage[dind].lake;
		    float hl;
		    switch (drainage[dind].type) {
		    	case 'L':
			    /*Transfer within lake nodes is performed in Lake_Fill*/
			    drainage[dind].masstr += drainage[ind].masstr;
			    /*If draining to an OPEN lake:*/
			    if (Lake[ild].n_sd) {
				/*Check: this can happen when a node transferring to a lake is eroded below the lake level or when the lake node was already deposited and became higher than the tributary node.*/
			    	if (IN_DOMAIN_1D(drainage[dind].dr))  if (topo[drainage[dind].dr] > topo[ind])
			    		PRINT_WARNING("node %d transferring to lake in %d is < than the saddle %d.", ind, dind, drainage[dind].dr);
				hl = topo[ind]-1.;  /*MIN_2 (Lake[ild].alt+1., topo[ind]-1.);??*/
				Lake_Fill (Lake, cfg, ctx, dind, hl, dt_fv);
				if (IN_DOMAIN_1D(drainage[dind].dr)) 
					/*Next line NOT commented in TISC!!*/
					/*Keep commented, otherwise it doubles the presence of masstr, not preserving mass*/
					/*drainage[drainage[dind].dr].masstr += drainage[dind].masstr*/;
				else 
					total_lost_sed_mass += drainage[dind].masstr * dt_fv;
			    }
			    /*If draining to a CLOSED lake:*/
			    else {
			    	hl = topo[ind]-1.;
			    	Lake_Fill (Lake, cfg, ctx, dind, hl, dt_fv);
				/*Check: should be no sediment left.*/
				/*Does not work as in TISC, because dind is changed in Lake_Fill*/
			    }
			    /*Next line NOT commented in TISC!! See comment a few lines above*/
			    //drainage[dind].masstr = 0;
			    break;
			case 'R':
			case 'E':
			    drainage[dind].masstr += drainage[ind].masstr;
			    break;
		     }
		  }
		  else {
			/*Transfers out of model.*/
			total_lost_sed_mass += drainage[ind].masstr * dt_fv;
		  }
		}
	    }
	if (cfg->verbose_level>=3) fprintf(stdout, "\b\b\b"); fflush(stdout);

	if (cfg->verbose_level>=3) {
		if (n_limited_sed)  PRINT_WARNING("sedimentation was limited by discretisation %d times.",  n_limited_sed);
		if (n_limited_eros) PRINT_WARNING("erosion was limited by discretisation %d times.", n_limited_eros);
	}

	return (1);
}



int Lake_Fill (
	struct LAKE_INFO_1D *Lake, ModelConfig *cfg, ModelContext *ctx, 
	int ind, 		/*Lake node receiving the sediments*/
	float hl, 		/*Limit of sedimentation heigh*/
	float dt_fv)
{
	/*
	  THIS ROUTINE TREATS THE PROBLEM OF LAKE SEDIMENT INFILLING.
	  The strategy is to call this routine when a lake node is
	  to receive sediments. Then, the sediments and water of
	  this node is equally distributed to all the lake neighbours,
	  depositing a part of the sediments which is limited by
	  the lake height.
	  No sedimentation occurs in the saddles, which are eroded
	  apart as rivers.
	  ind is index number of the lake node receiving sediments.
	*/

	int	i, j, k, m, il, incr=0;
	float 	d_mass, Dhsedmax, Dhsed, l_fluv_sedim_aux=l_fluv_sedim; /* Global from tao.h */

	if (!drainage[ind].masstr) return(0);

	il = drainage[ind].lake;

	/*Sediments very efficiently if the lake is small and endorheic (!!)*/
	if (!Lake[il].n_sd && Lake[il].n<cfg->Nx/300) l_fluv_sedim_aux = .1;

	/*<Determine the incr direction of filling*/
	if (IN_DOMAIN_1D(ind-1)) if (drainage[ind-1].lake==drainage[ind].lake) incr=-1;
	if (IN_DOMAIN_1D(ind+1)) if (drainage[ind+1].lake==drainage[ind].lake) incr=+1;

	/*Sediment following the incr direction*/
	for (i=0;;) {
	    if (drainage[ind].masstr<.00001) break;
	    d_mass = cfg->dx / l_fluv_sedim_aux * drainage[ind].masstr * dt_fv;
	    d_mass = MIN_2(d_mass, drainage[ind].masstr*dt_fv);
	    /*Sediment in the next lake node*/
	    /*Limit d_mass with the lake height, with a little slope from the river mouth*/
	    Dhsedmax = MAX_2(hl+1.-i*cfg->dx*.0001-ctx->topo[ind], 0);
	    Dhsed = MASS2SEDTHICK_1D(cfg, d_mass);
	    Dhsed = MIN_2(Dhsed, Dhsedmax);
	    d_mass = THICK2SEDMASS_1D(cfg, Dhsed);
	    Sediment (cfg, ctx, d_mass, ind);
	    i++;
	    drainage[ind].masstr -= d_mass/dt_fv;
	    if (IN_DOMAIN_1D(ind+incr)) {
	    	if (drainage[ind+incr].lake==drainage[ind].lake) {
			drainage[ind+incr].masstr = drainage[ind].masstr;
			ind+=incr;
		}
		else break;
	    }
	    else break;
    	    if (!incr) break;
	}

	/*If the lake is endorheic and there is still some transported mass left, then deposit the rest evenly in the lake.*/
	if (!Lake[il].n_sd && drainage[ind].masstr > .01) {
		d_mass = drainage[ind].masstr*dt_fv;
		for (i=0; i<Lake[il].n; i++) {
			//factor = (fabs(((float) Lake[il].cell[i]-(float) ind)/((float) Lake[il].n)) - .5) * .01 + 1;
			//d_mass_node = d_mass/Lake[il].n * factor;
			Sediment (cfg, ctx, d_mass/Lake[il].n, Lake[il].cell[i]);
		}
		if (MASS2SEDTHICK_1D(cfg, d_mass)/Lake[il].n>10.) PRINT_WARNING("filling closed lake %d in [%d] in rough way: %.1f kg/s.", il, ind, drainage[ind].masstr);
		drainage[ind].masstr = 0;
	}

	return(1);
}




int Landslide_Transport (ModelConfig *cfg, ModelContext *ctx, float critical_slope)
{
	float	Dheros, dl;
	int 	i, j, k, n_iters;
	
	/*
	  COMPUTES THE LOAD, HEIGHT & Block-THICKNESS INCREMENTS DUE TO 
	  LANDSLIDING. This process approaches short scale 
	  transport processes. Assumes that slope cannot exceed a critical 
	  slope.
	*/

	if (!critical_slope) return (0);
	
	if (cfg->verbose_level>=3) fprintf(stdout, " ls");	fflush(stdout);

	calculate_topo(cfg, ctx, ctx->topo);

	n_iters = MAX_2(floor(ctx->dt/ctx->dt_eros+.5), 1); 
	PRINT_INFO("n_iters=%3d", n_iters);

	for (k=0; k<n_iters; k++) {
 	    /*Adds results to the height and the next load Dq and removes material from the Blocks*/
 	    for(j=0; j<cfg->Nx; j++)  {
	    	int imaxslope;
		float maxdiff;

		/*Calculate max slope in the 8 sorrounding directions*/
		imaxslope=-1; maxdiff=0;
		if (IN_DOMAIN_1D(j-1)) {
			if (maxdiff < ctx->topo[j]-ctx->topo[j-1]) {
				maxdiff = ctx->topo[j]-ctx->topo[j-1];
				imaxslope = j-1;
			}
		}
		if (IN_DOMAIN_1D(j+1)) {
			if (maxdiff < ctx->topo[j]-ctx->topo[j+1]) {
				maxdiff = ctx->topo[j]-ctx->topo[j+1];
				imaxslope = j+1;
			}
		}

		if (imaxslope>=0) {
    		    if (critical_slope < maxdiff/cfg->dx) {
			    /*it should preserve mass instead of using denscrust!!*/
    			    Dheros = maxdiff/4; /*reduces slope by 1/2*/
			    Erode    (cfg, ctx, THICK2SEDMASS_1D(cfg, Dheros), j);
			    Sediment (cfg, ctx, THICK2SEDMASS_1D(cfg, Dheros), imaxslope);
    		    }
		}
	    }
	}

	return (1);
}




int Sediment (ModelConfig *cfg, ModelContext *ctx, double d_mass, int ind) 
{
	/*
	  Adds dh_sed to the highest Block.
	  d_mass > 0   (mass of deposited seds., subtracted from the river load, does not include the water deposited with the seds)
	  dh_sed > 0   (thickness of sediments deposited)
	*/

	float dh_sed;

	dh_sed = MASS2SEDTHICK_1D(cfg, d_mass);
	if (dh_sed < -10) PRINT_WARNING("trying to sediment negative mass: %f m", dh_sed);
	/*Increment load, Blocks and topo*/
	Dq[ind] +=  dh_sed * g * (cfg->denssedim-cfg->densenv);
	Blocks[ctx->numBlocks-1].thick[ind] += dh_sed;
	ctx->topo[ind] += dh_sed;
	/*record of eros/sed is performed in kg*/
	eros_now[ind]      -= d_mass ;
	total_erosion[ind] -= d_mass ;

	ctx->total_sed_mass += d_mass;
	return (1);
}


int Erode (ModelConfig *cfg, ModelContext *ctx, double d_mass, int ind) 
{
	/*
	  Erode a certain amount of rock mass from the uppermost Blocks.
	  Increment load, Blocks and topo.
	  d_mass  > 0   (eroded mass)
	  dh_eros > 0   (thickness of material eroded)
	*/

	int k;
	float dh_eros=0, dh_eros_Block, mass_per_m2=d_mass/cfg->dx/cfg->riverbasinwidth;

	for (k=ctx->numBlocks-1; mass_per_m2>0 && k>=0; k--) {
		if (Blocks[k].density == cfg->denssedim) {
			dh_eros_Block = MIN_2(Blocks[k].thick[ind], MASS2SEDTHICK_1D(cfg, mass_per_m2)*cfg->dx*cfg->riverbasinwidth);
			mass_per_m2 -= THICK2SEDMASS_1D(cfg, dh_eros_Block) / cfg->dx/cfg->riverbasinwidth;
			ctx->total_sed_mass -= THICK2SEDMASS_1D(cfg, dh_eros_Block);
		}
		else {
			dh_eros_Block = MIN_2(Blocks[k].thick[ind], fabs(mass_per_m2/Blocks[k].density));
			mass_per_m2 -= dh_eros_Block * Blocks[k].density;
			ctx->total_bedrock_eros_mass += dh_eros_Block * Blocks[k].density * cfg->dx*cfg->riverbasinwidth;
		}
		Blocks[k].thick[ind] -= dh_eros_Block;
		dh_eros += dh_eros_Block;
		Dq[ind] -=  g * dh_eros_Block * (Blocks[k].density - cfg->densenv);
	}
	/*Erode basement*/
	if (mass_per_m2>0) {
		Blocks_base[ind] -= mass_per_m2/cfg->denscrust;
		dh_eros += mass_per_m2/cfg->denscrust;
		Dq[ind] -=  g * mass_per_m2*(cfg->denscrust-cfg->densenv)/cfg->denscrust;
		ctx->total_bedrock_eros_mass +=  mass_per_m2 * cfg->dx*cfg->riverbasinwidth;
	}
	ctx->topo[ind] -= dh_eros;
	/*record of eros/sed is performed in kg*/
	eros_now[ind]      += d_mass ;
	total_erosion[ind] += d_mass;

	return (1);
}




int Add_Node_To_Lake (int ncell, int i_lake)
{
	/*	
	  Allocates memory for a new lake node.
	  New node is added at the end of the lake's node list.
	*/
	int i, il;

	il = fabs((float) i_lake);

	/*fprintf (stderr, "\nADD node   %d to lake %d", ncell, il);*/
	for (i=0; i<Lake[il].n; i++) 
		if (ncell==Lake[il].cell[i]) {
		  	PRINT_WARNING("tryed to repeat lake node <%d>.", ncell);
			return(0);
		}

	Lake[il].n++;
	Lake[il].cell = realloc(Lake[il].cell, Lake[il].n*sizeof(int));
	Lake[il].cell[Lake[il].n-1] = ncell;

	/*
	  Add information to the drainage grid.
	  Use the same sign as the nodes in the lake already had.
	  Negative sign notes the undefined character.
	*/
	drainage[ncell].lake = -il;

	return (1);
}



int Add_Saddle_To_Lake (int sd, int tr, int i_lake)
{
	/*	
	  Allocates memory for a new lake node.
	  New saddle is added at the end of the lake's saddle list.
	*/
	int i, il;

	il = fabs((float) i_lake);

	/*fprintf (stderr, "\nADD saddle %d to lake %d", sd, il);*/
	for (i=0; i<Lake[il].n_sd; i++) 
		if (sd==Lake[il].sd[i]) {
		  	PRINT_WARNING("tryed to repeat lake saddle <%d>.", sd);
			return(0);
		}

	Lake[il].n_sd++;
	Lake[il].sd = realloc(Lake[il].sd, Lake[il].n_sd*sizeof(int));
	Lake[il].sd[Lake[il].n_sd-1] = sd;
	drainage[sd].dr = tr;

	return (1);
}




int Attempt_Delete_Node_From_Lake (ModelConfig *cfg, ModelContext *ctx, int ncell)
{
	/*
	  Checks if it is necessary to divide lake when deleting node by evaporation.
	  Returns the number of subdivided lakes (>=2, 0 means no subdivision).
	*/
	int 	i, j, k, il, splitting=0;

	il = drainage[ncell].lake;
	/*Locate neighbours of the same lake*/

	if (!il) {PRINT_ERROR("the node is not in a lake"); return (0);}
	if (Lake[il].n<=1) return (0);
	/*Don't delete the node if there are still higher saddles in the lake*/
	if (Lake[il].n_sd) if (topo[Lake[il].sd[0]]>topo[ncell]) return (0);
	/*Don't delete sea nodes*/
	if (Lake[il].n_sd) if (topo[Lake[il].sd[0]]<sea_level && AT_BORDER_1D(Lake[il].sd[0])) return (0);

	/*
	  If the lake of this node was subdivided then don't delete this node.
	  If necessary, it will be deleted in Divide_Lake.
	*/
	splitting = Divide_Lake(cfg, ctx, ncell);
	if (splitting) {
		if (drainage[ncell].lake<0) {
			Define_Lake(drainage[ncell].lake);
		}
		return (splitting);
	}

	/*The node is now a river. Remove it from lake.*/
	Delete_Node_From_Lake (cfg, ctx, ncell);

	return (splitting);
}


int Deallocate_Lake (ModelContext *ctx, int i_lake)
{
	/*	
	  Frees memory of a lake
	*/

	int i, j, il, kjhjhfk;

	il = fabs((float) i_lake);

	/*fprintf(stderr, "\nDeleting lake %d (of %d).", i_lake, ctx->nlakes);*/
	free(Lake[il].cell);
	free(Lake[il].sd);

	for (i=il+1; i<=ctx->nlakes; i++) {
		for (j=0; j<Lake[i].n; j++) 
			drainage[Lake[i].cell[j]].lake = 
				(i-1) * fabs((float) drainage[Lake[i].cell[j]].lake)/drainage[Lake[i].cell[j]].lake;
		Lake[i-1].n = Lake[i].n;
		Lake[i-1].n_sd = Lake[i].n_sd;
		Lake[i-1].cell = Lake[i].cell;
		Lake[i-1].sd = Lake[i].sd;
	}

	ctx->nlakes--;
	nlakes = ctx->nlakes;

	Lake = (struct LAKE_INFO_1D *) realloc(Lake, (ctx->nlakes+1)*sizeof(struct LAKE_INFO_1D));
	
	return(1);
}




int Define_Lake (int i_lake)
{
	/*	
	  Marks a lake as defined, i.e., all its nodes are already known.
	*/
	int i, il;

	il = fabs((float) i_lake);

	for (i=0; i<Lake[il].n; i++) {
		if ((int) fabs((float) drainage[Lake[il].cell[i]].lake) != il)
			PRINT_ERROR("lake number %d was expected instead of %d at node [%d]",
				il, drainage[Lake[il].cell[i]].lake, Lake[il].cell[i]);
		drainage[Lake[il].cell[i]].lake = il;
	}
	return(1);
}



int Delete_Node_From_Lake (ModelConfig *cfg, ModelContext *ctx, int ln)
{
	/*
	  ACTUALLY REMOVES THE NODE FROM THE LAKE STRUCTURE 
	  Deallocates memory for a lake node.
	  The rest of nodes are shifted in the lake's node list.
	  Finds new drainage for this node.
	*/

	int 	i, j, k, il, ild, i_node, i_saddle;
	int	imaxderneg=SIGNAL;
	int 	ls;
	float	maxderneg=1e9;

	il = drainage[ln].lake;
	i_node = Lake_Node_Number(ln);
	i_saddle = Lake_Saddle_Number(ln);

	if (il<=0) PRINT_ERROR("[%d] is not a lake (%d)", ln, il);

	/*Calculate derivates*/
	if (IN_DOMAIN_1D(ln+1)) { 
		float deriv = (ctx->topo[ln+1]-ctx->topo[ln])/cfg->dx;
		if ((deriv<0 && deriv<maxderneg) || (deriv==0 && deriv<=maxderneg && drainage[ln+1].lake)) 
			{imaxderneg=ln+1; maxderneg=deriv;}
	}
	if (IN_DOMAIN_1D(ln-1)) {
		float deriv = (ctx->topo[ln-1]-ctx->topo[ln])/cfg->dx;
		if ((deriv<0 && deriv<maxderneg) || (deriv==0 && deriv<=maxderneg && drainage[ln-1].lake)) 
			{imaxderneg=ln-1; maxderneg=deriv;}
	}

	/*Determine new drainage by draining through the maximum slope*/
	if (IN_DOMAIN_1D(drainage[ln].dr)) 
		ild = drainage[drainage[ln].dr].lake;
	else 	ild = 0;
	if (maxderneg==0 && (drainage[ln].type=='E' || ild)) {
		/*Leave it the same*/
	}
	else {
		if (imaxderneg != SIGNAL) {
			drainage[ln].dr = imaxderneg;
		}
		else {
			drainage[ln].dr = SIGNAL;
		}
	}
	drainage[ln].type = 'R';
	drainage[ln].lake = 0;
	/*Deallocates node memory*/
	Lake[il].n--;
	for (i=i_node; i<Lake[il].n; i++) {
		Lake[il].cell[i] = Lake[il].cell[i+1];
	}
	Lake[il].cell = realloc(Lake[il].cell, Lake[il].n*sizeof(int));
	/*Arrange everything if the node to be deleted is a saddle*/
	if (i_saddle>=0) {
		/*Remove and deallocate saddle*/
		Lake[il].n_sd--;
		for (i=i_saddle; i<Lake[il].n_sd; i++) {
		    Lake[il].sd[i] = Lake[il].sd[i+1];
		}
		Lake[il].sd = realloc(Lake[il].sd, Lake[il].n_sd*sizeof(int));
		/*Rearrange drainage in lake*/
		if (Lake[il].n_sd) {
		    /*Lake nodes draining to this saddle should drain to the another (nearest) saddle:*/
		    for (i=0; i<Lake[il].n; i++) {
		    	ls = Lake[il].cell[i]; 
		    	if (drainage[ls].dr == ln) {
		    	    int imindist=-1;
			    float mindist=1e24, dist;
		    	    for (k=0; k<Lake[il].n_sd; k++) {
		    		dist = fabs(Lake[il].sd[k]-ls) * cfg->dx;
		    		if (dist < mindist) {imindist=k; mindist=dist;}
		    	    }
		    	    drainage[ls].dr = Lake[il].sd[imindist];
		    	    drainage[ln].discharge -= drainage[ls].discharge;
		    	    drainage[drainage[ls].dr].discharge += drainage[ls].discharge;
		    	}
		    	if (drainage[ls].lake != il)
		    	    PRINT_ERROR("Delete_Node_From_Lake: 'Lake' %d (%dth of %d) and 'drainage' %d don't match in node [%d].", il, i, Lake[il].n, drainage[ls].lake, ls);
		    }
		}
 		/*If this was the last saddle then change drainage of all lake nodes*/
		else {
		    for (i=0; i<Lake[il].n; i++) {
		    	ls = Lake[il].cell[i]; 
		    	drainage[ls].dr = SIGNAL;
		    	drainage[ln].discharge -= drainage[ls].discharge;
		    }
		}
	}

	/*Delete lakes with only saddles*/
	if (Lake[il].n == Lake[il].n_sd && Lake[il].n>0) {
		for (j=0; j<Lake[il].n_sd; j++) {
			drainage[Lake[il].cell[j]].lake = 0;
			drainage[Lake[il].cell[j]].type = 'R';
		}
		Deallocate_Lake(ctx, il);
	}

	if (drainage[ln].lake<0) {
		PRINT_ERROR("lake node [%d]: lake %d", ln, drainage[ln].lake);
		Define_Lake(drainage[ln].lake);
	}

	return (1);
}



int Divide_Lake (ModelConfig *cfg, ModelContext *ctx, int ind) 
{
	/*
	  This routine divides a lake if the node that is to be 
	  removed from the lake (by evaporation) is 
	  the only connexion between two or more parts of the lake.
	*/
	int 	i, j, k, lake_delete, il, i_saddle, 
		indi, indj, 
		imaxderneg, 
		local_num_lakes=0;
	float 	lake_evaporation, maxderneg=0;
	bool 	became_overburden_1, became_overburden_2;

	il=drainage[ind].lake;
	if (Lake[il].n<3) return(0);

	/*Find if two local lakes result from division*/
	if (IN_DOMAIN_1D(ind-1)) {if (drainage[ind-1].lake == drainage[ind].lake) local_num_lakes++;}
	if (IN_DOMAIN_1D(ind+1)) {if (drainage[ind+1].lake == drainage[ind].lake) local_num_lakes++;}
	if (local_num_lakes<2) {return(0);}

	/*So, there IS splitting*/

	if (IN_DOMAIN_1D(i_saddle=Lake_Saddle_Number(ind))) {
		PRINT_ERROR("splitting node %d shouldn't be a saddle!", ind);
	}

	/*Distribute the nodes among the 2 new lakes and delete the original lake*/
	{
	    int new_lake1, new_lake2;
	    new_lake1 = New_Lake(cfg, ctx);  new_lake2 = New_Lake(cfg, ctx);
	    for (j=0; j<Lake[il].n; j++) {
		if (Lake[il].cell[j]<ind) { /*left side*/
		    Add_Node_To_Lake   (Lake[il].cell[j], new_lake1);
		    if (drainage[Lake[il].cell[j]].type=='E') {
			    Add_Saddle_To_Lake (Lake[il].cell[j], drainage[Lake[il].cell[j]].dr, new_lake1);
		    }
 		}
 		if (Lake[il].cell[j]>ind) { /*right side*/
		    Add_Node_To_Lake   (Lake[il].cell[j], new_lake2);
		    if (drainage[Lake[il].cell[j]].type=='E') {
			    Add_Saddle_To_Lake (Lake[il].cell[j], drainage[Lake[il].cell[j]].dr, new_lake2);
		    }
 		}
	    }
	    Define_Lake(new_lake1); Define_Lake(new_lake2);
	}
	Deallocate_Lake(ctx, il);

	/*
	  Check overburdening of the new lakes after separation 
	  and determine the role of the node to be deleted.
	*/
	became_overburden_1=became_overburden_2=false;
    	lake_evaporation=0;
	for (i=0; i<Lake[ctx->nlakes-1].n; i++) lake_evaporation += evaporation[Lake[ctx->nlakes-1].cell[i]] * cfg->riverbasinwidth * cfg->dx;
	if (lake_evaporation < Lake_Input_Discharge(cfg, ctx->nlakes-1)) became_overburden_1=true;
    	lake_evaporation=0;
	for (i=0; i<Lake[ctx->nlakes].n; i++) lake_evaporation += evaporation[Lake[ctx->nlakes].cell[i]] * cfg->riverbasinwidth * cfg->dx;
	if (lake_evaporation < Lake_Input_Discharge(cfg, ctx->nlakes))   became_overburden_2=true;
	if (became_overburden_1 && became_overburden_2) 
		PRINT_WARNING("no lake remained underburden among 2.");

	/*Find the maximum negative or null slope at ind*/
	maxderneg=0; imaxderneg = SIGNAL;
	/*Calculate derivates. I include deriv==0 because it can happen that the removed node is in a plane and needs to drain somewhere in that plane.*/
	if (IN_DOMAIN_1D(ind-1)) {
		float deriv = (ctx->topo[ind-1]-ctx->topo[ind])/cfg->dx;
		if (deriv<=0 && deriv<=maxderneg && !became_overburden_1) {
			imaxderneg=ind-1;
			maxderneg=deriv;
		}
	}
	if (IN_DOMAIN_1D(ind+1)) {
		float deriv = (ctx->topo[ind+1]-ctx->topo[ind])/cfg->dx;
		if (deriv<=0 && deriv<=maxderneg && !became_overburden_2) {
			imaxderneg=ind+1;
			maxderneg=deriv;
		}
	}
	/*no drainage may result only if ind is a node in the middle of the lake*/
	if (imaxderneg==SIGNAL) 
		PRINT_WARNING("Divide_Lake: no available drainage direction!");

	/*Change drainage in ind*/
	if (became_overburden_1 || became_overburden_2) {
	    /*Keep it as a lake node: add it to one of the new overburden lake.*/
	    /*Drain the new overburden lake towards the saddle.*/
	    if (became_overburden_1) {
		for (i=0; i<Lake[ctx->nlakes-1].n; i++) {
			drainage[Lake[ctx->nlakes-1].cell[i]].dr = ind;
			drainage[ind].discharge += drainage[Lake[ctx->nlakes-1].cell[i]].discharge;
		}
		Add_Node_To_Lake   (ind, ctx->nlakes-1);
		Add_Saddle_To_Lake (ind, imaxderneg, ctx->nlakes-1);
		drainage[ind].type = 'E';
		drainage[ind].dr = imaxderneg;
/*fprintf(stderr, "\nDIVIDING1 @ %d:  disch:%f, to  %d  nextlake: %d/%d  topos:%f,,%f", ind, drainage[ind].discharge, imaxderneg, drainage[ind+1].lake, nlakes, topo[ind], topo[ind+1]);
*/	    }
	    else {
		for (i=0; i<Lake[ctx->nlakes].n; i++) {
			drainage[Lake[ctx->nlakes].cell[i]].dr = ind;
			drainage[ind].discharge += drainage[Lake[ctx->nlakes].cell[i]].discharge;
		}
		Add_Node_To_Lake   (ind, ctx->nlakes);
		Add_Saddle_To_Lake (ind, imaxderneg, ctx->nlakes);
		drainage[ind].type = 'E';
		drainage[ind].dr = imaxderneg;
/*fprintf(stderr, "\nDIVIDING2 @ %d:  disch:%f, to  %d", ind, drainage[ind].discharge, imaxderneg);
*/	    }
	}
	else {
		/*No overburden sublakes => Not a lake node anymore: convert in river*/
		drainage[ind].lake = 0;
		drainage[ind].type = 'R';
		drainage[ind].dr = imaxderneg;
	}
	/*Define again both lakes and delete those having only saddles.*/
	for (i=ctx->nlakes-1; i<=ctx->nlakes; i++) {
		Define_Lake(i); 
		if (Lake[i].n == Lake[i].n_sd) {
			for (j=0; j<Lake[i].n; j++) {
				drainage[Lake[i].cell[j]].lake = 0;
				drainage[Lake[i].cell[j]].type = 'R';
			}
			Deallocate_Lake(ctx, i);
		}
	}

	return(2);
}



float Lake_Input_Discharge (ModelConfig *cfg, int ilake) 
{
	/*
	  The input discharge of a lake is the discharge of its nodes (except saddles).
	  This means that the inputs from tributaries in saddles do not contribute,
	  but the rain in the whole lake does.
	*/

	int 	i;
	float 	total_input_discharge=0;

	for (i=0; i<Lake[ilake].n; i++) {
		total_input_discharge += drainage[Lake[ilake].cell[i]].discharge;
		total_input_discharge += precipitation[Lake[ilake].cell[i]] * cfg->riverbasinwidth * cfg->dx;
	}
	for (i=0; i<Lake[ilake].n_sd; i++) 
		total_input_discharge -= drainage[Lake[ilake].sd[i]].discharge;

	return (total_input_discharge);
}



int Lake_Node_Number(int ind)
{
	int i, il, i_node;

	/*
	  Returns the order number of the node in Lake structure (>=0).
	  Returns -1 if not in a lake or not in the spected lake.
	*/

	il = drainage[ind].lake;
	if (il<=0) {PRINT_ERROR("node [%d] is not in a lake.", ind); return(-1);}
	for (i=0,i_node=-1; i<Lake[il].n; i++) 
		if (Lake[il].cell[i] == ind) 
			{i_node=i; break;}
	if (i_node==-1) {PRINT_ERROR("node is not in the expected lake."); return 0;}

	return(i_node);
}


int Lake_Saddle_Number (int ind)
{
	int i, il, i_saddle;

	/*
	  Returns the order number of the saddle in Lake structure (>=0).
	  Returns -1 if not in a lake or not a saddle of its lake.
	*/

	il = drainage[ind].lake;
	if (il<=0) {PRINT_ERROR("potential saddle [%d] is not in a lake.", ind); return(-1);}
	for (i=0,i_saddle=-1; i<Lake[il].n_sd; i++) 
		if (Lake[il].sd[i] == ind) 
			{i_saddle=i; break;}

	return(i_saddle);
}



int New_Lake (ModelConfig *cfg, ModelContext *ctx)
{
	/*	
	  Allocates memory for a new lake.
	  New lake is added at the end of the lake list.
	  Returns the number of the created lake.
	*/

	ctx->nlakes++;
	
	if (ctx->nlakes>cfg->Nx/3)
		PRINT_WARNING("Lots of lakes (%d)", ctx->nlakes);

	Lake = (struct LAKE_INFO_1D *) realloc(Lake, (ctx->nlakes+1)*sizeof(struct LAKE_INFO_1D));
	Lake[ctx->nlakes].n = 0;
	Lake[ctx->nlakes].n_sd = 0;
	Lake[ctx->nlakes].cell = NULL;
	Lake[ctx->nlakes].sd = NULL;

	nlakes = ctx->nlakes;
	return (ctx->nlakes);
}


int Calculate_Precipitation_Evaporation (ModelConfig *cfg, ModelContext *ctx)
{
	/*
	  Calculates the precipitation array in m/s
	*/
	int i;

	switch (cfg->hydro_model) {
	    case 1:
	    {
		float altitude;     int il;
		for (i=0; i<cfg->Nx; i++) {
		    altitude = ctx->topo[i];
		    if ((il = drainage[i].lake)) {
		    	/*!!*/
		    	altitude = ctx->topo[Lake[il].cell[Lake[il].n-1]];
		    }
		    precipitation[i] = MAX_2((rain+Krain*altitude), 0);
		    evaporation[i] = evaporation_ct;
		    if (CXrain) precipitation[i] *= MAX_2 (0, 1 + (cfg->x0+i*cfg->dx-(cfg->x0+cfg->xf)/2)/CXrain);
//precipitation[i] *= (1.5+cos(2*pi*(Time/Matosec/1.+.25)));
		}
		break;
	    }
	    case 2:
	    {
		float DX=CXrain, windvel=Krain /*[m/s]*/;
		/*smooth out precipitation*/
		for (i=0; i<cfg->Nx; i++) {
		     if (DX>cfg->dx) {
		     	 int j, jl, n_smooth;
		     	 float factor=2/DX/sqrt(3.1415927)*cfg->dx, weight, tweight=0, adimdist;
		     	 precipitation[i]=0;
		     	 n_smooth=ceil(2*DX/cfg->dx);
		     	 for (j=i-n_smooth; j<=i+n_smooth; j++) {
		     	     adimdist=(j-i)*cfg->dx/DX;
		     	     jl = LIMIT(j,0,cfg->Nx-1);
		     	     weight=exp(-pow(adimdist, 2));
		     	     tweight += weight;
		     	     precipitation[i] += Orographic_Precipitation(cfg, ctx, jl, windvel) * weight;
		     	 }
		     	 precipitation[i] /= tweight;
		     }
		     else {
		     	 precipitation[i] = Orographic_Precipitation(cfg, ctx, i, windvel);
		     }
		     evaporation[i] = evaporation_ct;
		}
		break;
	    }
    	    case 3:
    	    {
    		float DX=CXrain, windvel=Krain /*[m/s]*/;
		float *precip_aux;
		precip_aux = (float *) calloc(cfg->Nx, sizeof(float));
		Orographic_Precipitation_Evaporation_conservative(cfg, ctx, precip_aux, evaporation, windvel);
		/*smooth out precipitation*/
		for (i=0; i<cfg->Nx; i++) {
		     if (DX>cfg->dx) {
 		     	 int j, jl, n_smooth;
 		     	 float factor=2/DX/sqrt(3.1415927)*cfg->dx, weight, tweight=0, adimdist;
		     	 precipitation[i]=0;
 		     	 n_smooth=ceil(3*DX/cfg->dx);
 		     	 for (j=i-n_smooth; j<=i+n_smooth; j++) {
 		     	     adimdist=(j-i)*cfg->dx/DX;
 		     	     jl = LIMIT(j,0,cfg->Nx-1);
 		     	     weight=exp(-pow(adimdist, 2));
 		     	     tweight += weight;
		     	     precipitation[i] += precip_aux[jl] * weight;
		     	 }
 		     	 precipitation[i] /= tweight;
 		     }
 		     else {
 		     	 precipitation[i] = precip_aux[i];
 		     }
		}
    		free(precip_aux);
    		break;
    	    }
    	}
	return (1);
}


float Orographic_Precipitation (ModelConfig *cfg, ModelContext *ctx, int i, float windvel)
{
	/*
	  Returns precipitation in m/s.
	  Mostly from Roe et al., 2003, JGR. See also Masek et al., 1994, JGR.; Purves & Hulton, 2000, IJC.
	  It calculates the local vapor convergence as a function of temperature, slope, and wind velocity
	*/
	float 	es0=6.1078 /*mb*/, a=17.27, b=35.86 /*centigrades*/, 
		vapor_convergence, alpha0, alpha1, Ts, 
		slope=0, esat, topoC, topoR, topoL;
	int 	il;
	/*Surface temperature*/
	topoC = ctx->topo[i];   if ((il = drainage[i].lake)) topoC = ctx->topo[Lake[il].cell[Lake[il].n-1]];
	Ts = TEMPERATURE_GROUND(topoC); 
		/*Tetens formula for Clausius-Clapeyron, giving the saturation vapor pressure in the surface*/
	esat = es0 * exp(a*(Ts-TEMP_FREEZE_WATER)/(Ts-b));
	if (i>0 && i<cfg->Nx-1) {
		topoR=ctx->topo[i+1];   if (il=drainage[i+1].lake) topoR = ctx->topo[Lake[il].cell[Lake[il].n-1]];
		topoL=ctx->topo[i-1];   if (il=drainage[i-1].lake) topoL = ctx->topo[Lake[il].cell[Lake[il].n-1]];
		slope = (topoR - topoL) / cfg->dx / 2;
	}
	else {
		if (i==0) {
			topoR = ctx->topo[i+1];   if (il=drainage[i+1].lake) topoR = ctx->topo[Lake[il].cell[Lake[il].n-1]];
			slope = (topoR - topoC) / cfg->dx;
		}
		if (i==cfg->Nx-1) {
			topoL = ctx->topo[i-1];   if (il=drainage[i-1].lake) topoL = ctx->topo[Lake[il].cell[Lake[il].n-1]];
			slope = (topoC - topoL) / cfg->dx;
		}
	}
	if (drainage[i].lake) slope = 0;
	alpha0 = rain / es0; /*arbitrary value so that convergence=rain at plains at 0 Celcius*/
	alpha1 = 110 / es0;  /*m yr-1 / (m s-1)*/;
	/*vapor budget*/
	vapor_convergence = (alpha0 + alpha1/secsperyr * windvel * slope) * esat;

	return ((vapor_convergence>0)? vapor_convergence : 0);
}



int Orographic_Precipitation_Evaporation_conservative (ModelConfig *cfg, ModelContext *ctx, float *precip_aux, float *evaporation, float windvel)
{
	/*
	  Returns precipitation and evaporation in m/s.
	  Mostly from Roe et al., 2003, JGR, modified by Garcia-Castellanos, 2007, EPSL. See also Masek et al., 1994, JGR.
	  It calculates the local vapor convergence as a function of temperature, slope, and wind velocity
	*/
	float 	es0=610.78 /*Pa*/, 
		L, Rv=461.5, /*J/kg/K*/
		alpha0, beta=.3,  /*## !! better beta = 0.3, pero en sumision de altiplano esta con 1*/
		esat, topoC;
	int 	i, il, istart=0, incr=0;
	float 	Wcol=0, Wmax, 	/*Water content in a column, and maximum water content*/
		rh=.9, 		/*Relative humidity at the upwind boundary (incoming air rel. humidity)*/
		z, dz=50, dtwind;
	dtwind = cfg->dx/fabs(windvel); /*time spent by the air in a cell*/
	alpha0 = rain/rh;   /*arbitrary value so that convergence=rain at plains with a saturated column*/
	if (windvel>0)	{istart=0;    incr=+1;}
	if (windvel<0)	{istart=cfg->Nx-1; incr=-1;}

	/*Start upwind (windward) and go downwind (leeward) */
	if (windvel) for (i=istart; i>=0 && i<cfg->Nx; i+=incr) {
	    /*calculate maximum possible water content (saturation) in column i, in m*/
	    Wmax=0;
	    topoC = ctx->topo[i];   if ((il = drainage[i].lake)) topoC = ctx->topo[Lake[il].cell[Lake[il].n-1]];
	    for (z=0; z<10000; z+=dz) {
		float temp_air;
		temp_air = TEMPERATURE_AIR(topoC, z);
		L = 2.4995e6+(temp_air-TEMP_FREEZE_WATER)*2359;
		esat = es0*exp(L/Rv*(1/TEMP_FREEZE_WATER - 1/temp_air));
		Wmax += esat/temp_air/Rv/cfg->denswater*dz; /*Blocks: m of water*/
		//fprintf(stderr, "\nz: %f m  temp=%f K  Wmax=%f m", z, temp_air, Wmax);
	    }
	    //fprintf(stderr, "\nz: %f m  Wmax=%f m  alpha0=%.2e", z, Wmax, alpha0*secsperyr);
	    if ((i==0 && windvel>0) || (i==cfg->Nx-1 && windvel<0)) {
		    /*calculate water content in incomming air column*/
		    Wcol = rh * Wmax;
		    PRINT_INFO("Incoming air water content: %.3f m => %.2e m2/s", Wcol, Wcol*windvel);
	    }
 	    /*calculate precipitation*/
 	    precip_aux[i] = alpha0 * Wcol / Wmax;
	    /*limit precipitation to at least 0 or the excess water and at most the available water in column*/
	    precip_aux[i] = MAX_2(0, precip_aux[i]);
	    precip_aux[i] = MAX_2((Wcol-Wmax)/dtwind, precip_aux[i]);
	    precip_aux[i] = MIN_2(Wcol/dtwind, precip_aux[i]);
 	    /*calculate change in water content in column (in m of water)*/
	    /*fprintf(stderr, "\nprecip: %e mm/yr   %e    %e m/s  %e", precip_aux[i], Wcol/dtwind, windvel, slope);*/
	    Wcol -= precip_aux[i] * dtwind;
	    evaporation[i] = evaporation_ct * (1+beta*fabs(windvel)) * (Wmax-Wcol)/Wmax; 
	    evaporation[i] = MAX_2(0, evaporation[i]); 
	    if (drainage[i].lake) {
	    	Wcol += evaporation[i] * dtwind;
	    }
	    //fprintf(stderr, "\nx=%f \tP,E: %.2f , %.2f m/yr   Wcol=%.3e m  Wmax=%.3e m   w=%.2f m/s", i*dx+x0, precip_aux[i]*secsperyr, evaporation[i]*secsperyr, Wcol, Wmax, windvel);
 	}

	return (1);
}



int Unify_Lakes (ModelConfig *cfg, ModelContext *ctx, int i_lake, int i_lake_to_delete)
{
	int i, il, ild;

	il = fabs((float) i_lake);
	ild = fabs((float) i_lake_to_delete);

	/*fprintf(stderr, "\nUnifying lakes %d and %d (to be deleted) of %d.", i_lake, i_lake_to_delete, nlakes);*/

	Lake[il].n    += Lake[ild].n;
	Lake[il].n_sd += Lake[ild].n_sd;
	Lake[il].cell    = realloc(Lake[il].cell, Lake[il].n*sizeof(int));
	Lake[il].sd = realloc(Lake[il].sd, Lake[il].n_sd*sizeof(int));
	for (i=0; i<Lake[ild].n; i++) {
		Lake[il].cell[Lake[il].n-Lake[ild].n+i] = Lake[ild].cell[i];
	}
	for (i=0; i<Lake[ild].n_sd; i++) {
		Lake[il].sd[Lake[il].n_sd-Lake[ild].n_sd+i] = Lake[ild].sd[i];
	}

	/*Changes the associated drainage lake signal*/
	for (i=0; i<Lake[ild].n; i++) drainage[Lake[ild].cell[i]].lake = i_lake;

	Deallocate_Lake(ctx, ild);

	return (1);
}
