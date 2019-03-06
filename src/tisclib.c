/*
GENERAL  SUBS  LIBRARY  FOR  tisc.c
Daniel Garcia-Castellanos
*/



int Allocate_Memory_for_external_use()
{
	/* Allocates dynamic memory for the arrays called in call_tisc.c and initializes them to zero*/
	int i, j;
		
	PRINT_INFO("TISC memo initialisation.");

	topo		= alloc_matrix(Ny, Nx);
	Blocks_base	= alloc_matrix(Ny, Nx);
	q		= alloc_matrix(Ny, Nx);
	D		= alloc_matrix(Ny, Nx);
	Dq		= alloc_matrix(Ny, Nx);
	w		= alloc_matrix(Ny, Nx);

	Blocks =    	(struct BLOCK *) calloc(NmaxBlocks, sizeof(struct BLOCK));


	if (hydro_model) {
		sortcell =    	(struct GRIDNODE *) calloc(Nx*Ny, sizeof(struct GRIDNODE));
		for (j=0; j<Nx; j++)  for (i=0; i<Ny; i++) {sortcell[i*Nx+j].row=i; sortcell[i*Nx+j].col=j;}
		evaporation      	= alloc_matrix(Ny, Nx);
		precipitation      	= alloc_matrix(Ny, Nx);
		precipitation_snow 	= alloc_matrix(Ny, Nx);
		precipitation_file  	= alloc_matrix(Ny, Nx); 
		for (j=0; j<Nx; j++) for (i=0; i<Ny; i++) precipitation_file[i][j] = -9999;
		drainage =    	(struct DRAINAGE **) calloc(Ny, sizeof(struct DRAINAGE *));
		for (i=0; i<Ny; i++) drainage[i] = (struct DRAINAGE *) calloc(Nx, sizeof(struct DRAINAGE));
		lake_former_step = alloc_matrix_int(Ny, Nx);
	}
	if (erosed_model) {
		eros_now	= alloc_matrix(Ny, Nx);
		accumul_erosion	= alloc_matrix(Ny, Nx);
	}

	/*Allocation for Lake #0, which is not used.*/
	Lake = (struct LAKE_INFO *) calloc(1, sizeof(struct LAKE_INFO));
	fflush (stdout);

	return(1);
}


int Allocate_Memory()
{
	/* Allocates dynamic memory for the arrays and initializes them to zero*/
	int i, j;
		
	Allocate_Memory_for_external_use();

	Dw		= alloc_matrix(Ny, Nx);
	h_last_unit	= alloc_matrix(Ny, Nx);
	h_water		= alloc_matrix(Ny, Nx);
	EET		= alloc_matrix(Ny, Nx);

	{int i,j; for (i=0;i<Ny;i++) for (j=0; j<Nx; j++) EET[i][j]=Te_default;}
	if (hydro_model) {
		if (K_ice_eros) {
			ice_thickness = alloc_matrix(Ny, Nx);
			ice_sedm_load = alloc_matrix(Ny, Nx);
			ice_velx_sl = alloc_matrix(Ny, Nx);
			ice_vely_sl = alloc_matrix(Ny, Nx);
			ice_velx_df = alloc_matrix(Ny, Nx);
			ice_vely_df = alloc_matrix(Ny, Nx);			
		}
	}

#ifdef GENSPD_SOLVER
	if (solver_type == 'g') Alloc_genspd();
#endif
	return(1);
}



int calculate_topo(float **topo_new)
{
	float mean=0;
	/*Calculates current topography based on Blocks, Blocks_base and deflection*/
	PRINT_DEBUG("Entering");
	for (int i=0; i<Ny; i++) for (int j=0; j<Nx; j++) {
		float thickness_above=0;
		for (int i_Block=0; i_Block<numBlocks; i_Block++) 
			thickness_above += Blocks[i_Block].thick[i][j];
		topo_new[i][j] = Blocks_base[i][j]-w[i][j];
		for (int i_Block=0; i_Block<numBlocks; i_Block++) {
			thickness_above -= Blocks[i_Block].thick[i][j];
			topo_new[i][j] += Blocks[i_Block].thick[i][j];
			if (Blocks[i_Block].density==denssedim) topo_new[i][j] -= compaction(sed_porosity, compact_depth, thickness_above, thickness_above+Blocks[i_Block].thick[i][j]);
		}
		mean += topo_new[i][j];
	}
	mean /= Nx*Ny;
	return (mean);
}




int Delete_Block(int i_Block)
{
	int  	k;

	/*Deallocates one Block*/

	PRINT_DEBUG("Block %d out of %d", i_Block, numBlocks);

	free_matrix (Blocks[i_Block].thick, Ny);
	if (Blocks[i_Block].type == 'V') {
		free_matrix (Blocks[i_Block].vel_x,    Ny);
		free_matrix (Blocks[i_Block].vel_y,    Ny);
		free_matrix (Blocks[i_Block].visc,     Ny);
		free_matrix (Blocks[i_Block].viscTer,  Ny);
	}
	if (Blocks[i_Block].type == 'S') {
		free_matrix (Blocks[i_Block].detr_ratio,  Ny);
		free_matrix (Blocks[i_Block].detr_grsize, Ny);
	}
	for (k=i_Block; k<numBlocks-1; k++) Blocks[k] = Blocks[k+1];
	numBlocks--;
	if (i_Block<i_Block_insert) i_Block_insert--;
	return(1);
}




int insert_new_Block(int num_new_Block)
{
	struct BLOCK	Block_aux;

	/*Creates a new Block and increments numBlocks by 1.
		num_new_Block ranges from 0 to numBlocks.
		Blocks number num_new_Block and above are shifted upwards.
		If num_new_Block == numBlocks then a new block is created on top of all. 
	*/

	if (verbose_level>=2) fprintf(stdout, "  b"); fflush(stdout);
	PRINT_DEBUG("New Block being created: %d ; numBlocks= %d ; i_first_Block_load = %d ; i_Block_insert = %d", num_new_Block, numBlocks, i_first_Block_load, i_Block_insert);
	if (numBlocks>NmaxBlocks-5) PRINT_WARNING("Lots of Blocks! "); 

	Blocks[numBlocks].thick = 	alloc_matrix (Ny, Nx);
	Blocks[numBlocks].vel_x = 	alloc_matrix (1, 1);
	Blocks[numBlocks].vel_y = 	alloc_matrix (1, 1);
	Blocks[numBlocks].visc  = 	alloc_matrix (1, 1);
	Blocks[numBlocks].viscTer = 	alloc_matrix (1, 1);

	Block_aux = Blocks[numBlocks];
	for (int j_Block=numBlocks; j_Block>num_new_Block; j_Block--) {
		Blocks[j_Block] = Blocks[j_Block-1];
	}
	Blocks[num_new_Block] = Block_aux;

	/*Default properties*/
	Blocks[num_new_Block].type = '-';
	Blocks[num_new_Block].age = Time;
	Blocks[num_new_Block].density = 0;
	Blocks[num_new_Block].erodibility = erodibility;
	Blocks[num_new_Block].vel_x[0][0] = 0;
	Blocks[num_new_Block].vel_y[0][0] = 0;
	Blocks[num_new_Block].last_vel_time = Time;
	Blocks[num_new_Block].time_stop = 9999*Matosec;
	Blocks[num_new_Block].shift_x = 0;
	Blocks[num_new_Block].shift_y = 0;
	Blocks[num_new_Block].last_shift_x = 0;
	Blocks[num_new_Block].last_shift_y = 0;

	numBlocks++;

	return (1);
}



int gradual_Block()
{
	float 	Dhl;

	/*Non-instantaneous loading of a file load (distributed along time).*/

	/*interpolation can only last until next load because then the original load shape is lost*/
	if (!switch_gradual || Time>Blocks[i_Block_insert].time_stop-dt*.1) return(0);

	for (int i=0; i<Ny; i++) for (int j=0; j<Nx; j++) { 
		Dhl = h_last_unit[i][j]*dt/(Blocks[i_Block_insert].time_stop-Blocks[i_Block_insert].age);
		/*Increments the load for this time interval*/
		if (Blocks[i_Block_insert].type == 'H') 
			Blocks[i_Block_insert].thick[i][j] = 0;
		else {
			if (Dhl>=0) {
				/*fprintf(stderr, "\nnumBlocks=%d, i_Block_insert=%d, i=%d, j=%d, Dhl=%f, n_image=%d, time_stop=%f, age=%f", numBlocks, i_Block_insert, i, j, Dhl, n_image, Blocks[i_Block_insert].time_stop/secsperyr, Blocks[i_Block_insert].age/secsperyr);*/
				Blocks[i_Block_insert].thick[i][j] += Dhl;
			}
			else {
				float 	h_load_aux, h_load_aux2;
				int	k;
				h_load_aux = fabs((double) Dhl);
				for (k=i_Block_insert-1; h_load_aux>0 && k>=0; k--) {
					h_load_aux2 = MIN_2(Blocks[k].thick[i][j], h_load_aux);
					h_load_aux -= h_load_aux2;
					Blocks[k].thick[i][j] -= h_load_aux2;
				}
				/*k is the deepest eroded Block*/
				if (k==-1) {
					Blocks_base[i][j] -= h_load_aux;
				}
			}
		}
		Dq[i][j] += (Blocks[i_Block_insert].density-densenv)*g*Dhl;
	}

	return(1);
}




int defineLESalmostdiagonalmatrix (double **A, double *b, float **q, float **Dq, float **w, BOOL doing_visco)
{
	/*
		DEFINE LA MATRIZ EN BANDA 'A' Y EL TERMINO INDEPENDIENTE 'b' DE
		COEFICIENTES DEL SISTEMA DE ECUACIONES  
			A暖 = b 
		A[][Ny] es un elemento de la diagonal en la matriz en banda.
	*/

	register int 	i, j, 
			nodo, 			/* nodo = columna*Ny + fila */
			NDi=2*Ny, NDs=2*Ny ;
	double		dx4=dx*dx*dx*dx,   dy4=dy*dy*dy*dy,   
			dx2dy2=dx*dx*dy*dy,  dx2=dx*dx,  
			dy2=dy*dy,  dxdy=dx*dy ,
			D0, Dx, Dy, Dx2, Dy2, Dxy, Krest=0, bound_deflect;
	char		filename[MAXLENFILE];

	for (i=0; i<Ny*Nx; i++) {
		for (j=0; j<NDi+1+NDs; j++) A[i][j] = 0.;
		b[i]= 0.;
	}
	
	/* non B.C. nodes */
	for (i=2; i<Ny-2; i++)  for (j=2; j<Nx-2; j++) {
				 /*'nodo' es el numero de orden asociado al nodo de malla (i, j)  	*/
		nodo = j*Ny + i; /*Starting with node (0,0) y siguiendo con (1,0). 			*/
				 /*Corresponde tambien a la fila de A[][] en curso			*/
		/*Define the restoring force value.*/
		GET_KREST(Krest, Dq, i,j)

		D0 = 	D[i][j];						/*			*/
		Dx = 	D[i][j+1] - D[i][j-1];					/*			*/
		Dy = 	D[i-1][j] - D[i+1][j];					/*   Derivates de D	*/
		Dx2 = 	D[i][j+1] - 2*D[i][j] + D[i][j-1];			/*			*/
		Dy2 = 	D[i-1][j] - 2*D[i][j] + D[i+1][j];			/*			*/
		Dxy = 	D[i-1][j+1] - D[i+1][j+1] - D[i-1][j-1] + D[i+1][j-1];	/*			*/

/*Dx=Dy=Dx2=Dy2=Dxy=0;*/
/*Dx=Dx2=Dxy=0;*/
/*Dy=Dy2=Dxy=0;*/

/*Nodo asociado    Elemento de la matriz    Valor*/

/* i  , j-2 */	A[nodo][0] 	= + 1*D0 / dx4     					- 2*Dx / 4/dx4 ;

/* i-1, j-1 */	A[nodo][Ny-1] 	= + 2*D0 / dx2dy2  		- 2*Dx / 4/dx2dy2 + 2*Dy / 4/dx2dy2 						- 2*(1-nu)*Dxy / 16/dx2dy2 				- 2*Pxy / 4/dxdy ;
/* i  , j-1 */	A[nodo][Ny] 	= - 4*D0 / dx4     		- 4*D0 / 1/dx2dy2  	+ 4*Dx / 4/dx4  	+ 4*Dx2 / 4/dx4 		+ nu*Dy2 / dx2dy2 			+ Dx / dx2dy2	+ Px / dx2 ;
/* i+1, j-1 */	A[nodo][Ny+1] 	= + 2*D0 / dx2dy2  		- 2*Dx / 4/dx2dy2 - 2*Dy / 4/dx2dy2						+ 2*(1-nu)*Dxy / 16/dx2dy2 				+ 2*Pxy / 4/dxdy ;

/* i-2, j   */	A[nodo][2*Ny-2]	= + 1*D0 / dy4 	   					+ 2*Dy / 4/dy4 ;
/* i-1, j   */	A[nodo][2*Ny-1]	= - 4*D0 / dy4 	   		- 4*D0 / 1/dx2dy2  	- 4*Dy / 4/dy4  	+ 4*Dy2 / 4/dy4 		+ nu*Dx2 / dx2dy2  			- Dy / dx2dy2 	+ Py / dy2 ;
/* i  , j   */	A[nodo][2*Ny] 	= + 6*D0 / dx4  + 6*D0 / dy4 	+ 8*D0 / 1/dx2dy2				- 2*Dx2 / 1/dx4 - 2*Dy2 / 1/dy4	- 2*nu*Dy2 / dx2dy2  - 2*nu*Dx2 / dx2dy2 		- 2*Px / dx2  - 2*Py / dy2   +  Krest ;
/* i+1, j   */	A[nodo][2*Ny+1]	= - 4*D0 / dy4 	   		- 4*D0 / 1/dx2dy2  	+ 4*Dy / 4/dy4  	+ 4*Dy2 / 4/dy4 		+ nu*Dx2 / dx2dy2  			+ Dy / dx2dy2 	+ Py / dy2 ;
/* i+2, j   */	A[nodo][2*Ny+2]	= + 1*D0 / dy4 	   					- 2*Dy / 4/dy4 ;

/* i-1, j+1 */	A[nodo][3*Ny-1]	= + 2*D0 / dx2dy2  		+ 2*Dx / 4/dx2dy2 + 2*Dy / 4/dx2dy2						+ 2*(1-nu)*Dxy / 16/dx2dy2 				+ 2*Pxy / 4/dxdy ;
/* i  , j+1 */	A[nodo][3*Ny] 	= - 4*D0 / dx4     		- 4*D0 / 1/dx2dy2 	- 4*Dx / 4/dx4  	+ 4*Dx2 / 4/dx4 		+ nu*Dy2 / dx2dy2 			- Dx / dx2dy2	+ Px / dx2 ;
/* i+1, j+1 */	A[nodo][3*Ny+1]	= + 2*D0 / dx2dy2  		+ 2*Dx / 4/dx2dy2 - 2*Dy / 4/dx2dy2 						- 2*(1-nu)*Dxy / 16/dx2dy2 				- 2*Pxy / 4/dxdy ;

/* i  , j+2 */	A[nodo][4*Ny] 	= + 1*D0 / dx4     					+ 2*Dx / 4/dx4 ;


		/*Independent term*/
		switch (doing_visco) {
		    /*Elastic.*/
		    case 0:
			b[nodo] =
			    Dq[i][j] - Px*(w[i][j+1]-2*w[i][j]+w[i][j-1])/dx2 - Py*(w[i+1][j]-2*w[i][j]+w[i-1][j])/dy2;
			break;
		    /*Viscoelastic.*/
		    case 1:
			b[nodo] =
			    (q[i][j] - w[i][j]*Krest - Px*(w[i][j+1]-2*w[i][j]+w[i][j-1])/dx2 - Py*(w[i+1][j]-2*w[i][j]+w[i-1][j])/dy2) / tau /*Termino que deberia estar:+ Dq[i][j]/dt*/;
			break;
		}
	}


	/*BOUNDARY CONDITIONS:*/
	/*North*/
	switch (boundary_conds[0]) {
	    case '6':
	    case '0':
		for (i=2;  i <= Nx-3  ; i++) {
			bound_deflect=0; nodo=i*Ny;
			if (boundary_conds[0]=='6')  bound_deflect=Dq[0][i]/Krest;
			A[nodo][2*Ny] = 1;    b[nodo]   = bound_deflect ;
			if (boundary_conds[0]=='6')  bound_deflect=Dq[1][i]/Krest;
			A[nodo+1][2*Ny] = 1;  b[nodo+1] = bound_deflect ;
		}
		if (boundary_conds[3]=='5') {
			bound_deflect=0;
			if (boundary_conds[0]=='6')  bound_deflect=Dq[0][1]/Krest;
			A[Ny][2*Ny] = 1;  b[Ny] = bound_deflect;
		}
		if (boundary_conds[2]=='5') {
			bound_deflect=0;
			if (boundary_conds[0]=='6')  bound_deflect=Dq[0][Nx-2]/Krest;
			A[Ny*(Nx-2)][2*Ny] = 1;  b[Ny*(Nx-2)] = bound_deflect;
		}
		break;
	    default:
		for (i=1;  i <= Nx-2  ; i++) {
			nodo = i*Ny; 		/* Derivada en y nula */
			A[nodo][2*Ny] = 1;
			A[nodo][2*Ny+1] = -1;
			b[nodo] = 0;
		}
		for (i=2;  i <= Nx-3  ; i++) {
			nodo = i*Ny+1; 		/* Momento en y nulo */
			A[nodo][2*Ny] = -2;
			A[nodo][2*Ny+1] = +1;
			A[nodo][2*Ny-1] = +1;
			b[nodo] = 0 ;
		}
		break;
	}
	
	/*South*/
	switch (boundary_conds[1]) {
	    case '6':
	    case '0':
		for (i=3;  i <= Nx-2  ; i++) {
			bound_deflect=0; nodo=i*Ny-2;
			if (boundary_conds[1]=='6')  bound_deflect=Dq[Ny-2][i-1]/Krest;
			A[nodo][2*Ny] = 1;    b[nodo]   = bound_deflect ;
			if (boundary_conds[1]=='6')  bound_deflect=Dq[Ny-1][i-1]/Krest;
			A[nodo+1][2*Ny] = 1;  b[nodo+1] = bound_deflect ;
		}
		if (boundary_conds[3]=='5') {
			bound_deflect=0;
			if (boundary_conds[1]=='6')  bound_deflect=Dq[Ny-1][1]/Krest;
			A[2*Ny-1][2*Ny] = 1;  b[2*Ny-1] = bound_deflect;
		}
		if (boundary_conds[2]=='5') {
			bound_deflect=0;
			if (boundary_conds[1]=='6')  bound_deflect=Dq[Ny-1][Nx-2]/Krest;
			A[Ny*(Nx-1)-1][2*Ny] = 1;  b[Ny*(Nx-1)-1] = bound_deflect;
		}
		break;
	    default:
		for (i=1;  i <= Nx-2  ; i++) {
			nodo = (i+1)*Ny-1; 	/* Derivada en y nula */
			A[nodo][2*Ny] = 1;
			A[nodo][2*Ny-1] = -1;
			b[nodo] = 0 ;
		}
		for (i=2;  i <= Nx-3  ; i++) {
			nodo = (i+1)*Ny-2; 	/* Momento en y nulo */
			A[nodo][2*Ny] = -2;
			A[nodo][2*Ny+1] = +1;
			A[nodo][2*Ny-1] = +1;
			b[nodo] = 0 ;
		}
		break;
	}

	/*East*/
	switch (boundary_conds[2]) {
	    case '6':
	    case '0':
		for (i=0;  i <= Ny-1 ; i++) { 
			bound_deflect=0; 
			if (boundary_conds[2]=='6')  bound_deflect=Dq[i][Nx-2]/Krest;
			nodo = i+(Nx-2)*Ny;
			A[nodo][2*Ny] = 1;   b[nodo] = bound_deflect ;
			if (boundary_conds[2]=='6')  bound_deflect=Dq[i][Nx-1]/Krest;
			nodo = i+(Nx-1)*Ny;
			A[nodo][2*Ny] = 1;   b[nodo] = bound_deflect ;
		}
		break;
	    default:
		for (nodo=(Nx-1)*Ny;  nodo <= Nx*Ny-1 ; nodo++) { 
			A[nodo][2*Ny] = 1;
			A[nodo][1*Ny] = -1;	/*Last column: derivada en x nula*/
			b[nodo] = 0 ;
		}
		for (nodo=(Nx-2)*Ny+1;  nodo <= (Nx-1)*Ny-2 ; nodo++) { 
			A[nodo][2*Ny] = -2;
			A[nodo][1*Ny] = +1;	/*2nd last column: momento en x nulo*/
			A[nodo][3*Ny] = +1;
			b[nodo] = 0 ;
		}
		break;
	}

	/*West*/
	switch (boundary_conds[3]) {
	    case '6':
	    case '0':
		/*B.C. fixed: null deflection at boundary nodes and their neighbors*/
		for (i=0;  i <= Ny-1 ; i++) { 
			bound_deflect=0;
			if (boundary_conds[3]=='6')  bound_deflect=Dq[i][0]/Krest; 
			A[i][2*Ny] =    1;   b[i]    = bound_deflect ;
			if (boundary_conds[3]=='6')  bound_deflect=Dq[i][1]/Krest; 
			A[i+Ny][2*Ny] = 1;   b[i+Ny] = bound_deflect ;
		}
		break;
	    default:
		/* Free B.C.: zero moment and derivate at boundaries*/
		for (nodo=0;  nodo <= Ny-1 ; nodo++) {
			A[nodo][2*Ny] = 1; 
			A[nodo][3*Ny] = -1;   	/* First column:  derivada en x nula*/
			b[nodo] = 0 ;
		}
		for (nodo=Ny+1;  nodo <= 2*Ny-2 ; nodo++) { 
			A[nodo][2*Ny] = -2;
			A[nodo][1*Ny] = +1;   	/* 2nd column: momento en x nulo*/
			A[nodo][3*Ny] = +1;
			b[nodo] = 0 ;
		}
		break;
	}


	if (verbose_level>=3 && switch_write_file && Nx*Ny<200) {
		sprintf(filename, "%s.mtrz", projectname); 
		WriteAlmostDiagonalMatrix(A, b, Nx*Ny, filename, NDs, NDi);
	}

	return(1);

/*index of grid nodes:

				NORTH


	0     Ny    2Ny    3Ny    .      .      .    (Nx-1)意y


	1    Ny+1  2Ny+1  3Ny+1   .      .      .    (Nx-1)意y+1


	2    Ny+2  2Ny+2  3Ny+2   .      .      .    (Nx-1)意y+2


WEST	3    Ny+3  2Ny+3  3Ny+3   .      .      .    (Nx-1)意y+3	EAST


	.      .      .      .    .      .      .      .


	.      .      .      .    .      .      .      .


	.      .      .      .    .      .      .      .


	Ny-1  2Ny-1  3Ny-1  4Ny-1 .      .      .      Nx意y-1


				SOUTH
*/
}




int defineLESmatrix_for_mathlib (float *A, int *IA, int *JA, float *b, float **q, float **w, int *nonzeroes, BOOL doing_visco)
{
	/* 
		DEFINE LA MATRIZ 'A' Y EL TERMINO INDEPENDIENTE 'b' DE 
		COEFICIENTES DEL SISTEMA DE ECUACIONES  
		A暖 = b 

		Uses mathlib library.
		Only tied BC (#0) available with this solver.
	*/

	register int 	i, j, nodo, NDi=2*Ny, NDs=2*Ny, nz ;
	float		dx4=dx*dx*dx*dx,   dy4=dy*dy*dy*dy,   dx2dy2=dx*dx*dy*dy,  dx2=dx*dx,  dy2=dy*dy,  dxdy=dx*dy ,
			D0, Dx, Dy, Dx2, Dy2, Dxy, Krest;

	nz=0;
	IA[0]=1;
	for (j=0; j<Nx; j++)  for (i=0; i<Ny; i++)  {
					/*'nodo' es el numero de orden asociado al nodo de malla (i, j)  	*/
		nodo = j*Ny + i + 1;	/*comenzando por (0,0) y siguiendo con (1,0). 				*/
					/*Corresponde tambien a la fila de A[][] en curso, osea, el no. de eq.	*/

		/*Interior: donde gobierna la equacion*/
		if (i>=2 && i<Ny-2 && j>=2 && j<Nx-2) {
			/*Decide the restoring force value.*/
			if (switch_topoest) {
				/*If the current i,j knot is below the load 
				then the compensation density is densinfill.
				*/
				if (q[i][j]) 	Krest = (densasthen-densinfill)*g ;
				/*Otherwise use the sediment density.*/
				else		Krest = (densasthen-denssedim)*g ;
			}
			/*When the load subsides with the basament then the 
			displaced mantle is subtituted with densenv density 
			(given at parameters file).
			*/
			else	Krest = (densasthen-densenv)*g ;

			D0 = 	D[i][j];						/*			*/
			Dx = 	D[i][j+1] - D[i][j-1];					/*			*/
			Dy = 	D[i-1][j] - D[i+1][j];					/*   Derivadas de D	*/
			Dx2 = 	D[i][j+1] - 2*D[i][j] + D[i][j-1];			/*			*/
			Dy2 = 	D[i-1][j] - 2*D[i][j] + D[i+1][j];			/*			*/
			Dxy = 	D[i-1][j+1] - D[i+1][j+1] - D[i-1][j-1] + D[i+1][j-1];	/*			*/


			/*Elemento de la matriz    Valor*/
	
			A[nz]	= + D0 / dx4  	- 2*Dx / 4 / dx4 ;
			JA[nz] = nodo-2*Ny;	nz ++;

			A[nz]	= + 2 * D0 / dx2dy2 	- 2*Dx / 4 / dx2dy2 	+ 2*Dy / 4 / dx2dy2 	- 2*(1-nu)*Dxy / 16 / dx2dy2 	- Pxy / 4 / dxdy ;
			JA[nz] = nodo-Ny-1;	nz ++;
			A[nz] 	= - 4 * D0 / dx4 	- 4 * D0 / dx2dy2  	+ Dx / dx4  	+ Dx2 / dx4  	+ nu*Dy2 / dx2dy2 	+ Dx / dx2dy2	+ Px / dx2 ;
			JA[nz] = nodo-Ny; 	nz ++;
			A[nz]	= + 2 * D0 / dx2dy2 	- 2*Dx / 4 / dx2dy2 	- 2*Dy / 4 / dx2dy2	+ 2*(1-nu)*Dxy / 16 / dx2dy2 	+ Pxy / 4 / dxdy ;
			JA[nz] = nodo-Ny+1;	nz ++;

			A[nz]	= + D0 / dy4 	+ 2 * Dy / 4 / dy4 ;
			JA[nz] = nodo-2;	nz ++;
			A[nz]	= - 4 * D0 / dy4 	- 4 * D0 / dx2dy2  	- Dy / dy4  	+ Dy2 / dy4  	+ nu*Dx2 / dx2dy2  	- Dy  / dx2dy2 	+ Py / dy2 ;
			JA[nz] = nodo-1;	nz ++;
			A[nz]	= + 6 * D0 / dx4	+ 6 * D0 / dy4	+ 8 * D0 / dx2dy2	- 2 * Dx2 / dx4	- 2 * Dy2 / dy4	- 2 * nu * Dy2 / dx2dy2	- 2 * nu * Dx2 / dx2dy2    -2*Px / dx2    - 2*Py / dy2	       + Krest ;
			JA[nz] = nodo;	nz ++;
			A[nz]	= - 4 * D0 / dy4 	- 4 * D0 / dx2dy2  	+ Dy / dy4  	+ Dy2 / dy4  	+ nu*Dx2 / dx2dy2  	+ Dy  / dx2dy2 	+ Py / dy2 ;
			JA[nz] = nodo+1;	nz ++;
			A[nz]	= + D0 / dy4 	- 2 * Dy / 4 / dy4 ;
			JA[nz] = nodo+2;	nz ++;

			A[nz]	= + 2 * D0 / dx2dy2 	+ 2*Dx / 4 / dx2dy2 	+ 2*Dy / 4 / dx2dy2	+ 2*(1-nu)*Dxy / 16 / dx2dy2 	+ Pxy / 4 / dxdy ;
			JA[nz] = nodo+Ny-1;	nz ++;
			A[nz]	= - 4 * D0 / dx4 	- 4 * D0 / dx2dy2 	- Dx / dx4  	+ Dx2 / dx4  	+ nu*Dy2 / dx2dy2 	- Dx / dx2dy2	+ Px / dx2 ;
			JA[nz] = nodo+Ny;	nz ++;
			A[nz]	= + 2 * D0 / dx2dy2 	+ 2*Dx / 4 / dx2dy2 	- 2*Dy / 4 / dx2dy2 	- 2*(1-nu)*Dxy / 16 / dx2dy2 	- Pxy / 4 / dxdy ;
			JA[nz] = nodo+Ny+1;	nz ++;

			A[nz] 	= + D0 / dx4  	+ 2*Dx / 4 / dx4 ;
			JA[nz] = nodo+2*Ny;	nz ++;
			
			/*Termino independiente*/
			switch (doing_visco) {
				case 0:		b[nodo-1] = q[i][j] ;				/* Elastico.	*/
								break;
				case 1:		b[nodo-1] = ( q[i][j] - w[i][j] * Krest ) / tau ;	/* Viscoelastico.*/
								break;
			}
		}
		else {
			/* C.C. fija: deflexion cero en el borde (puntos exteriores y sus contiguos).	*/
			if (boundary_conds[0] != '5') {
				A[nz] = 1; 	b[nodo-1] = 0;
				JA[nz] = nodo;	nz ++;
			}
			/* C.C. libre.	*/
			else {
				if (j==0) {
					A[nz] = 1;	JA[nz] = nodo;		nz ++;
					A[nz] = -1;	JA[nz] = nodo+Ny;	nz ++;
					b[nodo-1] = 0 ;
				}
				if (j==Nx-1) {
					A[nz] = 1;	JA[nz] = nodo-Ny;	nz ++;
					A[nz] = -1;	JA[nz] = nodo;  	nz ++;
					b[nodo-1] = 0 ;
				}
				if (i==0 && j!=0 && j!=Nx-1) {
					A[nz] = 1;	JA[nz] = nodo;		nz ++;
					A[nz] = -1;	JA[nz] = nodo+1;	nz ++;
					b[nodo-1] = 0 ;
				}
				if (i==Ny-1 && j!=0 && j!=Nx-1) {
					A[nz] = 1;	JA[nz] = nodo-1;	nz ++;
					A[nz] = -1;	JA[nz] = nodo;  	nz ++;
					b[nodo-1] = 0 ;
				}
				
				if (j==1 && i!=0 && i!=Ny-1) {
					A[nz] = 1;	JA[nz] = nodo-Ny;	nz ++;
					A[nz] = -2;	JA[nz] = nodo;  	nz ++;
					A[nz] = 1;	JA[nz] = nodo+Ny;	nz ++;
					b[nodo-1] = 0 ;
				}
				if (j==Nx-2 && i!=0 && i!=Ny-1) {
					A[nz] = 1;	JA[nz] = nodo-Ny;	nz ++;
					A[nz] = -2;	JA[nz] = nodo;  	nz ++;
					A[nz] = 1;	JA[nz] = nodo+Ny; 	nz ++;
					b[nodo-1] = 0 ;
				}
				if (i==1 && j!=0 && j!=1 && j!=Nx-2 && j!=Nx-1) {
					A[nz] = 1;	JA[nz] = nodo-1;	nz ++;
					A[nz] = -2;	JA[nz] = nodo;  	nz ++;
					A[nz] = 1;	JA[nz] = nodo+1;	nz ++;
					b[nodo-1] = 0 ;
				}
				if (i==Ny-2 && j!=0 && j!=1 && j!=Nx-2 && j!=Nx-1) {
					A[nz] = 1;	JA[nz] = nodo-1;	nz ++;
					A[nz] = -2;	JA[nz] = nodo;  	nz ++;
					A[nz] = 1;	JA[nz] = nodo+1;  	nz ++;
					b[nodo-1] = 0 ;
				}
			}
		}
		IA[nodo] = nz + 1;
	}

	*nonzeroes = nz;
	return (1);
}



int Matrix_Info(float **xarxa, 						/*Matrix.*/
		int nx, int ny, 					/*Numbers of x & y knots.*/
		float dx, float dy, 					/*x & y interval.*/
		float *max, float *min, 				/*Pointers to the returned maximum & minimum values.*/
		float *maxcurvx, float *maxcurvy, float *maxcurvxy)	/*Curvatures.*/
{
	/*RETURNS CARACTERISTIC VALUES OF A GRID OF FLOATS:  EXTREMOS Y CURVATURAS	*/

	int 	i, j;
	float 	curvx, curvy, curvxy;

	*max=-1e10; *min=1e10;
	if (maxcurvx) *maxcurvx = 0;
	if (maxcurvy) *maxcurvy = 0;
	if (maxcurvxy) *maxcurvxy=0;
	for (i=0; i<ny; i++)  for (j=0; j<nx; j++)  {  
		*max = MAX_2(*max, xarxa[i][j]);
		*min = MIN_2(*min, xarxa[i][j]);
		if (j>0 && j<nx-1 && maxcurvx)   { 
			curvx = (xarxa[i][j+1] -2*xarxa[i][j] +xarxa[i][j-1]) /dx/dx;
			*maxcurvx = MAX_2(*maxcurvx, fabs(curvx));
		}
		if (i>0 && i<ny-1 && maxcurvy)   { 
			curvy = (xarxa[i-1][j] -2*xarxa[i][j] +xarxa[i+1][j]) /dy/dy;
			*maxcurvy = MAX_2(*maxcurvy, fabs(curvy));
		}
		if (i>0 && i<ny-1 && maxcurvxy) { 
			curvxy= (xarxa[i-1][j+1] +xarxa[i+1][j-1] -xarxa[i-1][j-1] + xarxa[i+1][j+1]) /4 /dx/dy;
			*maxcurvxy= MAX_2(*maxcurvxy, fabs(curvxy));
		}	
	}
	return(1);
}




int Perfil_info(float *perfil, int n, float *max, float *min)
{
	int 	i;

	*max=-1e10; *min=1e10; 
	for (i=0; i<n; i++) { 
		*max = MAX_2(*max, perfil[i]);
		*min = MIN_2(*min, perfil[i]);
	}
	return(1);
}




int Repare_Blocks()
{
	int  	i_Block_max_arrange;

	/*Avoids Blocks to have negative thickness or 0 volume*/

	PRINT_DEBUG("");

	i_Block_max_arrange = (switch_topoest)? i_first_Block_load : numBlocks;
	for (int i_Block=1; i_Block < i_Block_max_arrange; i_Block++) {
		for (int i=0; i<Ny; i++) for (int j=0; j<Nx; j++) {
			Blocks[i_Block-1].thick[i][j] = MAX_2(Blocks[i_Block-1].thick[i][j], 0);
		}
	}
	/*Delete empty Blocks*/
	for (int i_Block=0; i_Block<numBlocks; i_Block++) {
		float Block_volume=0;
		for (int i=0; i<Ny; i++) for (int j=0; j<Nx; j++)  Block_volume += Blocks[i_Block].thick[i][j];
		Block_volume *= dx*dy;
		if (Block_volume<1e4 && Blocks[i_Block].type != 'H' && Blocks[i_Block].type != 'I') {
			PRINT_DEBUG("will remove Block %d (type %c) out of %d", i_Block, Blocks[i_Block].type, numBlocks);
			Delete_Block(i_Block); i_Block--;
		}
	}
	PRINT_DEBUG("Repare finished");
	return(1);
}



int solveLESalmostdiagonal (double **A, double *b, float **x)
{
	register int i, j, numcorr, NDi=2*Ny, NDs=2*Ny;
	float *xcorr ;

	/*THIS ROUTINE CALLS THE ONES WHICH SOLVE THE EQUATION SYSTEM */
	xcorr = alloc_array(Nx*Ny);

	if ( TriangularizeAlmostDiagonalEquationSystem(A, b, Nx*Ny, NDs, NDi) ) {
			PRINT_ERROR("\aUNDETERMINED EQ. SYSTEM !!!"); getchar(); return 1;
	}

	SolveAlmostDiagonalTriangularEquationSystem(A, b, Nx*Ny, NDs, NDi, xcorr);

	for (i=0; i<Ny; i++) for (j=0; j<Nx; j++) {
			numcorr = Ny*j + i ;
			x[i][j] = xcorr[numcorr] ;
	}

	return 0;
}



float Sort_Matrix (float **matrix, struct GRIDNODE *orden, int Nx, int Ny)
{
	int	numorden, i, j, fil, col;
	BOOL	**switch_done;
	float	maxmatrix;

	switch_done = (BOOL **) calloc (Ny, sizeof(BOOL *));
	for (i=0; i<Ny; i++) switch_done[i] = (BOOL *) calloc (Nx, sizeof(BOOL));

	for (numorden=0; numorden < Nx*Ny; numorden++) {
		maxmatrix=-1e12;
		for (i=0; i<Ny; i++) for (j=0; j<Nx; j++) {
			if ((maxmatrix < matrix[i][j])  &&  (switch_done[i][j] == NO)) {
				maxmatrix = matrix[i][j];
				fil=i;	col=j;
			}
		}
		orden[numorden].row=fil;
		orden[numorden].col=col;
		switch_done[fil][col]=YES;
		if (verbose_level>=1) fprintf(stdout, "\b\b\b%2d%%", (int) 100*numorden/(Nx*Ny-1) );
	}

	for (i=0; i<Ny; i++) free(switch_done[i]);
	free(switch_done);

	return (matrix[orden[0].row][orden[0].col]) ;
}



float ReSort_Matrix (float **matrix, struct GRIDNODE *orden, int Nx, int Ny)
{
	register int 	numorden, i, j, aux1, aux2;

	/* SORTS MATRIX NODES GIVEN A PREVIOUS SORT */

	for (numorden=0; numorden < Nx*Ny-1; numorden++) {
		if ( matrix[orden[numorden].row][orden[numorden].col] < matrix[orden[numorden+1].row][orden[numorden+1].col] ) {
			aux1 = orden[numorden+1].row;  
			aux2 = orden[numorden+1].col;
			for (i=numorden; i>=0 ; i--) {
				if ( matrix[orden[i].row][orden[i].col] >= matrix[aux1][aux2] )  break;
				orden[i+1].row = orden[i].row;
				orden[i+1].col = orden[i].col;
			}
			orden[i+1].row = aux1;
			orden[i+1].col = aux2;
		}
	}
	/*Check*/
	for (numorden=0; numorden < Nx*Ny-1; numorden++) {
		if ( matrix[orden[numorden].row][orden[numorden].col] < matrix[orden[numorden+1].row][orden[numorden+1].col] ) {
			PRINT_ERROR("UNSORTED!");
		}
	}

	/*Returns the maximum value:*/
	return (matrix[orden[0].row][orden[0].col]) ;
}



float calculate_sea_level()
{
	/*
	  Calculates the sea level
	*/

	if (!water_load) return (0);

	/*Calculates sea level*/
	if (n_sea_level_input_points) {
	    int i;
	    for (i=0; i<n_sea_level_input_points; i++) {
		if (var_sea_level[i][0]>=Time) break;
    	    }
    	    if (i!=0 && i!=n_sea_level_input_points) {
		sea_level = 	( (Time-var_sea_level[i-1][0])*var_sea_level[i][1] + 
				  (var_sea_level[i][0]-Time)*var_sea_level[i-1][1] ) 
				     / (var_sea_level[i][0]-var_sea_level[i-1][0]);
	    }
 	    else {
		if (i==0)			 sea_level = var_sea_level[0][1];
		if (i==n_sea_level_input_points) sea_level = var_sea_level[n_sea_level_input_points-1][1];
	    }
	}
	else sea_level = 0;
	return(sea_level);
}



int calculate_water_load()
{
	/*
	  Calculates the load related to changes 
	  in the water column (sea and lakes).
	*/
	float	water_volume=0;

	calculate_topo(topo);

	if (!water_load) return(0);

	for (int i=0; i<Ny; i++) for (int j=0; j<Nx; j++) {
		int il;
		float Dq_water, h_water_now=0;
		if (hydro_model) {
		    if ((il=drainage[i][j].lake)) {
		    	/*sea lake already has its proper level defined*/
			h_water_now = MAX_2(0, Lake[drainage[i][j].lake].alt-topo[i][j]);
		    }
		}
		else {
		    h_water_now = MAX_2(0, sea_level-topo[i][j]);
		}
		Dq_water = (h_water_now-h_water[i][j]) * g * (denswater-densenv);
		h_water[i][j] = h_water_now;
		/*Don't load the initial water column*/
		if (Time>Timeini) Dq[i][j] += Dq_water;
		water_volume += h_water[i][j];		
	}

	if (n_sea_level_input_points) {
		PRINT_SUMLINE("sea_level: %8.1f m   sea_volume = %.1f km3", sea_level, water_volume*dx*Nx/(Nx-1)*Ny/(Ny-1)/1e9);
	}
	return(1);
}




