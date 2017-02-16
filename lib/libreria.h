/*Declaration of functions of libreria.c */
float 	*alloc_array		(int num_fil);
double	*alloc_array_dbl	(int num_fil);
float 	**alloc_matrix  	(int num_fil, int num_col);
double 	**alloc_matrix_dbl 	(int num_fil, int num_col);
int 	**alloc_matrix_int 	(int num_fil, int num_col);
int 	free_matrix 		(float **matrix, int num_fil);
int 	free_matrix_dbl 	(double **matrix, int num_fil);
int 	free_matrix_int		(int **matrix, int num_fil);

char 	*replace_word	(char *s, char *old, char *new);

float dist(
	float x0, 		/* x coordinate of first point of segment	*/
	float y0, 		/* y coordinate of first point of segment	*/
	float x1, 		/* x coordinate of ending point of segment	*/
	float y1);		/* y coordinate of ending point of segment	*/
float dist_xy_segment		(float x, float y, float X0, float Y0, float X1, float Y1);
float distVincenty(float lat1, float lon1, float lat2, float lon2);
int interpol2D (
	float 	**xarxa, 		/*matriz bidimensional de retorno*/
	int	Nx, int Ny, 		/*Dimensions of the return matrix*/
	float	xmin, float ymin, 	/*Geometric left-down boundary of matrix domain*/
	float	dx, float dy, 		/*Spatial increments or grid arms (constant)*/
	float 	**datinput, 		/*vector de 3 componentes con los datos a interpolar*/
	int 	numinputs, 		/*numeros de datos*/
	int 	mode_interp);		/*#1: peso = inverso de la distancia
					  #2: peso = inverso del cuadrado de la distancia
					  #3: toma el valor mas cercano
					  #4: altura constante dentro de un poligono
					  #5: binary skyline
					  #6: text skyline (1 float number per row)
					*/
float interpol_point_in_mesh (
	float **M,		/*Matrix of the mesh: M[i-y][j-x].*/
	int Nx, int Ny, 	/*Grid x and y knots.*/
	float xmin, float dx,	/**/
	float ymin, float dy,	/**/
	float x_point,		/*X of Point in wich to interpolate.*/
	float y_point);		/*Y of Point in wich to interpolate.*/

float project_xy_segment(
	float X, 		/* x coordinate of the point to be inspected	*/
	float Y, 		/* y     "    					*/
	float x0, 		/* x coordinate of first point of segment	*/
	float y0, 		/* y coordinate of first point of segment	*/
	float x1, 		/* x coordinate of ending point of segment	*/
	float y1,		/* y coordinate of ending point of segment	*/
	float *xp, 		/* x coordinate of profected point 		*/
	float *yp,		/* y coordinate of profected point 		*/
	float *len);		/* length in line to the projection 		*/
int outin (
	float x, 
	float y, 
	float *poligonX, 
	float *poligonY, 
	int number_of_points_in_poligon);
int readinterp2D (
	FILE *file, 
	float **z_matrix, 	/*returned matrix*/
	int mode_interp,      	/*
				#0: no interpolation needed ('Nx x Ny' file)
				#1: inverse distance
				#2: inv. square dist.
				#3: nearest given point
				#4: poligons
				#5: binary (short int) skyline
				#6: text skyline (1 float number per row)*/
	float z_default,	/*undefined nodes will have this value*/
	float xmin, float xmax, 
	float ymin, float ymax, 
	int Nx, int Ny); 
