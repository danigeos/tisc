/*
	GENERAL SUBRUTINES LIBRARY

		Daniel GARCIA-CASTELLANOS, since II-1995
*/


#include "types_n_defs.h"
#include "libreria.h"

extern float 	
	pi, 			/*PI number*/
	CGU, 			/*Constant of Universal Gravitation [m3·s-2·kg-1]*/
	Rearth;			/*Earth's radius*/

extern int	
	verbose_level;

extern BOOL	
	switch_geograph_coor;			/*1 if x-y are geographycal coordinates in decimal degrees*/




float *alloc_array (int n)
{
	/* ALLOCATES MEMORY FOR A FLOAT ARRAY */

	float 	*newarr;
	if (!n) return (NULL);
	if ((newarr = (float *) calloc(n, sizeof(float))) == NULL) {
		printf("\nERROR: alloc_array: There isn't memory enough for an array of %d knots of 16 bits.\n", n);
		exit(0);
	}
	return (newarr);
}


double *alloc_array_dbl (int n)
{
	/* ALLOCATES MEMORY FOR A DOUBLE ARRAY */

	double 	*newarr;
	if (!n) return (NULL);
	if ((newarr = (double *) calloc(n, sizeof(double))) == NULL) {
		printf("\nERROR: alloc_array_dbl: There isn't memory enough for an array of %d knots of 32 bits.\n", n);
		exit(0);
	}
	return (newarr);
}



float **alloc_matrix (int n_fil, int n_col)
{
	/* ALLOCATES MEMORY FOR A FLOAT MATRIX */

	float 	**new_matrix;
	if (!n_fil || !n_col) return (NULL);
	if ((new_matrix = (float **) calloc(n_fil , sizeof(float *))) == NULL) {
		printf("\nERROR: alloc_matrix: There isn't memory enough for a net of %dx%d knots of 16 bits.\n", n_col, n_fil);
		exit(0);
	}
	for (int i=0; i<n_fil; i++) {
		if ((new_matrix[i] = (float *) calloc(n_col , sizeof(float))) == NULL) {
			printf("\nERROR: There isn't memory enough for a net of %dx%d knots of 16 bits.\n", n_col, n_fil);
			exit(0);
		}
	}	
	return (new_matrix);
}


double **alloc_matrix_dbl (int n_fil, int n_col)
{
	/* ALLOCATES MEMORY FOR A DOUBLE MATRIX */

	double 	**new_matrix;
	if (!n_fil || !n_col) return (NULL);
	if ((new_matrix = (double **) calloc(n_fil , sizeof(double *))) == NULL) {
		printf("\nERROR: alloc_matrix_dbl: There isn't memory enough for a net of %dx%d knots of 32 bits.\n", n_col, n_fil);
		exit(0);
	}
	for (int i=0; i<n_fil; i++) {
		if ((new_matrix[i] = (double *) calloc(n_col , sizeof(double))) == NULL) {
			printf("\nERROR: There isn't memory enough for a net of %dx%d knots of 32 bits.\n", n_col, n_fil);
			exit(0);
		}

	}
	return (new_matrix);
}



int **alloc_matrix_int (int n_fil, int n_col)
{
	/* ALLOCATES MEMORY FOR A MATRIX OF INTEGERS*/

	int 	**new_matrix;
	if (!n_fil || !n_col) return (NULL);
	if ((new_matrix = (int **) calloc(n_fil , sizeof(int *))) == NULL) {
		printf("\nERROR: alloc_matrix_dbl: There isn't memory enough for a net of %dx%d knots of 32 bits.\n", n_col, n_fil);
		exit(0);
	}
	for (int i=0; i<n_fil; i++) {
		if ((new_matrix[i] = (int *) calloc(n_col , sizeof(int))) == NULL) {
			printf("\nERROR: There isn't memory enough for a net of %dx%d knots of 8? bits.\n", n_col, n_fil);
			exit(0);
		}

	}
	return (new_matrix);
}



float anombloq(register float x1, register float x2, float z1, float z2, float contrdens)
{
	register float 	anomaliaaux, x1cuad, x2cuad, z1cuad, z2cuad;
	//CALCULATES GRAVITY ANOMALY DUE TO A RECTANGLE 2D  (3D HORIZONTAL 
	//INFINITE PRISM). INTERNATIONAL SYSTEM UNITS USED.
	//Z1 and Z2 are positive vertical distances if the rectangle is under measure point.
	//Always z1<z2, even if they are negative (block over measure level). 
	//If z2<z1 then this routine giives the anomaly with changed sign.

	if (z1==0) 	z1=.001 ;
	if (z2==0) 	z2=.001 ;
	if (z1==z2) 	z2+=.001 ;	

	x1cuad=x1*x1; x2cuad=x2*x2; z1cuad=z1*z1; z2cuad=z2*z2;

	anomaliaaux = x2*log((x2cuad+z2cuad)/(x2cuad+z1cuad));
	anomaliaaux += x1*log((x1cuad+z1cuad)/(x1cuad+z2cuad));
	anomaliaaux += 2*z2*(atan(x2/z2)-atan(x1/z2));
	anomaliaaux += 2*z1*(atan(x1/z1)-atan(x2/z1));
	anomaliaaux *= CGU*contrdens;
	return (anomaliaaux);
}



float compaction(float phi0, float comp_depth, float z1, float z2) 
{
	//Returns the decrease in thickness by compaction or decompaction 
	//due to the overlying load. Always positive. Thickness changes 
	//when bringing the layer z2-z1 to depth 0 or viceversa. 
	//z2>z1>0 (positive downwards) 
	//z1 is the top
	float z;
	if (!comp_depth) return(0);
	if (z1<0) z1=0; 
	if (z2<z1) z2=z1+1.; 
	float dh=(z2-z1)/4; /*Temptative value for iteration*/
	for (int i=0; i<100; i++) {
		float dha=dh; 
		dh=phi0*comp_depth*(exp(-z2/comp_depth)-exp(-z1/comp_depth)+1-exp(-(dh+z2-z1)/comp_depth)); 
		if (fabs(dh-dha)<.1) break;
	} 
	dh = MIN_2(z2-z1, dh);
	return(dh); 
}



int diffusion_2D(float **Matrix, float **d_Matrix, int Nx, int Ny, float Kdiff, float dx, float dy, float dteros)
{
	int 	i, ii, j, jj, k;
	float	deriv2x, deriv2y, 
		mean_height, change, change_max;

	for (i=0; i<Ny; i++)  for (j=0; j<Nx; j++) {
		ii=i; jj=j;
		if (ii==0) 	ii=1;
		if (ii==Ny-1) 	ii=Ny-2;
		if (jj==0) 	jj=1;
		if (jj==Nx-1) 	jj=Nx-2;
		if (j>0 && j<Nx-1) 
			deriv2x = ( Matrix[ii][jj+1] - 2*Matrix[ii][jj] + Matrix[ii][jj-1] ) / (dx*dx);
		else
			deriv2x = 0;
		if (i>0 && i<Ny-1) 
			deriv2y = ( Matrix[ii-1][jj] - 2*Matrix[ii][jj] + Matrix[ii+1][jj] ) / (dy*dy);
		else
			deriv2y = 0;
		change = Kdiff * (deriv2x + deriv2y) * dteros ;
		/*Mean matrix value*/
		mean_height = ( Matrix[ii+1][jj] + Matrix[ii-1][jj] + Matrix[ii][jj+1] + Matrix[ii][jj-1] ) / 4 ;
		/*
		Maximum change is that making flat the matrix.
		ojo, esto no funciona si es un punto de ensilladura
		*/
		change_max = mean_height - Matrix[i][j] ;
		if (change >= 0) {
			change = MIN_2(change,  fabs(change_max));
		}
		else {
			change = MAX_2(change, -fabs(change_max));
		}
		d_Matrix[i][j] = change;
	}

	return (1);
}




float dist_2D(
	float x0, 		/* x coordinate of first point	*/
	float y0, 		/* y coordinate of first point	*/
	float x1, 		/* x coordinate of ending point	*/
	float y1)		/* y coordinate of ending point	*/
{
	float 	distance, dx, dy;

	/*
		RETURNS THE DISTANCE BETWEEN 2 POINTS X-Y, EVEN WHEN GIVEN IN 
		GEOGRAPHYCAL COORDINATES (geodesic)
	*/

	if (switch_geograph_coor) {
		/*returns meters*/
		/*Spherical Earth*/
		/*distance = Rearth * acos (sin(y0*pi/180)*sin(y1*pi/180) + cos(y0*pi/180)*cos(y1*pi/180)*cos((x1-x0)*pi/180));*/
		/*Ellipsoidal Earth*/
		distance = distVincenty(y0, x0, y1, x1);
	}
	else {
		dx = x1-x0;  
		dy = y1-y0;
		distance = sqrt(dx*dx + dy*dy);
	}

	/*printf("\n#x0=%.2f   y0=%.2f   x1=%.2f   y1=%.2f   d=%.2f\n", x0, y0, x1, y1, distance);*/

	return (distance);
}


float distVincenty(float lat1, float lon1, float lat2, float lon2) {
	/* LGPL license © 2002-2008 Chris Veness
	 * Calculate geodesic distance (in m) between two points specified by 
	 * latitude/longitude (in decimal degrees)
	 * using Vincenty inverse formula for ellipsoids
	 */
        char aux[40];
        double cosSqAlpha, sinSigma, cos2SigmaM, cosSigma, sigma, sinAlpha, sinLambda, cosLambda, C;
        double a = 6378137, b = 6356752.3142,  f = 1/298.257223563;  // WGS-84 ellipsiod
        double L = (lon2-lon1)*pi/180;
        double U1 = atan((1-f) * tan(lat1*pi/180));
        double U2 = atan((1-f) * tan(lat2*pi/180));
        double sinU1 = sin(U1), cosU1 = cos(U1);
        double sinU2 = sin(U2), cosU2 = cos(U2);

        double lambda = L, lambdaP = 2*pi;
        double iterLimit = 20;
        while (fabs(lambda-lambdaP) > 1e-12 && --iterLimit>0) {
          sinLambda = sin(lambda);
          cosLambda = cos(lambda);
          sinSigma = sqrt((cosU2*sinLambda) * (cosU2*sinLambda) +
            (cosU1*sinU2-sinU1*cosU2*cosLambda) * (cosU1*sinU2-sinU1*cosU2*cosLambda));
          if (sinSigma==0) return 0;  // co-incident points
          cosSigma = sinU1*sinU2 + cosU1*cosU2*cosLambda;
          sigma = atan2(sinSigma, cosSigma);
          sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
          cosSqAlpha = 1 - sinAlpha*sinAlpha;
          cos2SigmaM = cosSigma - 2*sinU1*sinU2/cosSqAlpha;
          sprintf(aux,"%f", cos2SigmaM);
          if (strstr(aux, "NaN")!=NULL) cos2SigmaM = 0;  // equatorial line: cosSqAlpha=0 (§6)
          C = f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha));
          lambdaP = lambda;
          lambda = L + (1-C) * f * sinAlpha *
            (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)));
        }
        if (iterLimit==0) return -9999e99;  // formula failed to converge

        double uSq = cosSqAlpha * (a*a - b*b) / (b*b);
        double A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)));
        double B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)));
        double deltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)-
          B/6*cos2SigmaM*(-3+4*sinSigma*sinSigma)*(-3+4*cos2SigmaM*cos2SigmaM)));
        double s = b*A*(sigma-deltaSigma);

        return s;
}



float dist_xy_line(
	float x, 		/* x coordinate of the point to be inspected	*/
	float y, 		/* y     "    					*/
	float *lineX, 		/* Polygon x coordinates (lineX[np])		*/
	float *lineY, 		/* Polygon y coordinates (lineY[np])		*/
	int   np)		/* Number of points in line.			*/
{
	int	i;
	float 	distance=1e25, this_distance, 
		x0, y0, xf, yf;
	/*
		RETURNS THE DISTANCE FROM X-Y POINT TO A LINE.
		BECAREFUL: returned distance is >0 at right, <0 at left of line
						Daniel G. C. 7-X-95
	*/

	for (i=0; i<np-1; i++){
		x0 = lineX[i];
		y0 = lineY[i];
		xf = lineX[i+1];
		yf = lineY[i+1];

		this_distance = dist_xy_segment(x, y, x0, y0, xf, yf);
		if (fabs(this_distance) < fabs(distance)) distance = this_distance ;
	}

	return ( distance );	
}



float dist_xy_pol(
	float x, 		/* x coordinate of the point to be inspected	*/
	float y, 		/* y     "    					*/
	float *lineX, 		/* Polygon x coordinates (lineX[np])*/
	float *lineY, 		/* Polygon y coordinates (lineY[np])*/
	int np)			/* Number of points in line.*/

{
	int	i;
	float 	distance=1e20, this_distance, 
		x0, y0, xf, yf;

	/*
	  RETURNS THE DISTANCE FROM X-Y POINT TO A POLYGON.
	  ATTENTION: returned distance is <0 negative! if the polygon is clockwise 
	  and >0 if counterclockwise.
					Daniel G. C. 7-X-95
	*/

	for (i=0; i<np; i++){
		x0 = lineX[i];
		y0 = lineY[i];

		if (i<(np-1)) {
			xf = lineX[i+1];
			yf = lineY[i+1];
		}	
		else {
			xf = lineX[0];
			yf = lineY[0];
		}

		this_distance = dist_xy_segment(x, y, x0, y0, xf, yf);
		if (fabs(this_distance) < fabs(distance)) distance = this_distance ;
	}

	return ( distance );	
}


float dist_xy_segment(
	float X, 		/* x coordinate of the point to be inspected	*/
	float Y, 		/* y     "    					*/
	float x0, 		/* x coordinate of first point of segment	*/
	float y0, 		/* y coordinate of first point of segment	*/
	float x1, 		/* x coordinate of ending point of segment	*/
	float y1)		/* y coordinate of ending point of segment	*/
{
	float 	dist_to_0, dist_to_1, length_segment, escy=1, escx=1,
		dx01, dy01, dxX0, dyY0, dxX1, dyY1, dx10, dy10, dx0X, dy0Y, dx1X, dy1Y, 
		cos_0to1_XYto1, cos_0to1_0toXY, cos_1to0_1toXY, 
		distance, dist_to_ends;

	/*
		RETURNS THE DISTANCE FROM A X-Y POINT TO A LINE SEGMENT
		distance is > 0 if the point is in the right side of the 
		segment (moving from 0 to 1).
		For the case of geographical coordinates, assumes that
		the earth is flat between the point and its projection, i. e., 
		that the distance is sqrt(Dx*Dx+Dy*Dy). This is zero if the 
		projection and the point are in the same parallel or meridian,
		but not in the other cases.
	*/

	if (switch_geograph_coor) {
		escy = pi/180 * Rearth;
		escx = escy * cos((2*Y+y0+y1)/4*pi/180);
	}
	dx01 = (x1-x0) * escx;
	dy01 = (y1-y0) * escy;
	dxX0 = (x0-X) * escx;
	dyY0 = (y0-Y) * escy;
	dxX1 = (x1-X) * escx;
	dyY1 = (y1-Y) * escy;
	dx10 = -dx01;
	dy10 = -dy01;
	dx0X = -dxX0;
	dy0Y = -dyY0;
	dx1X = -dxX1;
	dy1Y = -dyY1;

	dist_to_0  = dist_2D(X, Y, x0, y0);
	dist_to_1  = dist_2D(X, Y, x1, y1);
	length_segment = dist_2D(x0, y0, x1, y1);

	if (dist_to_0 == 0) return (0);
	if (dist_to_1 == 0) return (0);
	if (length_segment == 0) return (dist_to_0);

	cos_0to1_XYto1  = (dx01*dxX1 + dy01*dyY1) / dist_to_1 / length_segment;
	cos_0to1_0toXY  = (dx01*dx0X + dy01*dy0Y) / dist_to_0 / length_segment;
	cos_1to0_1toXY  = (dx10*dx1X + dy10*dy1Y) / dist_to_1 / length_segment;

	distance = dist_to_1 * sqrt(fabs(1 - cos_0to1_XYto1*cos_0to1_XYto1));

	if (dist_to_0 < dist_to_1) 	dist_to_ends = dist_to_0;
	else				dist_to_ends = dist_to_1;	

	if (cos_0to1_0toXY<=0 || cos_1to0_1toXY<=0) distance = dist_to_ends;
	if (dist_to_ends==0) distance=0;

	if (dx0X*dy01-dy0Y*dx01 < 0 ) distance *= -1;

	/*printf("\n#x0=%.2f  y0=%.2f    x1=%.2f  y1=%.2f  d=%.2f\n", x0,y0, x1,y1, distance);*/

	return (distance);
}




float evaluate_xy_points(
	FILE *file_true, 	/*File with well known points*/
	FILE *file_guess,	/*File with model points -- will be linearly interpolated*/
	char mode)		/*'W' if you want to weight with 3rd column of true file*/
{
	/*
		Routine to evaluate a set of x-y guessed points in respect to 
		another set of true x-y points (optionally weighted). 
		x-y points are read from two opened files. Returns the mean 
		error value in the true points. 
		Guess file is linearly interpoled at true points and then error 
		is computed.
		Repetition of true points leads to increment in its associated 
		weight. Weights for true points can be set in the third column.
	*/

#define NMAXPOINTS 20000

	int 	i, j, num_true, num_guess, num_fields, n_weight=0;
	float	a, b, w, 
		y_guess_interp, 
		error, error2, sum_weights, 
		x_guess[NMAXPOINTS], y_guess[NMAXPOINTS], 
		x_true[NMAXPOINTS], y_true[NMAXPOINTS], w_true[NMAXPOINTS];
	char	linea[MAXLENLINE], *lin ;

	for (num_true=n_weight=0; ; ) {
		num_fields = sscanf(lin=fgets(linea, MAXLENLINE-1, file_true), "%f %f %f", &a, &b, &w);
		if (num_fields >= 2) {
			x_true[num_true] = a;
			y_true[num_true] = b;
			if (num_fields >= 3) {
				w_true[num_true] = w;
				n_weight ++;
			}
			else {
				w_true[num_true] = 1;
			}
			num_true++;
		}
		if (lin==NULL) break;
	}

	if (mode != 'W') n_weight = 0;

	for (num_guess=0; ; ) {
		num_fields = sscanf(lin=fgets(linea, MAXLENLINE-1, file_guess), "%f %f", &a, &b);
		if (num_fields == 2) {
			x_guess[num_guess] = a;
			y_guess[num_guess] = b;
			num_guess++;
		}
		if (lin==NULL) break;
	}

	if (verbose_level>=1) {
		if (n_weight!=num_true && n_weight) {
			PRINT_WARNING("in evaluate_xy_points: Only part of the true points are weighted. Weight will be equal for all true points.\n");
		}
		if (n_weight==num_true && num_true) {
			fprintf(stderr, "\n\nWeight will be applied using third column.\n");
		}
		fprintf(stderr, "\nNum_True=%d, Num_Weight=%d, Num_Guess=%d", num_true, n_weight, num_guess);
		if (n_weight!=num_true && mode == 'W') PRINT_WARNING("in evaluate_xy_points: Not all weigths in true file using mode 'W'. Weight will be equal for all true points.\n");
	}
	if (n_weight!=num_true) {
		for (i=0; i<num_true; i++) {
			w_true[i] = 1;
		}
	}
	for (i=0, sum_weights=0; i<num_true; i++) {
		sum_weights += w_true[i];
	}
	PRINT_INFO("Sum_weights=%f", sum_weights);

	for (i=0, error=error2=0; i<num_true; i++) {
		for (j=0; j<num_guess; j++) {
			if (x_guess[j]>=x_true[i]) break;
		}
		if (j<=0) 		y_guess_interp = y_guess[0];
		if (j>=num_guess) 	y_guess_interp = y_guess[num_guess-1];
		if (j>0 && j<num_guess)	{
			y_guess_interp = y_guess[j-1] + 
				(x_true[i]-x_guess[j-1]) * 
				(y_guess[j]-y_guess[j-1]) /
				(x_guess[j]-x_guess[j-1]) ;
		}
		error  += fabs(y_guess_interp - y_true[i]) * w_true[i] / sum_weights ;
		error2 +=     (y_guess_interp - y_true[i]) * w_true[i] / sum_weights ;
	}

	/*Returns the mean error in the true points ('error')*/
	return(error/*+fabs(error2)*/);
}



int free_matrix (float **matrix, int n_fil)
{
	/* FREES THE MEMORY ALLOCATED FOR A FLOAT MATRIX */

	int 	i;

	for (i=0; i<n_fil; i++) {
		if (matrix[i]) free(matrix[i]);
		else PRINT_DEBUG("invalid address");
	}	
	if (matrix) free(matrix);
	else PRINT_DEBUG("invalid address (2)");

	return (1);
}



int free_matrix_dbl (double **matrix, int n_fil)
{
	/* FREES THE MEMORY ALLOCATED FOR A FLOAT MATRIX */


	for (int i=0; i<n_fil; i++) {
		free(matrix[i]);
	}	
	free(matrix);

	return (1);
}




int free_matrix_int (int **matrix, int n_fil)
{
	/* FREES THE MEMORY ALLOCATED FOR A FLOAT MATRIX */

	for (int i=0; i<n_fil; i++) {
		free(matrix[i]);
	}	
	free(matrix);

	return (1);
}




float geoidanompolyg(
	float *y_pol, 			/*horizontal coordinate of line points clock-wise sorted*/
	float *z_pol, 			/*vertical coordinate of line points clock-wise sorted. z>0 downwards*/
	int numpoints, 			/*Number of line points*/
	float ym, float zm, 		/*position of measurement*/
	float dens_contrast)		/*Density of the body*/
{
/*
		 CALCULATES GEOID HEIGHT ANOMALY DUE TO A PolygON 2D 
	(3D HORIZONTAL INFINITE POLIEDRUS). INTERNATIONAL SYSTEM UNITS USED.

	Polygon must be clock-wise sorted, otherwise anomaly sign is changed.
	z are positive vertical distances if the line is under measure point, 
	y is the horizontal distance.
	The geoid anomaly generated by a clockwise sorted body with positive
	density constrast is not necesarly positive, but it has a maximum at 
	the minimum distance (convexe shape).
	Theory and notation from Chapman, 1979; Ayala et al., 1996.
*/

	double 	N=0, a2, mi, bi, ci, di, mi2, ci2,
		ni=0, aux1, aux2, aux5, 
		zi, zh, yi, yh, 
		Dyi, Dzi, Dyh, Dzh, Dyi2, Dzi2, Dyh2, Dzh2,
		t1,  t2,  t3,  t4,  t5,  t6,  t7,  t8,  t9,  t10, 
		t11, t12, t13, t14, t15, t16, t17, t18, t19, t20 ;
	int	h, i;
	char 	aux[40];

	if (!dens_contrast) return (0);

	for (i=0; i<numpoints; i++) {
		/*h is the point previous to i */
		h = (i==0)? numpoints-1 : i-1;

		zi = z_pol[i];
		zh = z_pol[h];
		yi = y_pol[i];
		yh = y_pol[h];

		/*What to do if this or the previous point coincide with the measuring point*/
		if ((ym==yi && zm==zi) || (ym==yh && zm==zh)) zm -= 1;

		/*Auxiliar variables*/
		Dzi = zi - zm;	Dzi2 = Dzi*Dzi;
		Dzh = zh - zm;	Dzh2 = Dzh*Dzh;
		Dyi = yi - ym;	Dyi2 = Dyi*Dyi;	
		Dyh = yh - ym;	Dyh2 = Dyh*Dyh;	

		/*Geometric variables if the segment h-i is a line z = mi*y + bi
		*/
		mi = bi = ci = di = NO_DATA;
		if (yi!=yh) mi = (zi-zh) / (yi-yh);	 mi2 = mi*mi;
		bi = zi - mi * yi;
		a2 = mi * ym + bi - zm;
		if (zi!=zh) ci = (yi-yh) / (zi-zh);	 ci2 = ci*ci;
		di = -bi * ci - ym;


		/*General case: */
		if (zm != mi*ym+bi && ci != 0 && mi != 0 && ((yi-yh) != 0 ||
		(zi-zh) != 0)) {
			/*OJO*/
			aux1 = ( Dzi2 + SQUARE(Dzi*ci + ci*zm+di) ) / SQUARE(ci*zm+di);
			aux2 = ( Dzh2 + SQUARE(Dzh*ci + ci*zm+di) ) / SQUARE(ci*zm+di);

			/*Evaluation in yi, zi*/
			t1 = 	+ mi/2 * Dyi2 * log(Dyi2 + Dzi2);
			t2 =	- mi/2 * Dyi2 ;
			t3 =	+ mi2 * a2 / (1 + mi2) * Dyi;
			t4 =	- mi * a2*a2 * (mi2 - 1) / 2 / SQUARE(1+mi2) * log(Dyi2 + Dzi2);
			t5 =	- 2 * mi2 * a2*a2 / SQUARE(1+mi2) * atan(((1+mi2)*Dyi + mi*a2) / a2);
			t6 =	- mi * Dyi2 ;
			t7 =	+ Dzi2 * atan(Dyi/ Dzi);
			t8 =	+ (ci * zm + di) / (1 + ci2) * Dzi;
			t9 =	- ci * SQUARE(ci*zm+di) / SQUARE(1+ci2) * log(aux1);
			t10=	+ SQUARE(ci*zm+di) * (1-ci2) / SQUARE(1+ci2) * atan(Dyi/Dzi) ;

			/*Evaluation in yh, zh*/
			t11 =	+ mi/2 * Dyh2 * log(Dyh2 + Dzh2);
			t12 =	- mi/2 * Dyh2 ;
			t13 =	+ mi2 * a2 / (1 + mi2) * Dyh;
			t14 =	- mi * a2*a2 * (mi2 - 1) / 2 / SQUARE(1+mi2) * log(Dyh2 + Dzh2);
			t15 =	- 2 * mi2 * a2*a2 / SQUARE(1+mi2) * atan(((1+mi2)*Dyh + mi*a2) / a2);
			t16 =	- mi * Dyh2 ;
			t17 =	+ Dzh2 * atan(Dyh/ Dzh);
			t18 =	+ (ci * zm + di) / (1 + ci2) * Dzh;
			t19 =	- ci * SQUARE(ci*zm+di) / SQUARE(1+ci2) * log(aux2);
			t20 =	+ SQUARE(ci*zm+di) * (1-ci2) / SQUARE(1+ci2) * atan(Dyh/Dzh) ;

			ni = t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8 + t9 + t10 - t11 - t12 - t13 - t14 -t15 - t16 - t17 - t18 - t19 - t20;
		}
 
		/*Special cases: */
		if (zm == mi*ym+bi) {
			t1 =  Dzi2 / 2 / mi * log((1+ci2) * Dzi2) 
				- 3*Dzi2 / 2 / mi + Dzi2 * atan(1/ mi);
			t2 = Dzh2 / 2 / mi * log((1+ci2) * Dzh2) 
				- 3*Dzh2 / 2 / mi + Dzh2 * atan(1/ mi);
			ni = t1 - t2;
			if (Dzi2 == 0 || Dzh2 == 0) ni = 0;		
			/*OJO*/
		}

		if (ci == 0) {
			t1 =  Dyi * ( Dzi * log(Dyi2 + Dzi2) - 2*Dzi +
			2*fabs(Dyi)*atan(Dzi/fabs(Dyi)) )
				- 2*Dyi*zi + (Dyi2 + Dzi2) * atan(Dyi/Dzi) +
				Dyi*Dzi ;
			t2 =  Dyi * ( Dzh * log(Dyi2 + Dzh2) - 2*Dzh +
			2*fabs(Dyi)*atan(Dzh/fabs(Dyi)) )
				- 2*Dyi*zh + (Dyi2 + Dzh2) * atan(Dyi/Dzh) +
				Dyi*Dzh ;
			ni = t1 - t2;
			if (Dyi2 + Dzi2 == 0 || Dyi2 + Dzh2 == 0 || Dzi==0 ||
			Dzh==0 || Dyi==0) ni = 0;	/*OJO*/
		}

		if (mi == 0) {
			ni = 0;
		}

		if ((yi-yh) == 0 && (zi-zh) == 0) {
			ni = 0;
		}

		sprintf(aux,"%f", ni);
		if (strstr(aux, "NaN")!=NULL) {
			fprintf(stderr, "\n!#! "); ni=0;
			fprintf(stderr, "\n# Dzi=%9.3e Dzh=%9.3e Dyi=%9.3e Dyh=%9.3e aux1=%9.3e aux2=%9.3e", 
				Dzi,Dzh,Dyi,Dyh, aux1,aux2);
			fprintf(stderr, "\n> %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e",
				t1,t2,t3,t4,t5,t6,t7,t8,t9,t10);
			fprintf(stderr, "\n< %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e",
				t11,t12,t13,t14,t15,t16,t17,t18,t19,t20);
		} 

		N += ni ;
	}

	N *= - CGU * dens_contrast / 9.8;

	return ((float) N);
}



int WriteAlmostDiagonalMatrix (double **A, double *b, int rango, char *filename, int NDs, int NDi)
{
	int	charsize;
	FILE	*fich;

	charsize = (int) 600/rango;
	fich = fopen(filename, "wt");
	if (fich) {
		PRINT_INFO("Writing matrix in file '%s'.", filename); 
	}
	else {
		PRINT_WARNING("Unable to open matrix file '%s'.", filename); 
		return(0);
	}
	for (int i=0; i < rango; i++) {
		for(int j=0; j <= NDs+NDi; j++) {
			fprintf(fich, "\n%3d\t%3d\t%23.16e", i, j, A[i][j]); 
			/*fprintf(fich, "%d\t%d\t%d 0 1 2\t%c\n", i, j,
			charsize, 
				( j-i <= NDs  &&  i-j <= NDi)?
					((A[i][j-i+NDi])? 
						((A[i][j-i+NDi]>0)? 'P' : 'N' )
						: ',' )  :  '.'
			);*/
		}
		fprintf(fich, "\nTerm. Indep.: %3d\t%23.16e\n", i, b[i]);
	}
	fclose(fich);

	return(1);
}



float gravanompolyg(
	float *x_pol, 				/*x-y points of polygon clock-wise sorted*/
	float *z_pol, 				/*x-y points of polygon clock-wise sorted. z>0 downwards*/
	int numpoints, 				/*Number of polygon points*/
	float x_measure, float z_measure, 	/*position of measurement*/
	float dens_contrast)			/*Density of the body*/
{
/*
	CALCULATES GRAVITY ANOMALY VERTICAL COMPONENT AND GEOID ANOMALY DUE 
	TO A PolygON 2D (3D HORIZONTAL INFINITE POLIEDRUS). 
	International System units used.
	Polygon must be clock-wise sorted, otherwise anomaly sign is changed.
	z are positive vertical distances if the polygon is under measure point.
	Theory and notation from Talwani et al., 1959.
*/

	float 	vert_anom=0, ai, thetai, thetaj, phii, Zi, 
		Zi_aux, Zi_aux2, Zi_aux3, Zi_aux4, Zi_aux5, 
		zi, zj, xi, xj;
	register int	i, j;
	char aux[40];

	if (!dens_contrast) return (0);

	for (i=0; i<numpoints; i++) {
		j = (i==numpoints-1)? 0 : i+1;

		zi  =	z_pol[i] - z_measure;
		zj  =	z_pol[j] - z_measure;
		xi  =	x_pol[i] - x_measure;
		xj  =	x_pol[j] - x_measure;

		if (zi != zj) ai = xj + zj*(xj-xi)/(zi-zj);
		thetai =  atan2(zi, xi);
		thetaj =  atan2(zj, xj);
		phii   =  atan2((zj-zi), (xj-xi));

		/*General case:*/
		if (xi != 0 && xj != 0 && zi != zj && xi != xj &&
		fabs(fabs(thetai)-fabs(thetaj))>=1e-6 && (xi!=0 || zi!=0) &&
		(xj!=0 || zj!=0)) {
			Zi_aux2 = cos(thetai)  * (tan(thetai)  - tan(phii));
			Zi_aux3 = cos(thetaj)  * (tan(thetaj)  - tan(phii));
			Zi_aux  = ( Zi_aux2 ) / ( Zi_aux3 );
			Zi_aux4 = thetai - thetaj-2*pi*floor((thetai - thetaj)/2/pi);
			if (Zi_aux4 > pi) Zi_aux4 -= 2*pi ;
			Zi_aux5 = Zi_aux4 + tan(phii) * log(Zi_aux);

			Zi      = ai * sin(phii) * cos(phii) * (Zi_aux5) ;
		}

		/*Special cases:*/
		if (xi == 0) {  	/* Case A */
			Zi_aux = cos(thetaj)*(tan(thetaj)-tan(phii)) / sin(thetai);
			Zi_aux2 = tan(phii)*log(fabs(Zi_aux));
			Zi_aux3 = ai*sin(phii)*cos(phii);
			Zi_aux4 = thetaj-thetai-2*pi*floor((thetaj-thetai)/2/pi);
			if (Zi_aux4 > pi) Zi_aux4 -= 2*pi ;
			Zi = -Zi_aux3 * (Zi_aux4 + Zi_aux2); 
		}

		if (xj == 0)	{	/* Case B */
			Zi_aux = cos(thetai)*(tan(thetai)-tan(phii)) / sin(thetaj);
			Zi_aux2 = tan(phii)*log(fabs(Zi_aux));
			Zi_aux3 = ai*sin(phii)*cos(phii);
			Zi_aux4 = thetai-thetaj-2*pi*floor((thetai-thetaj)/2/pi);
			if (Zi_aux4 > pi) Zi_aux4 -= 2*pi ;
			Zi =  Zi_aux3 * (Zi_aux4 + Zi_aux2); 
		}

		if (zi == zj)		/* Case C */
			Zi = zi*(thetaj-thetai);
		if (xi == xj)  		/* Case D */ 	/*sin(phii)* ??*/
			Zi = xi*log(cos(thetai)/cos(thetaj));
		if (fabs(fabs(thetai)-fabs(thetaj))<1e-6) 
			Zi = 0;		/* Case E */
		if (xi == 0 && zi == 0) /* Case F */
			Zi = 0;
		if (xj == 0 && zj == 0) 
			Zi = 0;		/* Case G */

		sprintf(aux,"%f", Zi);
		if (strstr(aux, "NaN")!=NULL) { fprintf(stderr, "*"); Zi=1e10; }
		if (strstr(aux, "NaN")!=NULL) {
				fprintf(stderr, "\n>%d  xi=%8.3e zi=%8.3e  xn=%8.3e zn=%8.3e \t a=%.3e thi=%.7f thn=%.7f fii=%.7f", 
					i+1, xi, zi,   xj, zj,   ai, thetai,
					thetaj, phii);
				fprintf(stderr, "\n>\t%.3e %.3e %.3e %.3e  Zi=%.3e", 
					Zi_aux, Zi_aux2, Zi_aux3, Zi_aux4, Zi);
		}	

		vert_anom += Zi ;
	}

	vert_anom *= 2*CGU*dens_contrast;

	return (vert_anom);
}




int interpol2D (
	float **xarxa, 		/*matriz bidimensional de retorno*/
	int Nx, int Ny, 	/*Dimensions of the return matrix*/
	float xmin, float ymin, /*Geometric left-down boundary of matrix domain*/
	float dx, float dy, 	/*Spatial increments or grid arms (constant)*/
	float **datinput, 	/*vector de 3 componentes con los datos a interpolar*/
	int numinputs, 		/*numeros de datos*/
	int mode_interp)	/*#1: peso = inverso de la distancia
				  #2: peso = inverso del cuadrado de la distancia
				  #3: toma el valor mas cercano
				  #4: altura constante dentro de un polygono
				  #5: binary skyline
				  #6: text skyline (1 float number per row)
				  #7: same as 4 but nodes out of all polygons are interpolated with the distance to each polygon (no defauklt value assigned)
				*/
{
	/* INTERPOLACION 2D de los ficheros de entrada*/
	int 	i, j, k, 
		closestdatinputnum, 
		npols, nmax_pols=200, *n_points_pol;
	float 	*interp_submode;
	float 	dist2,
		dist=0, 
		dist_old, 
		*height, 
		*peso, 
		**Xpol, 
		**Ypol, 
		sumapesos, 
		DX,DY, 
		x, y, 
		mindistance, distance ;

	if (!numinputs) return (0);
	n_points_pol = (int *) calloc(nmax_pols, sizeof(int));
	interp_submode = (float *) calloc(nmax_pols, sizeof(float));
	height = alloc_array(nmax_pols);
	peso = alloc_array(numinputs);
	Xpol = alloc_matrix(nmax_pols, numinputs);
	Ypol = alloc_matrix(nmax_pols, numinputs);

	PRINT_DEBUG("Interpoling from %d points in mode %d.", numinputs, mode_interp);
	switch (mode_interp) {
	  case 1:
	  case 2:
		for (i=0; i<Ny; i++) for (j=0; j<Nx; j++) {
			xarxa[i][j]=0;	sumapesos=0;
			for (k=0; k<numinputs; k++) {
				DX = datinput[k][0]-j*dx-xmin ;
				DY = datinput[k][1]-(Ny-1-i)*dy - ymin ;
				dist2 = SQUARE(DX) + SQUARE(DY) ;
				if (mode_interp==1) 	peso[k] = 1 / sqrt((dist2+1e-5)) ;
				if (mode_interp==2) 	peso[k] = 1 / (dist2+1e-5) ;
				sumapesos += peso[k];
			}
			for (k=0; k<numinputs; k++) {
				xarxa[i][j] += datinput[k][2] * peso[k];
			}
			xarxa[i][j] /= sumapesos;
		}
		break;
	  case 3:
		for (i=0; i<Ny; i++)  for (j=0; j<Nx; j++) {
			mindistance = 1e20;
			for (k=0; k<numinputs; k++) {
				DX = datinput[k][0]-j*dx-xmin;
				DY = datinput[k][1]-(Ny-1-i)*dy - ymin;
				distance = sqrt(DX*DX+DY*DY);
				if (mindistance > distance) {
					mindistance = distance;
					closestdatinputnum = k;
				}
			}
			xarxa[i][j] = datinput[closestdatinputnum][2] ;
		}
		break;
	  case 4:
	  case 7:
	  case 9:
	  case 10:
		/*Defines the polygons in variables*/
		for (k=0,npols=0; k<numinputs; k++) {
			if (datinput[k][2]!=NO_DATA) {
				if (datinput[k][2]!=111111 && datinput[k][2]!=222222) {
					npols++;
					height[npols-1] = datinput[k][2];
				}
				else {
					interp_submode[npols-1] = datinput[k][2];
				}
				
			}
			if (!npols) PRINT_ERROR("interpoling file in mode %d: npols=0.", mode_interp);
			Xpol[npols-1][n_points_pol[npols-1]]=datinput[k][0];
			Ypol[npols-1][n_points_pol[npols-1]]=datinput[k][1];
			n_points_pol[npols-1]++;
		}
		/*
		Interpolates: with 111111 signal, if a node is inside of the 
		i-th polygon but outside the (i+1)th polygon then weights with 
		the distance to each one.
		*/
		for (i=0; i<Ny; i++)  for (j=0; j<Nx; j++) {
			float interp_z=SIGNAL;  /*initializes this variable at every loop*/
			BOOL switch_in;
			x = xmin + dx*j;
			y = ymin + dy*(Ny-i-1);
			switch_in = NO;
			for (k=0; k<npols; k++) {
			    BOOL switch_in_old;
			    switch_in_old = switch_in;
			    dist_old = dist;
			    switch_in = outin(x, y, Xpol[k], Ypol[k], n_points_pol[k]);
			    if (mode_interp == 4 || mode_interp == 7)  dist = fabs(dist_xy_pol (x, y, Xpol[k], Ypol[k], n_points_pol[k]));
			    if (mode_interp == 9 || mode_interp == 10) dist = fabs(dist_xy_line(x, y, Xpol[k], Ypol[k], n_points_pol[k]));
			    if (mode_interp == 9 || mode_interp == 10) switch_in = 1;
			    if (switch_in)
			    	interp_z = height[k];
			    else 
			    	if (switch_in_old)
				    if (interp_submode[k-1] == 111111) {
					if (!dist || !dist_old) {
					  if (!dist)	 interp_z = height[k];
					  if (!dist_old) interp_z = height[k-1];
					}
					else {
					  interp_z =
					      (height[k]/dist + height[k-1]/dist_old) /
						  (1/dist+1/dist_old);
					}
				    }
			}
			if (interp_z != SIGNAL) xarxa[i][j] = interp_z;
			else if (mode_interp == 7 || mode_interp == 10) {
				float divisor = 0, dist, thick=0;
				for (k=0; k<npols; k++) {
				    if (mode_interp == 7)  dist = fabs(dist_xy_pol (x, y, Xpol[k], Ypol[k], n_points_pol[k]));
				    if (mode_interp == 10) dist = fabs(dist_xy_line(x, y, Xpol[k], Ypol[k], n_points_pol[k]));
				    if (!dist) {thick = height[k];  divisor = 1;  break;}
				    thick += height[k]/dist;
				    divisor += 1/dist;
				}
				interp_z = thick/divisor;
				xarxa[i][j] = interp_z;
			}
		}
		break;
	  default:
		PRINT_ERROR("Interpolation mode '%d' does not exist.\n", mode_interp);
		exit(0);
	}

	free(n_points_pol);
	free_matrix(Xpol, nmax_pols);
	free_matrix(Ypol, nmax_pols);
	free(height);
	free(peso);
	free(interp_submode);

	return(1);
}



float interpol_point_in_mesh (
	float **M,		/*Matrix of the mesh: M[i-y][j-x].*/
	int Nx, int Ny, 	/*Grid x and y knots.*/
	float xmin, float dx,	/*In TISC dx = (xmax-xmin)/(Nx-1)*/
	float ymin, float dy,	/*In TISC dy = (ymax-ymin)/(Ny-1)*/
	float x_point,		/*X of Point in which to interpolate.*/
	float y_point)		/*Y of Point in which to interpolate.*/
{
	/*INTERPOLATES THE Z VALUE OF A POINT WITH A GIVEN x-y-z MESH.*/	

	int	i1, j1, i2, j2, i3, j3, i4, j4 , i0,j0;
	float	result, ymax=ymin+(Ny-1)*dy, 
		d1, d2, d3, d4;

	/*i0 = (int) floor((ymax - y_point) / dy +.5);
	j0 = (int) floor((x_point - xmin) / dx +.5);
	DOMAIN_LIMIT(i0,j0);
	return(M[i0][j0]);*/

	/*upper-left node*/
	i1 = (int) floor((ymax - y_point) / dy );
	j1 = (int) floor((x_point - xmin) / dx );
	DOMAIN_LIMIT(i1,j1);
	d1 = sqrt(SQUARE(xmin+j1*dx-x_point) + SQUARE(ymin+(Ny-1-i1)*dy-y_point));

	/*upper-right node (except if in right border)*/
	i2 = i1 ;		
	j2 = j1+ ((j1<Nx-1) ? 1 : -1) ;
	d2 = sqrt(SQUARE(xmin+j2*dx-x_point) + SQUARE(ymin+(Ny-1-i2)*dy-y_point));

	/*lower-left node (except if in lower border)*/
	i3 = i1+ ((i1<Ny-1) ? 1 : -1) ;	
	j3 = j1 ;
	d3 = sqrt(SQUARE(xmin+j3*dx-x_point) + SQUARE(ymin+(Ny-1-i3)*dy-y_point));

	/*lower-right node*/
	i4 = i3 ;		
	j4 = j2 ;
	d4 = sqrt(SQUARE(xmin+j4*dx-x_point) + SQUARE(ymin+(Ny-1-i4)*dy-y_point));

	result =  ( 
		M[i1][j1] * d2*d3*d4 + 
		M[i2][j2] * d1*d3*d4 + 
		M[i3][j3] * d1*d2*d4 + 
		M[i4][j4] * d1*d2*d3 ) 
		/ (d2*d3*d4 + d1*d3*d4 + d1*d2*d4 + d1*d2*d3) ;

/*printf("\n ######## %.1f,%.1f: %f", x_point, y_point, result);*/
	return (result);
}


float interpol_in_xy_data (
	float *x_array,		/*x position of y data*/
	float *y_array,		/*y data*/
	int   n_x, 		/*number of x-y data*/
	float x)  		/*x in wich to interpolate*/
{
	/*  
	    RETURNS THE INTERPOLATED y IN AN x VALUE 
		GIVEN A SORTED SET OF x-y DATA
	*/
	int 	i;
	float 	y;

	if (x<=x_array[0])  	return(y_array[0]);
	if (x>=x_array[n_x-1]) 	return(y_array[n_x-1]);

	for (i=1; i<n_x; i++) if (x<=x_array[i]) break;

	y = (y_array[i-1]*(x_array[i]-x)+y_array[i]*(x-x_array[i-1]))/
			(x_array[i]-x_array[i-1]);

	return (y);
}




int outin(
	float x, 		/* x coordinate of the point to be inspected	*/
	float y, 		/* y     "    					*/
	float *polygonX, 	/* Polygon x coordinates (polygonX[np])*/
	float *polygonY, 	/* Polygon y coordinates (polygonY[np])*/
	int np)			/*Number of points in polygon.*/

{
	/*
		RETURNS 0 IF x,y IS OUTSIDE OF polygon, 1 IF INSIDE, -1 IF 
		BELONGING TO IT.
		'Inside' includes any finite closed region of 2D cartesian space, 
		(i.e., if the polygon cuts itself n times then n+1 closed 
		regions are considered).
		The way this routine works can be understood using the analogy
		of Gauss Theorem for the electric field (in 2D) of a point
		load: total flux along a  closed polygon including the load is
		non-zero; whereas flux  along a polygon which does not include
		the load is zero.
		I made this routine based on the magnetic field analogy using
		the definition of circulation: circulation around a current is 
		non-zero. 
					D. Garcia-Castellanos.  After 15-III-95
	*/

	int	i;
	float 	circulation=0, partial_circ, L, L2, D, D2, E, E2, x0, y0, xf, yf, 
		b, K, sq, disc, pol_length=0;

	for (i=0; i<np; i++){
		x0 = polygonX[i];
		y0 = polygonY[i];
		if (i<(np-1)) xf = polygonX[i+1];  else xf = polygonX[0];
		if (i<(np-1)) yf = polygonY[i+1];  else yf = polygonY[0];

		D2 = (x0-x )*(x0-x ) + (y0-y )*(y0-y);
		E2 = (xf-x )*(xf-x ) + (yf-y )*(yf-y);	
		L2 = (xf-x0)*(xf-x0) + (yf-y0)*(yf-y0);
		L = sqrt(L2);					/*Length of this polygon side.*/
		D = sqrt(D2);					/*Distance from the point to the first point of the side*/
		E = sqrt(E2);					/*Distance from the point to the secnd point of the side*/
		pol_length += L;
		K = ((x0-x)*(yf-y0) - (y0-y)*(xf-x0)) / L;	/*(x0-x,y0-y) x (xf-x0,yf-y0) / L */
		b = ((x0-x)*(xf-x0) + (y0-y)*(yf-y0)) * 2 / L;	/*(x0-x,y0-y) * (xf-x0,yf-y0) * 2 / L */
		disc = b*b-4*D2;
		sq = sqrt(fabs(disc));

		/*If the point coincides with the point of polygon:*/
		if (D == 0) return(-1);
		if (E == 0) return(-1);
		if (fabs(D+E-L) < L/1e6) return(-1);

		if (disc > 0) partial_circ = 
				 K/sq * (log((2*L+b-sq)/(2*L+b+sq)) - log((b-sq)/(b+sq)));
		if (disc < 0) partial_circ = 
				 K/sq * (atan((2*L+b)/sq) - atan(b/sq));
		/*If x,y is alineated with this polygon side:*/
		if (disc==0)  partial_circ =  + 1/D - 1/E;  
		/*If a point is repeated in the polygon:*/
		if (L == 0)   partial_circ = 0;

		circulation += partial_circ;
	}

	/*printf("\nLength of polygon = %.2f\nCirculation = %.5e\n", pol_length, circulation);*/

	/*Circulation results 0 outside, pi inside*/
	return ( (fabs(circulation)>1)? 1:0 );	
}



float polygon_area(
	float *X, 	/* Polygon x coordinates (polygonX[np])*/
	float *Y, 	/* Polygon y coordinates (polygonY[np])*/
	int np)		/* Number of points in polygon.*/

{
	/*
		RETURNS THE AREA OF A PolygON. EVEN NON CONVEX POLYGONS.
		Positive area if the polygon is given clockwise; negative 
		otherwise.
		You dont need to close the polygon.
			Daniel G. C. 25-VIII-98
	*/

	int	i;
	float 	area=0;

	for (i=0; i<np-1; i++){
		area += (Y[i]+Y[i+1]) * (X[i+1]-X[i]) / 2;
	}
	area += (Y[np-1]+Y[0]) * (X[0]-X[np-1]) / 2;

	return (area);	
}




double prism_vert_grav (
	float x1, float x2,	/*X horizontal distance from measure to prism x-sides in m [x1<x2]*/
	float y1, float y2,	/*Y horizontal distance from measure to prism y-sides in m [y1<y2]*/
	float z1, float z2, /*Z vertical   distance from measure to prism top and bottom faces, respectively, in m [z positive downwards, z1<z2]*/
	float dens) 		/*density of prism*/
{
	double vert_anom=0;

	/*
		VERTICAL GRAVITATIONAL ATTRACTION PRODUCED BY A VERTICAL parallelepiped of constant density

		x2>x1 ; y2>y1 ; z2>z1
		x2 = x_right_prism - x_measure
		z1 = z_top_prism - z_measure
		z is positive downwards
		Returns gravity in m/s2; a positive value for downwards attraction, i.e., for bodies of positive density below the measurement (z2>z1>0)
	*/

	if (!x1) x1=1e-2;
	if (!x2) x2=1e-2;
	if (!y1) y1=1e-2;
	if (!y2) y2=1e-2;

	if (!(x2-x1)) return(0);
	if (!(y2-y1)) return(0);
	if (!(z2-z1)) return(0);
	if ((x2-x1)<0) {PRINT_ERROR("prism_vert_grav: x1>x2."); return(0);}
	if ((y2-y1)<0) {PRINT_ERROR("prism_vert_grav: y1>y2."); return(0);}
	if ((z2-z1)<0) {PRINT_ERROR("prism_vert_grav: z1>z2. %.2e > %.2e", z1, z2); return(0);}
	if (!dens)    return(0);

#define VERTANOMCONTRIB(s,x,y,z) {\
	double R;\
	R=sqrt(x*x+y*y+z*z);\
	vert_anom += s * ( x*log(y+R) + y*log(x+R) + z*atan2(z*R,x*y) );\
}

	VERTANOMCONTRIB(+1, x2,y2,z2);
	VERTANOMCONTRIB(-1, x2,y2,z1);
	VERTANOMCONTRIB(-1, x2,y1,z2);
	VERTANOMCONTRIB(+1, x2,y1,z1);
	VERTANOMCONTRIB(-1, x1,y2,z2);
	VERTANOMCONTRIB(+1, x1,y2,z1);
	VERTANOMCONTRIB(+1, x1,y1,z2);
	VERTANOMCONTRIB(-1, x1,y1,z1);

	vert_anom *= - CGU * dens;

	if (vert_anom<0) {
		//Negative values may result from bodies above the measure!
		//PRINT_WARNING("Downwards gravity attraction is negative! x1/x2/y1/y2/z1/z2= %.2f/%.2f/%.2f/%.2f/%.3f/%.3f\tdens=%.1f\tanom=%.2e\nSet to 0.", x1, x2, y1, y2, z1, z2, dens, vert_anom);
		//vert_anom=0;
	}

	return (vert_anom);
}



float project_xy_line(
	float x, 		/* x coordinate of the point to be inspected	*/
	float y, 		/* y     "    					*/
	float *lineX, 		/* line x coordinates (lineX[np])		*/
	float *lineY, 		/* line y coordinates (lineY[np])		*/
	int np,			/* Number of points in line.			*/
	float *xp,		/* Return value of projected point: x		*/
	float *yp,		/* Return value of projected point: y		*/
	float *len)		/* Length in line to the projection 		*/

{
	int	i, i_mindist;
	float 	distance=1e20, this_distance, len_seg, 
		x0, y0, xf, yf;


	/*
		RETURNS THE PROJECTION OF AN X-Y POINT ON A LINE.
						Daniel G.-C. 7-X-95
	*/

	for (i=0; i<np-1; i++){
		x0 = lineX[i];
		y0 = lineY[i];
		xf = lineX[i+1];
		yf = lineY[i+1];

		this_distance = dist_xy_segment(x, y, x0, y0, xf, yf);
		if (fabs(this_distance) < fabs(distance)) {
			distance = this_distance ;
			i_mindist = i ;
		}
	}

	for (i=0, *len=0; i<i_mindist; i++){
		x0 = lineX[i];
		y0 = lineY[i];
		xf = lineX[i+1];
		yf = lineY[i+1];

		*len += dist_2D(x0, y0, xf, yf);
	}

	x0 = lineX[i_mindist];
	y0 = lineY[i_mindist];
	xf = lineX[i_mindist+1];
	yf = lineY[i_mindist+1];

	project_xy_segment (x, y, x0, y0, xf, yf, xp, yp, &len_seg);

	
	/*printf("\n%d %.2f  %.2f %.2f   %.2f \n", i_mindist, distance, *xp, *yp, len_seg);*/
	

	*len += len_seg;

	return ( distance );	
}


float project_xy_segment(
	float X, 		/* x coordinate of the point to be inspected	*/
	float Y, 		/* y     "    					*/
	float x0, 		/* x coordinate of first point of segment	*/
	float y0, 		/* y coordinate of first point of segment	*/
	float x1, 		/* x coordinate of ending point of segment	*/
	float y1,		/* y coordinate of ending point of segment	*/
	float *xp, 		/* x coordinate of profected point 		*/
	float *yp,		/* y coordinate of profected point 		*/
	float *len_seg)		/* length in line to the projection 		*/
{
	float 	dist_to_0, dist_to_1, length_segment, escy=1, escx=1,
		dx01, dy01, dxX0, dyY0, dxX1, dyY1, dx10, dy10, dx0X, dy0Y, dx1X, dy1Y, 
		cos_0to1_XYto1, cos_0to1_0toXY, cos_1to0_1toXY, 
		distance, dist_to_ends;

	/*RETURNS THE PROJECTION OF AN X-Y POINT TO A LINE SEGMENT (x0,y0)-(x1,y1)*/

	if (switch_geograph_coor) {
		escy = pi/180 * Rearth;
		escx = escy * cos(pi/180*(y0+y1)/2);
	}
	dx01 = (x1-x0) * escx;
	dy01 = (y1-y0) * escy;
	dxX0 = (x0-X) * escx;
	dyY0 = (y0-Y) * escy;
	dxX1 = (x1-X) * escx;
	dyY1 = (y1-Y) * escy;
	dx10 = -dx01;
	dy10 = -dy01;
	dx0X = -dxX0;
	dy0Y = -dyY0;
	dx1X = -dxX1;
	dy1Y = -dyY1;

	dist_to_0  = dist_2D(X, Y, x0, y0);
	dist_to_1  = dist_2D(X, Y, x1, y1);
	length_segment = dist_2D(x0, y0, x1, y1);

	if (dist_to_0 == 0)      {*xp=x0; *yp=y0; *len_seg=0;              return (0);}
	if (dist_to_1 == 0)      {*xp=x1; *yp=y1; *len_seg=length_segment; return (0);}
	if (length_segment == 0) {*xp=x0; *yp=y0; *len_seg=length_segment; return (dist_to_0);}

	cos_0to1_XYto1  = (dx01*dxX1 + dy01*dyY1) / dist_to_1 / length_segment;
	cos_0to1_0toXY  = (dx01*dx0X + dy01*dy0Y) / dist_to_0 / length_segment;
	cos_1to0_1toXY  = (dx10*dx1X + dy10*dy1Y) / dist_to_1 / length_segment;

	distance = dist_to_1 * sqrt(fabs(1 - SQUARE(cos_0to1_XYto1)));
	*len_seg = dist_to_0 * cos_0to1_0toXY;
	if (!dist_to_0) *len_seg=0;

	if (dist_to_0 < dist_to_1) 	dist_to_ends = dist_to_0;
	else				dist_to_ends = dist_to_1;	

	if (cos_0to1_0toXY<=0 || cos_1to0_1toXY<=0) distance = dist_to_ends;
	if (dist_to_ends==0) distance=0;
	if (cos_0to1_0toXY<=0) *len_seg=0;
	if (cos_1to0_1toXY<=0) *len_seg=length_segment;

	*xp = x0 + dist_to_0 * cos_0to1_0toXY * (x1-x0) / length_segment;
	*yp = y0 + dist_to_0 * cos_0to1_0toXY * (y1-y0) / length_segment;

	if (cos_0to1_0toXY<=0 || dist_to_0 == 0) { *xp = x0; *yp = y0; }
	if (cos_1to0_1toXY<=0 || dist_to_1 == 0) { *xp = x1; *yp = y1; }

	
	/*fprintf(stderr, "\nX =%8.2f  Y =%8.2f"
			"\nx0=%8.2f  y0=%8.2f \tx1=%8.2f  y1=%8.2f"
			"\nxp=%8.2f  yp=%8.2f \td=%8.2f l=%8.2f  \td1=%8.2f  d2=%8.2f", 
		X, Y, x0, y0, x1, y1, *xp, *yp, distance, *len_seg, dist_to_0, dist_to_1);
	*/

	return (distance);
}




int readinterp2D (
	FILE *file, 
	float **z_matrix, 	/*returned matrix*/
	int mode_interp,      	/*see interpol2D for available modes*/
	float z_default,	/*undefined nodes will have this value*/
	float xmin, float xmax, 
	float ymin, float ymax, 
	int Nx, int Ny) 
{
	/*LOADS AND INTERPOLATES THE INPUT FILES*/

	int 	i, j, num_fields, n_input_points, 
		nmax_input_points, **n_data_this_cell, n_not_used=0;
	short int z_binary;
	float 	x=0, y=0, z, dx=(xmax-xmin)/(Nx-1), dy=(ymax-ymin)/(Ny-1),
		**input_points, add_random=0;
	char	auxstr[MAXLENLINE], *lin ;

	{
	    int nlines=0, nread;
	    char str1[MAXLENLINE], str2[MAXLENLINE], line[MAXLENLINE+200], *lineptr;
#ifndef NMAXHEADERLINES
#define NMAXHEADERLINES 100
#endif
#ifndef Match_Param_int
#define Match_Param_int(x, y)  if (!strcasecmp(str1, x)) {y=atoi(str2); if (verbose_level>=3) fprintf(stdout, "\n"x"\t%d", y);}
#endif
#ifndef Match_Param_flt
#define Match_Param_flt(x, y)  if (!strcasecmp(str1, x)) {y=atof(str2);; if (verbose_level>=3) fprintf(stdout, "\n"x"\t%f", y);}
#endif
#ifndef Match_Param_char
#define Match_Param_char(x, y) if (!strcasecmp(str1, x)) {strcpy(y,str2); if (verbose_level>=3) fprintf(stdout, "\n"x"\t%s", y);}
#endif
	    rewind(file); 
 	    while ((lineptr=fgets(line, MAXLENLINE+200-1, file)) != NULL && nlines<NMAXHEADERLINES) {
 	    	nlines++; nread=sscanf(lineptr, "%s %s", str1, str2);
 	    	if (nread == 2) {
	    	    Match_Param_int ( "mode_interp",	    mode_interp )
	    	    Match_Param_flt ( "add_random", 	    add_random )
	    	    Match_Param_flt ( "z_default",	    z_default )
		}
	    }
	    rewind(file);
	}

	for (i=0; i<Ny; i++) for (j=0; j<Nx; j++) z_matrix[i][j] = z_default;

	nmax_input_points = MAX_2(3*Nx*Ny + 3, 5000); /*mode 3 needs more than Nx*Ny*/
	input_points = 	alloc_matrix(nmax_input_points, 3);

	n_input_points = 0;
	switch (mode_interp) {
	  case 1:
	  case 2:
	  case 3:
		for (i=0; ; i++) {
			while (1) {
				lin=fgets(auxstr, MAXLENLINE-1, file);
				if (lin==NULL) {num_fields=EOF; break;}
				num_fields = sscanf(lin, "%f %f %f", &x, &y, &z);
				if (num_fields>=3) break;
			}
			if (num_fields == EOF) break;
			input_points[i][0]=x;
			input_points[i][1]=y;
			input_points[i][2]=z;
		}
		n_input_points = i;
		PRINT_WARNING("Nx=%d, Ny=%d, xmin=%.1f, ymin=%.1f, dx=%.1f, dy=%.1f, mode_interp=%d", Nx, Ny, xmin, ymin, dx, dy, mode_interp);
		interpol2D(z_matrix, Nx, Ny, xmin, ymin, dx, dy, input_points, n_input_points, mode_interp);
		break;
	  case 4:
	  case 7:
	  case 9:
	  case 10:
		input_points[0][2]=NO_DATA;
		input_points[1][2]=NO_DATA;
		for (i=0; ;) {
			input_points[i+2][2]=NO_DATA;
			while (1) {
				lin=fgets(auxstr, MAXLENLINE-1, file);
				if (lin==NULL) {num_fields=EOF; break;}
				num_fields = sscanf(lin, "%f %f", &x, &y);
				if (num_fields>=1) break;
			} 
			if (num_fields == EOF) break;
			if (num_fields == 1) {
				int num_fields_height;
				char signals_word[20]="";
				num_fields_height = sscanf(lin, "%f %s", &x, signals_word);
				input_points[i][2]=x;
				input_points[i+1][2]= 111111;	/*default is 'n'*/
				if (num_fields_height>=2) {
					if (strchr(signals_word, 'n')) input_points[i+1][2]= 111111;
					if (strchr(signals_word, 'c')) input_points[i+1][2]= 222222;
				}
			}
			if (num_fields == 2) {
				input_points[i][0]=x;
				input_points[i][1]=y;
				i++;
			}
		}
		n_input_points = i;
		interpol2D(z_matrix, Nx, Ny, xmin, ymin, dx, dy, input_points, n_input_points, mode_interp);
		break;
	  case 0:
		for (i=0; i<Ny; i++) for (j=0; j<Nx; j++) {
			while (1) {
				lin=fgets(auxstr, MAXLENLINE-1, file);
				if (lin==NULL) {PRINT_ERROR("not Nx*Ny points in file."); exit(0);}
				if ((num_fields=sscanf(lin, "%f %f %f", &x, &y, &z)) >= 3 ) break;
			}
			z_matrix[i][j] = z;
		}
		n_input_points = Nx*Ny;
		break;
	  case 5:
		for (i=0; i<Ny; i++) for (j=0; j<Nx; j++) {
			if (!fread(&z_binary, sizeof(short int), 1, file)) {
				PRINT_ERROR("not Nx*Ny points in file."); 
				exit(0);
			}
/*			input_points[i*Nx+j][2] = z_binary;
*/			z_matrix[i][j] = z_binary;
		}
		n_input_points = Nx*Ny;
		break;
	  case 6:
		for (i=0; i<Ny; i++) for (j=0; j<Nx; j++) {
			while (1) {
				lin=fgets(auxstr, MAXLENLINE-1, file);
				if (lin==NULL) {PRINT_ERROR("not Nx*Ny values in the input file to be interpolated in mode %d.", mode_interp); exit(0);}
				if ((num_fields=sscanf(lin, "%f", &z)) >= 1 ) break;
			}
			z_matrix[i][j] = z;
		}
		n_input_points = Nx*Ny;
		break;
	  case 8:
		n_data_this_cell = calloc(Ny, sizeof(int *));
		for (i=0; i<Ny; i++) n_data_this_cell[i] = calloc(Nx, sizeof(int));
		while (1) {
			lin=fgets(auxstr, MAXLENLINE-1, file);
			if (lin==NULL) {break;}
			if ((num_fields=sscanf(lin, "%f %f %f", &x, &y, &z)) >= 3) {
			    j = floor((x-xmin)/dx + .5);
			    i = Ny - 1 - floor((y-ymin)/dy + .5);
			    if (IN_DOMAIN(i,j)) {
			    	n_input_points++;
				/*fprintf(stderr, "\n>>>%d %d   %d   %f    %d", i,j, n_data_this_cell[i][j], z_matrix[i][j], n_input_points);*/
			    	z_matrix[i][j] = (z_matrix[i][j]*n_data_this_cell[i][j] + z) / (n_data_this_cell[i][j]+1);
			    	n_data_this_cell[i][j]++;
			    }
			    else n_not_used++;
			}
		}
		PRINT_WARNING("s%d values not used", n_not_used);
		free(n_data_this_cell);
		for (i=0; i<Ny; i++) free(n_data_this_cell[i]);
		break;
	  default:
		PRINT_ERROR("Interpolation mode '%d' not available.\n", mode_interp);
		break;
	}

	for (i=0; i<Ny; i++)  for (j=0; j<Nx; j++)  z_matrix[i][j] += add_random * ((((float) rand()) / ((float) RAND_MAX)) -.5);

	free_matrix (input_points, nmax_input_points);
	PRINT_DEBUG("Interpoled %d points in mode %d.", n_input_points, mode_interp);
	return (1);
}






int readinterplin (
	FILE *file, 	/*input x-y file*/
	float *vector, 	/*ouput y-array*/
	int nx, 	/*number of alements in array*/
	float xpos0, 	/*x min*/
	float xposf)	/*x max*/
{
	int 	i ;
	float	xpos, dxpos, a, b, aa, bb;
	char	linea[MAXLENLINE], *lin ;

	/*LINEARLY INTERPOLATES DATA FROM AN 1D x-y FILE.
	UNDERSTANDS ANY LINE STARTING WHITH TWO FLOATS AND SKIPS THE
	OTHER LINES (so, comments in files are allowed)*/

	dxpos = (xposf-xpos0) / (nx-1) ;
	while ((sscanf(lin=fgets(linea, MAXLENLINE-1, file), "%f %f", &aa, &bb)) < 2) 	if (lin==NULL) {aa=bb=0; break;}
	while ((sscanf(lin=fgets(linea, MAXLENLINE-1, file), "%f %f", &a, &b)) < 2)  	if (lin==NULL) {a=aa; b=bb; break;}
	for (i=0; i<nx; i++)
	{
		xpos = i*dxpos+xpos0;
		while (xpos>a) {
			aa=a, bb=b;
			if (fgets(linea, MAXLENLINE-1, file)==NULL) {a=xposf+.001*dxpos, b=bb; }
			else 	sscanf(linea, "%f %f", &a, &b) ;
		}
		if (xpos>=aa) {
			if (a!=aa) 	vector[i] = (xpos-aa)*(b-bb)/(a-aa)+bb;
			else		vector[i] = bb;
		} 
		else	vector[i]=bb;
	}
	return (1) ;
}





float ReSort_Array (float *array, int *orden, int Nx)
{
	/* SORTS ARRAY NODES GIVEN A PREVIOUS SORT */
	register int 	numorden, i, j, aux1;
	for (numorden=0; numorden < Nx-1; numorden++) {
		if ( array[orden[numorden]] < array[orden[numorden+1]] ) {
			aux1 = orden[numorden+1];
			for (i=numorden; i>=0 ; i--) {
				if ( array[orden[i]] >= array[aux1] )  break;
				orden[i+1] = orden[i];
				orden[i+1] = orden[i];
			}
			orden[i+1] = aux1;
		}
	}
	if  (verbose_level>=3) {
		for (numorden=0; numorden < Nx-1; numorden++) {
			if ( array[orden[numorden]] < array[orden[numorden+1]] ) {
				PRINT_ERROR("UNSORTED!");
	}}}

	/*Returns the maximum value:*/
	return (array[orden[0]]) ;
}





int SolveAlmostDiagonalTriangularEquationSystem (
		double **A, 		/* Almost diagonal matrix of coeficients of the sistem A·x = b */
		double *b, 		/* Independent term */
		int num_ecs, 		/* Number of ecuations (rows of the matrix A)*/
		int NDsre, int NDire, 	/* Number of upper and lower diagonals */
		float *x)		/* Returning solution vector */
{
	/* RESUELVE UN SISTEMA DE ECUACIONES CASI DIAGONAL Y TRIANGULARIZADO:  A·x = b */

	register int i, j;

	for (i=num_ecs-1; i>=0; i--)
	{
		x[i] = b[i];
		for (j=NDire+1; (j<NDire+NDsre+1)&&(i+j-NDire<num_ecs); j++)
		{
			x[i] -= A[i][j] * ( x [i+j-NDire] );
		}
		x[i] /= A[i][NDire];
	}
	return(1);
}



double surface_topo(
	double **array, 	/*Array of z-values*/
	double ymax, 
	double dx, 		/*x increment in minutes*/
	double dy, 		/*y increment in minutes*/
	int nx, 		/*x dimension of the array*/
	int ny)			/*y dimension of the array*/
{
	/*
	Calculates the surface of an array in meters or degrees (according 
	to if switch_geograph_coor) asuming triangles between points
	*/

	int	i, j, nt=0;
	double 	dx_grad, dy_grad, 
		dxtr,		/*length of x side of triangle*/
		dytr,		/*length of y side of triangle*/
		dxdy2, 
		dxdyside, 
		surface,	/*Return value of surface.*/
		surfaceSW_NE=0,	/*Surface with one direction of truncation.*/
		surfaceSE_NW=0;	/*Surface with the other direction of truncation.*/

	if (switch_geograph_coor) {
		dx_grad = dx;
		dy_grad = dy;
		dy *= Rearth *pi/180 ;
	}
	dxdy2=dx*dx+dy*dy;

#define AREA_TRIANGLE(area, a, b, c)	{/*Heron's formula*/ double semiperimeter; semiperimeter = (a + b + c)/2; area += sqrt(semiperimeter*(semiperimeter-a)*(semiperimeter-b)*(semiperimeter-c));} 
	/*Southwest-Northeast triangle truncation:*/
	for (i=0; i<ny; i++) {
		if (switch_geograph_coor) {
			double y;
			y = ymax - i * dy_grad;
			dx = dx_grad * Rearth*pi/180 * cos(y*pi/180);
			dxdy2=dx*dx+dy*dy;
		}
		for (j=0; j<nx; j++) {
			/*South-West triangle*/
			if (i<ny-1 && j>0) {
				dxtr = sqrt(SQUARE(dx) + SQUARE(array[i][j]-array[i][j-1]));
				dytr = sqrt(SQUARE(dy) + SQUARE(array[i][j]-array[i+1][j]));
				dxdyside = sqrt(dxdy2 + SQUARE(array[i+1][j]-array[i][j-1]));
				AREA_TRIANGLE(surfaceSE_NW, dxtr, dytr, dxdyside);
				nt++;
/*fprintf(stderr, "\n>1 %d %d %.1f %.1f %.1f %.1f %.1f %.1f %d", i,j, surfaceSE_NW, array[i][j], dxtr, dytr, dxdy2, dy, nt);*/
			}
			/*North-East triangle*/
			if (i>0 && j<nx-1) {
				dxtr = sqrt(SQUARE(dx) + SQUARE(array[i][j]-array[i][j+1]));
				dytr = sqrt(SQUARE(dy) + SQUARE(array[i][j]-array[i-1][j]));
				dxdyside = sqrt(dxdy2 + SQUARE(array[i-1][j]-array[i][j+1]));
				AREA_TRIANGLE(surfaceSE_NW, dxtr, dytr, dxdyside);
				nt++;
/*fprintf(stderr, "\n>2 %d %d %.1f %.1f %.1f %.1f %d", i,j, surfaceSE_NW, array[i][j], dxtr, dytr, nt);*/
			}
		}
	}


	/*Southeast-Northwest triangle truncation:*/
	for (i=0; i<ny; i++) {
		if (switch_geograph_coor) {
			double y;
			y = ymax - i * dy_grad;
			dx = dx_grad * Rearth*pi/180 * cos(y*pi/180);
			dxdy2=dx*dx+dy*dy;
		}
		for (j=0; j<nx; j++) {
			/*South-East triangle*/
			if (i<ny-1 && j<nx-1) {
				dxtr = sqrt(SQUARE(dx) + SQUARE(array[i][j]-array[i][j+1]));
				dytr = sqrt(SQUARE(dy) + SQUARE(array[i][j]-array[i+1][j]));
				dxdyside = sqrt(dxdy2 + SQUARE(array[i+1][j]-array[i][j+1]));
				AREA_TRIANGLE(surfaceSW_NE, dxtr, dytr, dxdyside);
				nt++;
/*fprintf(stderr, "\n>3 %d %d %.1f %.1f %.1f %.1f %d", i,j, surfaceSW_NE, array[i][j], dxtr, dytr, nt);*/
			}
			/*North-West triangle*/
			if (i>0 && j>0) {
				dxtr = sqrt(SQUARE(dx) + SQUARE(array[i][j]-array[i][j-1]));
				dytr = sqrt(SQUARE(dy) + SQUARE(array[i][j]-array[i-1][j]));
				dxdyside = sqrt(dxdy2 + SQUARE(array[i-1][j]-array[i][j-1]));
				AREA_TRIANGLE(surfaceSW_NE, dxtr, dytr, dxdyside);
				nt++;
/*fprintf(stderr, "\n>4 %d %d %.1f %.1f %.1f %.1f %d", i,j, surfaceSW_NE, array[i][j], dxtr, dytr, nt);*/
			}
		}
	}

	surface = (surfaceSE_NW + surfaceSW_NE) / 2; 

	if ((surfaceSE_NW-surfaceSW_NE)/surface > 1e-4) 
		PRINT_ERROR("SE-NW=%.5e [m2]\tSW-NE=%.5e [m2]\n", 
				surfaceSE_NW, surfaceSW_NE);
	PRINT_ERROR("Number of triangles (x2, both directions): %.1f\nSurface: %.5e  SE_NW=%.5e SW_NE=%.5e m2", 1.0*nt*1.0, surface, surfaceSE_NW, surfaceSW_NE);

	return (surface);
}



int TriangularizeAlmostDiagonalEquationSystem(
		double **A,		/* Almost diagonal matrix of coeficients of the sistem A·x = b */
		double *b, 		/* Independent term */
		int num_rows, 		/* Number of ecuations (rows of the matrix A)*/
		int NDsre, int NDire) 	/* Number of upper and lower diagonals */
{
	/*Toma una matriz en banda y la triangulariza devolviendo 
	la matriz sin las NDire diagonales inferiores.
	El valor de retorno es 1 si ha habido error.	*/

	register	int 	i, j, k, l, cont = 0, contant=0, ndiagonales=NDire+NDsre+1  ;
	register	double 	aux, fac;

	for (i=0 ;  i< num_rows-1 ; i++)
	{
		if (!A[i][NDire])	/*si el de la diagonal principal es 0*/
		{			/*    cambiar esa l¡nea por otra*/
			for (j=i+1; j<i+1+NDire; j++) {
				if (j<num_rows) {
				if (A[j][i+NDire-j]) {
					for (k=i+NDire-j ; k < ndiagonales ; k++) {	/*  la cambia por esta  */
						aux=A[j][k] ;
						if ((k-i+j) < ndiagonales) {A[j][k]=A[i][k-i+j]; A[i][k-i+j]=aux ;} 
						else 		A[j][k]=0 ;
					}
					break ;
				}
				}
			}
			if (j==(i+1+NDire)) return 1;		/*  si no ha podido => sist. indeterminado   */
		}
		for (j=i+1; (j<i+1+NDire) && (j<num_rows) ; j++)  {	/*  resta a las siguientes filas esta por un factor tal que les haga cero el primer termino  */
			fac=A[j][i+NDire-j]/A[i][NDire] ;
			for (k=i+NDire-j;k<ndiagonales;k++)  {
				if ((k+j-i)<ndiagonales) aux=A[i][k+j-i]; else aux=0;
				A[j][k] -= aux*fac;
			}
			b[j] -= b[i]*fac;
		}
	}
	return 0;
}


void sistecslintriang(float **a, float *b, int nx, float *x)
{
	register int i, j;

	for (i=nx-1; i>=0; i--)
	{
		x[i] = b[i];
		for (j=i+1; j<nx; j++)
		{
			x[i] -= a[i][j]*x[j];
		}
		x[i] /= a[i][i];
	}
}



int triangularizar(float **a, float *b, int nx)
{
	register int i, j, k, l;
	float aux, fac;

	for (i=0; i<nx-1; i++)
	{
		for (j=i+1; j<nx; j++)
		{
			if (!(a[i][i]))
			{
				for (k=i+1; k<nx; k++) {
					if (a[k][i]) {
						for (l=i;l<nx;l++) { aux=a[i][l]; a[i][l]=a[k][l]; a[k][l]=aux; }
						aux=b[i]; b[i]=b[k]; b[k]=aux;
						break; }
					if (k==(nx-1)) return 1;		/* error: sistema indeterminado */
				}
			}
			if (!a[j][i]) goto etiq;
			fac=a[j][i]/a[i][i];
			for (k=i;k<nx;k++)
			{
				a[j][k] -= a[i][k]*fac;
			}
			b[j] -= b[i]*fac;
			etiq:;
		}
	}
	return 0;
}




char *replace_word(char *s, char *old, char *new)
{
	char *ret;
	int i, count = 0;
	size_t newlen = strlen(new);
	size_t oldlen = strlen(old);

	for (i = 0; s[i] != '\0'; i++) {
		if (strstr(&s[i], old) == &s[i]) {
			count++;
			i += oldlen - 1;
		}
	}

	ret = malloc(i + count * (newlen - oldlen));
	if (ret == NULL)
	exit(EXIT_FAILURE);

	i = 0;
	while (*s) {
		if (strstr(s, old) == s) {
			strcpy(&ret[i], new);
			i += newlen;
			s += oldlen;
		} 
		else
			ret[i++] = *s++;
	}
	ret[i] = '\0';

	return ret;
}
