/*
	Uses the function 'outin' to know whether a point is inside or outside
	of a given  polygon.
*/


#include "geomodel.h"		/*Most general definitions and types*/


int 	outin (float x, float y, float *polygonX, float *polygonY, int number_of_points_in_polygon);
float 	gravanompolyg (float *, float *, int, float, float, float);
float 	geoidanompolyg (float *, float *, int, double, double, double);
float 	dist_xy_pol(
	float x, 		/* x coordinate of the point to be inspected	*/
	float y, 		/* y     "    					*/
	float *polygonX, 	/* polygon x coordinates (polygonX[np])*/
	float *polygonY, 	/* polygon y coordinates (polygonY[np])*/
	int np);		/*Number of points in polygon.*/


main(int argc, char **argv)
{
	int	np=10000, i, result, iarg;
	float	area=0, 
		*polygonX, *polygonY, 
		a, b, x, y, 
		distance, gravanom, geoideanom;
	char 	inputfilename[MAXLONFICH];
	FILE 	*file=NULL;

	/*Takes memory for polygon*/
	polygonX = (float *) calloc(np, sizeof(float));
	polygonY = (float *) calloc(np, sizeof(float));

	verbose_level = 1;
	
	/*Interpreting command line*/
	for (iarg=1; iarg<argc; iarg++) {
		if (argv[iarg][0] == '-') {
			int 	ilet;
			float 	value;
			char 	prm[MAXLONLINE];

			for (ilet=2; ilet < strlen(argv[iarg])+2; ilet++) prm[ilet-2] = argv[iarg][ilet];
			value=atof(prm);

			switch (argv[iarg][1]) {
				case 'h':
					syntax();
				case 'V':
					verbose_level = value;
					break;
				default:
					fprintf(stderr, "\nWarning in %s: incomprehensible parameter '%s'.", argv[0], argv[iarg]);
					break;
			}
		}
		else {
			if ((file = fopen(argv[iarg], "rt")) == NULL) {
				fprintf(stderr, "\nERROR: Input polygon file '%s' not found.", argv[iarg]);
			}
		}
	}

	if (file == NULL) {
		fprintf(stderr, "\nERROR: Input x-y polygon file not found.");
		syntax();
		exit(0);
	}
	for (i=0;;i++) {
		TAKE_LINE_2(polygonX[i], polygonY[i])
	}

	fin_lectura:;
	np = i;

	for (i=0;i<np;i++) {
		float xi, yi, xn, yn;
		xi = polygonX[i]; yi = polygonY[i];
		if (i!=np-1) {	xn = polygonX[i+1]; yn = polygonY[i+1]; }
		else         {	xn = polygonX[0];   yn = polygonY[0];   }
		area += xi*yn-xn*yi;
	}
	area = .5*fabs(area);

	if (verbose_level>=1) fprintf(stderr, "\nPoints in polygon: %d", np);
	if (verbose_level>=2) fprintf(stderr, "\nPolygon area: %.3e", area);
	if (verbose_level>=2) fprintf(stderr, "\nEnter points to be inspected <x, y>: \n");

	for (i=0; ;) {
		char	linea[MAXLONLINE], *lin, zstr[MAXLONLINE] ;
		int nfields;
		
		lin=fgets(linea, MAXLONLINE-1, stdin);
		if ( lin == NULL ) { 
			/*Stops with ^D or end of file */
			if (verbose_level>=3) fprintf(stderr, "\nEnd of file."); 
			if (verbose_level>=2) fprintf(stderr, "\t%d points were checked.\n", i);
			exit(0) ;
		}
		nfields=sscanf(lin, "%f %f %s", &x, &y, zstr);
		if (nfields >=2) {
			result = outin(x, y, polygonX, polygonY, np);
			distance = dist_xy_pol(x, y, polygonX, polygonY, np);
			fprintf (stdout, "%d\t%s", result, linea);
			i++;
			if (verbose_level>=2) {
				if (result==1) 	fprintf(stderr, "%3d Inside;", i);
				else 		fprintf(stderr, "%3d Outside;", i);
				fprintf(stderr, " Dist.=%.2f\n", distance);
			}
		}	
	}
}


syntax () {
	fprintf(stderr, "\n\nSyntax: \t outin <polygonfile>  [-h] [-V<level>]");
	fprintf(stderr, "\n\nDetermines whether x,y points fall inside or outside a given polygon.");
	fprintf(stderr, "\nReads polygon from <polygonfile> in two columns (x,y); "
			"\nReads points (x,y,...) from stdin.");
	fprintf(stderr, "\nWrites [0|1], x, y, ... in stdout. 0 means outside; 1 inside.");
	fprintf(stderr, "\n-V provides additional information such as the polygon area.\n");
	AUTHORSHIP;
}
