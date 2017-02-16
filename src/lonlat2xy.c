/*
	Converts x, y, z from geographic degrees to meters.
	Assumes a flat Earth.
*/


#include "geomodel.h"

int outin(
	float x, 		/* x coordinate of the point to be inspected	*/
	float y, 		/* y     "    					*/
	float *poligonX, 	/* Poligon x coordinates (poligonX[np])		*/
	float *poligonY, 	/* Poligon y coordinates (poligonY[np])		*/
	int   np);		/* Number of points in poligon.			*/



main(int argc, char **argv)
{
	int	i, j, iarg, np=10000;
	float	mean_lat=37,	
		x0=-7, y0=36, xx, yy, zz, 
		escala_mediax, 		/* m/º */ 
		escala_mediay, 
		*poligonX, *poligonY,
		poligon_value=0, 
		xm, ym, zm, alt_shift = 100;
	char	zstr[MAXLONLINE];
	BOOL	switch_reverse=NO, switch_outside_inside=NO, 
		switch_change_z_val=NO, switch_km_io=NO, 
		switch_divide_at_360=NO;
	FILE 	*filein=NULL, *filelimit=NULL;

	if (argc<2) syntax(argc, argv); 

	poligonX = (float *) calloc(np, sizeof(float));
	poligonY = (float *) calloc(np, sizeof(float));

	for (iarg=1; iarg<argc; iarg++) {
		int 	ilet;
		float 	value, value2;
		char 	prm[MAXLONLINE], prm2[MAXLONLINE];
		
		if (argv[iarg][0] == '-') {
			for (ilet=2; ilet < strlen(argv[iarg])+2; ilet++) 
				prm[ilet-2] = argv[iarg][ilet];
			for (ilet=3; ilet < strlen(argv[iarg])+2; ilet++) 
				prm2[ilet-2] = argv[iarg][ilet];
			value  = atof(prm);
			value2 = atof(prm2);
			switch (argv[iarg][1]) {
				case 'd':
					switch_divide_at_360 = YES;
					break;
				case 'I':
					if ((filelimit = fopen(prm, "rt")) == NULL) {
						fprintf(stderr, "\nInput poligon file '%s' not found.", prm);
						fprintf(stderr, "\nTwo columns (x-y) will be read from file.\n");
						exit(0);
					}
					switch_change_z_val = YES;
					for (i=0;;i++) {
						float  a, b;
						int    nread;
						char   linea[MAXLONLINE], *lin ;
						do {
							lin=fgets(linea, MAXLONLINE-1, filelimit);
							if (lin==NULL) break;
							nread=sscanf(lin, "%f %f", &a, &b); 
						} while (nread<2);
						if (lin==NULL) break;
						poligonX[i]=a;
						poligonY[i]=b;
					}
					np = i;
					break;
				case 'i':
					poligon_value = value;
					break;
				case 'k':
					switch_km_io = YES;
					break;
				case 'L':
					mean_lat = value;
					break;
				case 'r':
					switch_reverse = YES;
					break;
				case 's':
					switch_outside_inside = YES;
					break;
				case 'V':
					verbose_level = 1;
					if (argv[iarg][2]) verbose_level = value;
					break;
				case 'x':
					x0 = value;
					break;
				case 'y':
					y0 = value;
					break;
				default: 
					fprintf(stderr, "\nOption '%s' not found.", argv[iarg]);
					syntax(argc, argv);
					exit(0);
			}
		}
		else {
			if ((filein = fopen(argv[iarg], "rt")) == NULL) {
				fprintf(stderr, "\nInput file '%s' not found.", argv[iarg]);
				fprintf(stderr, "\nAt least two columns (x-y) will be read from file.\n");
				exit(0);
			}
		}
	}

	escala_mediay = 4e7 / 360;
	escala_mediax = escala_mediay * cos(mean_lat*pi/180);
	if (filein == NULL) {
		filein = stdin;
		fprintf(stderr, "\nReading point from stdin...");
	}

	for (j=0; ; j++) {
		char	linea[MAXLONLINE], *lin ;
		int 	nfields;
		sprintf(zstr, "");
		lin=(fgets(linea, MAXLONLINE-1, filein));
		if ( lin == NULL ) { 
			/*Stops with ^D or end of file */
			if (verbose_level>=2) fprintf(stderr, 
						"\nEnd of file."
						"\nPoints converted: %d\n", j);
			if (filein != stdin) fclose(filein);
			exit(0) ;
		}
		nfields=sscanf(lin, "%f %f %s", &xx, &yy, zstr) ;
		if (switch_km_io && switch_reverse) {
			xx *= 1e3; yy *= 1e3;
		}
		if (nfields >=2) {
			if (!switch_reverse) {
				if (!switch_divide_at_360) {
					if (xx>180 && xx<360) xx = xx-360;
				}
				else {
					if (xx>-180 && xx<0)  xx = xx+360;
				}
				xm = (xx - x0) * escala_mediax ;
				ym = (yy - y0) * escala_mediay ;
				zm = zz ;
			}
			else {
				xm = xx / escala_mediax + x0 ;
				ym = yy / escala_mediay + y0 ;
			}
			if (switch_change_z_val && outin(xx, yy, poligonX, poligonY, np)!=switch_outside_inside) 
				sprintf(zstr, "%f", poligon_value);
			if (switch_km_io && !switch_reverse) 
				fprintf(stdout, "%f\t%f\t%s\n", xm/1e3, ym/1e3, zstr);
			else 
				fprintf(stdout, "%f\t%f\t%s\n", xm, ym, zstr);
		}
		else {
			fprintf(stdout, "%s", linea);
		}
	}
}


syntax(int argc, char **argv)
{
	fprintf(stderr, "\nSyntax:  %s\t<xy_file> [-I<xy_pol_file>] [-i<value>]", argv[0]);
	fprintf(stderr, "\n                 \t[-L<mean_lat>] [-r] [-s] [-V[<level>]] [-x<lon0>] [-y<lat0>]\n");
	fprintf(stderr, "\nReads lon-lat (in decimal degrees) from stdin or [xy_file] and writes to stdout their conversion to meters."
			"\n-r reverts the conversion (xy to lonlat).");
	fprintf(stderr, "\nThe values lon0 and lat0 define the lon-lat of the point x,y=0,0.");
	fprintf(stderr, "\nmean_lat is the central latitude or your region of conversion.");
	fprintf(stderr, "\nWrites in stdout: (x,y,z)  or (x,y,value) when inside the poligon defined ");
	fprintf(stderr, "\nin -I (outside if -s is switched).");
	fprintf(stderr, "\n-d means that longitudes of -180 and 180 (instead of 0 and 360, which is the "
			"\ndefault) will have the same converted x value.");
	fprintf(stderr, "\n-k for input/output in km instead of m.\n");
	fprintf(stderr, "\n");
}
