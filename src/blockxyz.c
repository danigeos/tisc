/*
  Makes a xy grid from xyz sparse points.
*/
#include <tisc.h>

int main(int argc, char **argv)
{
	int	i, iarg, j, n_input_points=0, **n_values, n_rad_x=0, n_rad_y=0;
	FILE 	*file;
	char 	mode, filename[MAXLONFICH]="--", zstr[MAXLONLINE];
	float	xmin, xmax, ymin, ymax;
	float	x, y, z, default_value=0, **values, radius=0;

	if (argc<=2) {argc=2; argv[1]="-h";}
	for (iarg=1; iarg<argc; iarg++) {
		int 	ilet;
		float 	value, value2;
		char 	prm[MAXLONLINE], prm2[MAXLONLINE], *ptr;
		if (argv[iarg][0] == '-') {
			for (ilet=2; ilet < strlen(argv[iarg])+2; ilet++) 
				prm[ilet-2] = argv[iarg][ilet];
			for (ilet=3; ilet < strlen(argv[iarg])+2; ilet++) 
				prm2[ilet-3] = argv[iarg][ilet];
			value  = atof(prm);
			value2 = atof(prm2);
			switch (argv[iarg][1]) {
				case 'd':
					default_value = atof(prm);
					break;
				case 'h':
					syntax (argc, argv);
					exit(0);
					break;
				case 'M':
					mode = prm[0];
					radius = value2;
					if (verbose_level>=2) fprintf(stderr, "\nradio : %s ; %s ; %f", prm, prm2, value2);
					break;
				case 'N':
					Nx = atoi(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) Ny = atoi(ptr);
					else Ny = Nx;
					break;
				case 'R':
					xmin = atof(strtok(prm, "/"));
					xmax = atof(strtok(NULL, "/"));
					ymin = atof(strtok(NULL, "/"));
					ymax = atof(strtok(NULL, "/"));
					break;
				case 'V':
					verbose_level = 1;
					if (argv[iarg][2]) verbose_level = value;
					break;
			}
		}
		else {
			strcpy(filename, argv[iarg]);
		}
	}

	dx = (xmax-xmin)/(Nx-1);
	dy = (ymax-ymin)/(Ny-1);
	if (radius) {
		n_rad_x = floor(radius/dx+.5); 
		n_rad_y = floor(radius/dy+.5); 
		if (verbose_level>=2) fprintf(stderr, "\nRadius in x,y: %d %d", n_rad_x, n_rad_y);
	}

	n_values = calloc (Ny, sizeof(int *));
	for (i=0; i<Ny; i++) n_values[i] = calloc (Nx, sizeof(int));
	values = calloc (Ny, sizeof(float *));
	for (i=0; i<Ny; i++) values[i] = calloc (Nx, sizeof(float));
	for (i=0; i<Ny; i++) for (j=0;j<Nx;j++) values[i][j]= default_value;

	if (verbose_level>=2) fprintf(stderr, "\nFile: '%s'; Nx=%d; Ny=%d; R=%.1f/%.1f/%.1f/%.1f;  mode='%c'\n", filename, Nx, Ny, xmin, xmax, ymin, ymax, mode);

	if ((file = fopen(filename, "rt")) == NULL) {
		file = stdin;
	}


	for (;;) {
		char	linea[MAXLONLINE], *lin, zstr[MAXLONLINE];
		int 	nfields;
		lin = fgets(linea, MAXLONLINE-1, file);
		if ( lin == NULL ) { 
			/*Stops with ^D or end of file */
			if (verbose_level>=2) fprintf(stderr, "\nEnd of file."); 
			if (verbose_level>=2) fprintf(stderr, "\t%d points were checked.\n", n_input_points);
			break;
		}
		nfields=sscanf(lin, "%f %f %f", &x, &y, &z);
		if (nfields >= 3) {
			int ii, jj;
			j = floor((x-xmin)/dx+.5);
			i = floor(Ny-1-(y-ymin)/dy+.5);
			if (verbose_level>=3) fprintf(stderr, "\t\t(i,j)= %d %d\n", i, j);
			if (INSIDE_BORDER(i,j)) {
				switch (mode) {
				    case 'a':
					for (ii=MAX_2(0,i-n_rad_y); ii<=MIN_2(Ny-1,i+n_rad_y);ii++) 
					for (jj=MAX_2(0,j-n_rad_x); jj<=MIN_2(Nx-1,j+n_rad_x);jj++) {
						values[ii][jj] = (values[ii][jj]*n_values[ii][jj]+z)/(n_values[ii][jj]+1);
						n_values[ii][jj]++;
					}
					break;
				    case 'm':
					for (ii=MAX_2(0,i-n_rad_y); ii<=MIN_2(Ny-1,i+n_rad_y);ii++) 
					for (jj=MAX_2(0,j-n_rad_x); jj<=MIN_2(Nx-1,j+n_rad_x);jj++) {
						values[ii][jj] = MIN_2(z, values[ii][jj]);
						n_values[ii][jj]++;
					}
					break;
				    case 'M':	
					for (ii=MAX_2(0,i-n_rad_y); ii<=MIN_2(Ny-1,i+n_rad_y);ii++) 
					for (jj=MAX_2(0,j-n_rad_x); jj<=MIN_2(Nx-1,j+n_rad_x);jj++) {
						values[ii][jj] = MAX_2(z, values[ii][jj]);
						n_values[ii][jj]++;
					}
					break;
				}
				if (n_values[i][j]==0) values[i][j] = z;
				n_input_points++;
			}
		}	
	}
	for (i=0;i<Ny;i++) for (j=0;j<Nx;j++)  { 
		x = xmin+j*dx; y = ymin+(Ny-i-1)*dy;
		fprintf (stdout, "%.2f\t%.2f\t%.2f\t%d\n", x, y, values[i][j], n_values[i][j]);
	}

	fclose(file);

	if (verbose_level>=3) fprintf(stderr, "\nUsed %d xyz points from input.\n", n_input_points);
}



syntax(int argc, char **argv)
{
	fprintf(stderr, 
		"\nSyntax:  %s [pointsfile] [-d<empty_val>] [-h] -M<mode>[<radius>] -N<Nx>[/<Ny>]\n\t-R<xmin/xmax/ymin/ymax> [-V<level>]", argv[0]);
	fprintf(stderr, 
		"\nReads x,y,z from pointsfile or stdin.\nWrites one x,y,z_value,n for each cell of the given grid (according to <mode>). "
		"\n<mode> is 'm' for minimum, 'M' for maximum, 'a' for average. "
		"\n<Nx> and <Ny> are used in this way:  dx = (xmax-xmin)/(Nx-1);  written x = xmin+j*dx "
		"\n<radius> is the distance along x & y axis to compute those (in input x,y units; but default is one cell)."
		"\n<empty_val> is the initial value set for the interpolated cells (default is 0).\n"
	);
}
