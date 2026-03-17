/*
	Program to find the pair of moment & shear force boundary conditions
	which best fit the observed bathymetry of a trench using tao software. 
	Evaluates the x-y resulting points from projectname.pfl tao output file
	in respect to the known x-y points from  projectname.FIT. This file can
	optionally have a third column  containing the weight for each point. It
	can use two different ways of inversion: Gridding (-ag option) and
	Searching (-as). 'Gridding' try points in a domain  which progressively
	makes smaller; 'Searching' tries next point varying parameters as a
	function of last errors. When convergence is reached, it invoques tao to
	make a PostScript with best fitting.

	taofit ignores the three force parameter values in the parameters  file
	*.PRM. if -m, -s or -p are not specified, taofit will assume  a zero for
	the boundary force passed to tao.

				Daniel Garcia-Castellanos, I-1996
*/

#include "tao.h"

float 	evaluate_xy_points (FILE *file_true, FILE *file_guess, char mode);
int match_parameter(char *str1, char *str2, int show, int replace, char *line);

int main(int argc, char **argv)
{
	int	iter, param, iarg, 
		num_since_min_err=0, max_since_min_err,
		n_gridding = 4, nmom, nshe;
	float	moment, d_mo, d_mo_ant, 
		shear, d_sh, d_sh_ant, 
		tecforce = 0,
		converg_fact_searching = 1, 	/*Default for searching method; if <<1					, >0 is faster but worse */
		converg_fact_gridding = .75, 	/*Default for gridding method; if <<1, >0 is faster but worse */
		domain_mo_min=0,	domain_sh_min=0,
		error=1e5, error_ant=1e5, d_err, 
		error_min=1e5, moment_err_min=1e19, shear_err_min=1e15;
	char	command[MAXLONLINE], mode='X', 
		projectname[MAXLONLINE],
		file_true_name[MAXLONLINE], file_guess_name[MAXLONLINE],
		algorithm = 'g';		/* s fo searching method; g for gridding */
	FILE	*file_true, *file_guess;


	param =	1;

	d_mo = 	-6e16;
	d_sh =  -4e12;
	moment=	-1e17;
	shear = -1e13;

	/*Interpreting command line*/
	for (iarg=1; iarg<argc; iarg++) { 
		if (argv[iarg][0] == '-') {
			int 	ilet;
			float 	value;
			char 	prm[MAXLONLINE];

			for (ilet=2; ilet < strlen(argv[iarg])+2; ilet++) prm[ilet-2] = argv[iarg][ilet];
			value=atof(prm);

			switch (argv[iarg][1]) {
				case 'a':
					algorithm = argv[iarg][2];
					break;
				case 'c':
					converg_fact_searching = value;
					converg_fact_gridding = value;
					break;
				case 'm':
					moment	= value;
					break;
				case 'p':
					tecforce = value;
					break;
				case 's':
					shear	= value;
					break;
				case 'M':
					d_mo = value;
					break;
				case 'S':
					d_sh = value;	/*Depth of cooling plate*/
					break;
				case 'W':
					mode = 'W';
					break;
			}
		}
		else {
			strcpy(projectname, argv[iarg]);
		}
	}
	if (argc<2) {
		fprintf(stderr, "Syntax:\t%s  <project_name> [-a<f|g>] [-m<moment>] [-p<tecforce>] [-s<shear>] [-M<d_mo>] [-S<d_sh>] [-W]\n", argv[0]);
		fprintf(stderr, "  Two columns (x-y) will be read from files.\n");
		fprintf(stderr, "  True-file can optionally be weighted in a third column.\n");
		exit(0);
	}

	fprintf(stdout, "\n%s: project name is '%s'. tecforce=%.3e", argv[0], projectname, tecforce);

	if (algorithm == 's') max_since_min_err = 16;
	if (algorithm == 'g') max_since_min_err = 1000;

	if (algorithm == 'g') {
		nmom = nshe = 0; 
		d_mo = fabs(d_mo);  d_sh = fabs(d_sh);
		domain_mo_min = moment - d_mo/2;
		domain_sh_min = shear - d_sh/2 ;
		moment = domain_mo_min + nmom / (n_gridding-1) * d_mo;
		shear = domain_sh_min + nshe / (n_gridding-1) * d_sh;
	}

	sprintf(file_true_name, "%s.FIT", projectname);
	sprintf(file_guess_name, "%s.pfl", projectname);
	fprintf(stdout, "\n\tIter.  \tmoment  \tshear   \tError\n");

	for (iter=0; ; iter++) {
		sprintf(command, "tao %s -S -m%.3e -s%.3e -p%.3e", argv[1], moment, shear, tecforce);
		system(command);

		if ((file_true = fopen(file_true_name, "rt")) == NULL) {
			fprintf(stderr, "\nInput file '%s' not found.\n", file_true_name);
			exit(0);
		}
		if ((file_guess = fopen(file_guess_name, "rt")) == NULL) {
			fprintf(stderr, "\nInput file '%s' not found.", file_guess_name);
			exit(0);
		}

		error = evaluate_xy_points(file_true, file_guess, mode);

		d_err = error-error_ant;

		fprintf(stdout, "\t%d\t%.3e\t%.3e\t", iter, moment, shear);
		/*fprintf(stdout, "  %.2e  %.2e  %.2f\t", d_mo, d_sh, d_err);*/
		fprintf(stdout, "%.2f\n", error);

		if (error<error_min) {
			error_min = error;
			num_since_min_err = 0;
			moment_err_min = moment;
			shear_err_min =  shear;
		}
		else 	num_since_min_err++;

		if (algorithm == 's' && num_since_min_err > max_since_min_err) {
			d_mo *= 1.3;
			d_sh *= 1.3;
			moment = moment_err_min;
			shear =  shear_err_min;
			num_since_min_err = 0;
			fprintf(stdout, " * ");
		}

		fclose(file_guess);
		fclose(file_true);

		/*Modify B.C.*/
		if (algorithm == 's' && iter>0){ switch(param) {
			case 1:
				d_mo_ant = d_mo;
				d_mo = -(d_err)/fabs(d_err)*d_mo;
				if (d_mo*d_mo_ant < 0) d_mo *= .9;
				d_mo *= converg_fact_searching;
				break;
			case 2:
				d_sh_ant = d_sh;
				d_sh = -(d_err)/fabs(d_err)*d_sh;
				if (d_sh*d_sh_ant < 0) d_sh *= .9;
				d_mo *= converg_fact_searching;
				break;
		}}

		if (algorithm == 'g') {
			nmom ++;
			if (nmom == n_gridding) { nmom=0; nshe++;}
			if (nshe == n_gridding) { 
				nmom = nshe = 0;
				d_mo *= converg_fact_gridding; 	d_sh *= converg_fact_gridding;
				domain_mo_min = moment_err_min - d_mo/2 ;
				domain_sh_min = shear_err_min - d_sh/2  ;
			}
			moment = domain_mo_min + (nmom+.0) / (n_gridding-1) * d_mo;
			shear = domain_sh_min + (nshe+.0) / (n_gridding-1) * d_sh;
		}

		error_ant = error;

		if (algorithm == 's') {
			switch(param) {
				case 1:
					shear += d_sh;
					param = 2;
					break;
				case 2:
					moment += d_mo;
					param = 1;
					break;
			}
		}


		if (error<20) {
			fprintf(stdout, "\nConverged!\n");
			break;
		}
		if (fabs(d_mo)<2e15 && fabs(d_sh)<2e10) {
			fprintf(stdout, "\nMin_Err\t%.3e\t%.3e\t%.3e\t%.2f\n", tecforce, moment_err_min, shear_err_min, error_min);
			fprintf(stderr, "\nLittle variations in parameters!\n");
			break;
		}
	}
	sprintf(command, "tao %s -S -m%.3e -s%.3e -p%.3e -P", argv[1], moment_err_min, shear_err_min, tecforce);
	system(command);
}
