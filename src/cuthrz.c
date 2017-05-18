/*
	Cuts a *.hrz TISC output file producing a *.pfl tao output file
*/

#include "tisc.h"
#include "tisclib.c"
#include "tiscio.c"

int 	read_file_output_Blocks ();


int main(int argc, char **argv)
{
	int i, iarg;
	FILE *commstdout;

	commstdout = stderr;
	zini = 0;
	switch_write_file = YES;

	syntax(argc, argv);
	for (iarg=1; iarg<argc; iarg++) {
		char prm[100], prm2[100];
		int ilet;
		float value, value2;
		if (argv[iarg][0] == '-') {
			for (ilet=2; ilet < strlen(argv[iarg])+2; ilet++) 
				prm[ilet-2] = argv[iarg][ilet];
			for (ilet=3; ilet < strlen(argv[iarg])+2; ilet++) 
				prm2[ilet-3] = argv[iarg][ilet];
			value  = atof(prm);
			value2 = atof(prm2);
			switch (argv[iarg][1]) {
				case 'V':
					verbose_level=value;
					break;
			}
		}
		else strcpy(projectname, argv[iarg]);
	}

	read_file_output_Blocks();

	w = alloc_matrix(Ny, Nx);

	write_file_cross_section();

	return(1);
}



int syntax (int argc, char **argv)
{
	if (argc<2) {
		fprintf(stderr, 
			"\nSyntax:   %s  prj  [-V<level>]"
			"\n\tReads 'prj.hrz' (a TISC horizon file) and cuts it along a cross-section path that is read from 'prj.PRFL'."
			"The cross section is written in 'prj.pfl'. Any pre-existing 'prj.pfl' file will be overwritten! "
			"\n\t'prj.PRFL' should contain two x,y columns in m.\n", argv[0]);
		AUTHORSHIP
		exit (1);
	}
}
