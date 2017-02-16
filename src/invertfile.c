/*
ESTE COMANDO YA EXISTE EN LINUX - SE LLAMA tac !!   (2007-11)
En Mac OSX no existe tac, pero se puede usar (2012-05):
| awk '{print NR,$0}' | sort -nr | sed 's/^[0-9]* //' |
O TAMBIEN:
tail -r

	INVIERTE EL ORDEN DE LAS FILAS DE UN FICHERO DE TEXTO Y LO ENVIA 
	A LA SALIDA STANDAR

		Daniel Garcia Castellanos	Feb. 1995
*/


#include "universal.h"

#define	MAXNUMLINEAS	50000
#define	MAX_LENGTH_LINE	2000

int main (int argc, char **argv)
{
	int 	cont, numlines, switch_info=0;
	FILE 	*filein;
	char 	fichinput[MAXLONFICH], linea[MAX_LENGTH_LINE], *vector_de_lineas[MAXNUMLINEAS] ;

        if (argc > 1 && strcmp(argv[1], "-h") == 0) switch_info=1;
	if (argc > 1) strcpy(fichinput, argv[1]);
	if ((filein = fopen(fichinput, "r")) == NULL) {
		if (switch_info) {
			fprintf(stderr, "\nInput file %s not found.", fichinput) ;
			fprintf(stderr, "\nSyntax: \t invertfile [<inputfile>] [-h]\n If <inputfile> is not given then stdin is read. Standard output.\n") ;
		}
		filein = stdin; 
	}
	for (cont=0; ;cont++) {
	 	if (fgets(linea, MAX_LENGTH_LINE, filein) == NULL) break;
		if (cont>MAXNUMLINEAS) {
			fprintf (stderr, "\nERROR in %s: Too many lines. \nChange MAXNUMLINEAS value at source code file and re-compile.\n", 
			argv[0]);
			exit(0);
		}
		vector_de_lineas[cont] = (char *) malloc ((sizeof(char)+1) * strlen(linea));
		strcpy(vector_de_lineas[cont] , linea);
	}
	numlines = cont ;
	close(filein);

	if (switch_info) fprintf(stderr, "\nInverts the order of the %i lines in '%s'.\n", numlines, fichinput);
	for (cont=numlines-1; cont>=0 ; cont--) {
	 	if (fprintf(stdout, "%s", vector_de_lineas[cont]) < 0) { fprintf(stderr, "\nERROR trying to print at output file."); exit(0); }
	}
	exit(numlines);
}
