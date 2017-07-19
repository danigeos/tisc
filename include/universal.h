/*
	include file for most of my programs.
*/

#include "types_n_defs.h"

#define	LENGTHVERS	40

float 	
	g	= 9.81, 			/*Earth's surface average gravity acceleration [m·s-2]*/
	pi	= 3.14159265, 		/*PI number*/
	sqrt2	= 1.4142136, 	/*Square root of 2*/
	CGU	= 6.6732E-11, 		/*Constant of Universal Gravitation [m3·s-2·kg-1]*/
	number_e = 2.7182818, 	/*Number e*/
	densice = 920, 			/*Ice density*/
	denswater = 1020, 		/*Water density*/
	E	= 7E+10, 			/*Young modulus [N/m2]; Gerbault, 2000. At what stress level is the central Indian Ocean lithosphere buckling? EPSL, 178; Panteleyev, A. N.  & M. Diament, 2007. GJI 114-220. Influence of Some Rheological Parameters On Flexure of the Oceanic Lithosphere; Other authors use 1e11 (e.g., Minshull)*/
	nu	= .25, 				/*Poisson coefficient*/
	Rearth	= 6378e3 ;		/*Earth's radius*/


int	
	Nx, Ny, Nz,  			/*Number of grid knots in x & z directions)*/
	verbose_level;

float
	dx, dy, dz, dxy;		/*Horizontal & vertical grid intervals [m]*/


char	
	version[LENGTHVERS],			/*Version of the program*/
	version_input[LENGTHVERS];		/*Version actually found in input files*/


BOOL	
	switch_geograph_coor,			/*1 if x-y are geographycal coordinates in decimal degrees*/
	switch_ps, 						/*SI if postscript file will have to be wroten*/
	switch_write_file;				/*SI if results file must be wroten*/

