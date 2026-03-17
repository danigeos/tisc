#ifndef GEOMODEL_H
#define GEOMODEL_H

/*
	GENERAL INCLUDE FILE FOR GEOPHYSICAL MODELS
	(Consolidated from types_n_defs.h and universal.h)
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <stdbool.h>
#include <time.h>

#define	MAXLENFILE	256			/*Max. length for filenames*/
#define	MAXLENLINE	1024			/*Max. length for character strings, input lines, ...*/
#define	MAXLONFICH	MAXLENFILE
#define	MAXLONLINE	MAXLENLINE
#define secsperMa 	(365.24e6*24*3600)	/*Converts Myr to seconds*/
#define Matosec 	secsperMa
#define secsperyr 	(365.24  *24*3600)	/*Converts years to seconds*/
#define NO_DATA		-9999
#define SIGNAL		-9999

#define PRINT_ERROR(...)	{if (verbose_level>=0) {fprintf(stderr, "\nERROR: In %s: ", __func__); fprintf(stderr, __VA_ARGS__);}} 	/**/
#define PRINT_SUMLINE(...)	{if (verbose_level>=1) {fprintf(stdout, "\n  "__VA_ARGS__);} fflush(stdout);} 		/**/
#define PRINT_INFO(...) 	{if (verbose_level>=2) {fprintf(stdout, "\nInfo: "__VA_ARGS__);} fflush(stdout);} 	/**/
#define PRINT_WARNING(...)	{if (verbose_level>=3) {fprintf(stderr, "\nWarning: In %s: ", __func__); fprintf(stderr, __VA_ARGS__);}}	/**/
#define PRINT_DEBUG(...)	{if (verbose_level>=4) {fprintf(stderr, "\nDebug: In %s: ", __func__); fprintf(stderr, __VA_ARGS__);}} 	/*Used to track position in the code*/
#define PRINT_DEBUGPLUS(...){if (verbose_level>=5) {fprintf(stderr, "\nDebug+: In %s: ", __func__); fprintf(stderr, __VA_ARGS__);}} 	/*Used to help on code debugging*/
#define PRINT_ARRAY_INFO(array, name, units) {\
	    float max=-1e19, min=1e19, vol=0;\
	    for (int i=0; i<Nx; i++) {\
	    	    vol += array[i]*dx;\
	    	    if (max<array[i])  max=array[i];\
	    	    if (min>array[i])  min=array[i];}\
	    PRINT_SUMLINE("%s:  max= %.1f %s\tmin= %.1f %s\tintegr=%.2e %s*m", name, max, units, min, units, vol, units);}
#define PRINT_GRID_INFO(grid, name, units) {\
	    float max=-1e19, min=1e19, vol=0;\
	    for (int i=0; i<Ny; i++)    for (int j=0; j<Nx; j++) {\
	    	    vol += grid[i][j]*dx*dy;\
	    	    if (max<grid[i][j])  max=grid[i][j];\
	    	    if (min>grid[i][j])  min=grid[i][j];}\
	    PRINT_SUMLINE("%s:  max= %.1f %s\tmin= %.1f %s\tintegr=%.2e %s*m2", name, max, units, min, units, vol, units);}

/*Function definitions:*/
#define MAX_2(x,y)	(((x)>(y))? (x) : (y))	/*Gives maximum of two values*/
#define MIN_2(x,y)	(((x)<(y))? (x) : (y))	/*Gives minimum of two values*/
#define LIMIT(x,y,z)	(((x)<(y))? (y) : ((x)>(z))? (z) : (x))	/*Limits x in an interval [y,z]*/
#define SIGN(x)		(((x)<0)? (-1) : (((x)>0)? (+1) : (0)))
#define SQUARE(x)	((x)*(x))
/*Linearly interpolates y values given at two x coordinates. Constant value assumed out of the interval*/
#define LININTERP(x, x1,x2, y1,y2)	((x<=x1)? y1 : ((x<x2)? (y1+(y2-y1)/(x2-x1)*(x-x1)) : y2))

/*DEFINITIONS FOR 2D GRID DOMAINS*/
#define OUT_DOMAIN(i,j) 	((i<0)  || (i>=Ny)   || (j<0)  || (j>=Nx)  )
#define OUT_BORDER(i,j) 	OUT_DOMAIN(i,j)
#define IN_DOMAIN(i,j)		((i>=0) && (i<Ny)    && (j>=0) && (j<Nx)   )
#define NOT_OUT_BORDER(i,j)  	IN_DOMAIN(i,j)
#define INSIDE_BORDER(i,j)  	IN_DOMAIN(i,j)
#define AT_BORDER(i,j)		((i==0) || (i==Ny-1) || (j==0) || (j==Nx-1))
#define NOT_AT_BORDER(i,j)	((i>0)  && (i<Ny-1)  && (j>0)  && (j<Nx-1) )
#define NEIGHBOURS(i,j,ni,nj)	(fabs(i-ni) <= 1 && fabs(j-nj) <= 1 && (i != ni || j != nj))
#define DOMAIN_LIMIT(i,j)	{if (i<0) i=0;  if (i>=Ny) i=Ny-1;  if (j<0) j=0;  if (j>=Nx) j=Nx-1;};
#define BORDER_LIMIT(i,j)	DOMAIN_LIMIT(i,j)
#define BORDER_INDEX(i,j)	(((j)==0)? 3 : ((j)==Nx-1)? 2 : ((i)==0)? 0 : ((i)==Ny-1)? 1 : -1)
/*DEFINITIONS FOR 1D ARRAY DOMAINS*/
#define OUT_DOMAIN_1D(i)	((i<0)  || (i>=Nx))
#define IN_DOMAIN_1D(i)		((i>=0) && (i<Nx))
#define AT_BORDER_1D(i)		((i==0) || (i==Nx-1))
#define NOT_AT_BORDER_1D(i)	((i>0)  && (i<Nx-1))
#define DOMAIN_LIMIT_1D(i)	{if (i<0) i=0;  if (i>Nx-1) i=Nx-1;};

#define TEMP_FREEZE_WATER	273.15	/*K*/
#define GAMMA_GROUND		.003 	/*!! .0065 "lapse rate constant" in K/m*/
#define GAMMA_AIR		.008 	/*"Environmental lapse rate constant" in K/m*/
#define GAMMA			.0065 	/*K/m*/
#define TEMPERATURE_GROUND(topo)	(temp_sea_level - GAMMA_GROUND*(topo)) /*[K] air temperature at ground level*/
#define TEMPERATURE_AIR(topo,height)	(TEMPERATURE_GROUND(topo) - GAMMA_AIR*(height)) /*[K] air temperature*/
#define TEMPERATURE(altitude)	(temp_sea_level - GAMMA*(altitude)) /*[K]; "lapse rate constant" in K/m*/
#define TEMPERATURE_ICE(altitude)	(temp_sea_level - GAMMA*(altitude)) /*in [K]; "lapse rate constant" in K/m*/
#define IF_LAKE_IS_SEA(il) 	if (il) if (Lake[il].n_sd) if (topo[Lake[il].row_sd[0]][Lake[il].col_sd[0]]<sea_level && AT_BORDER(Lake[il].row_sd[0], Lake[il].col_sd[0]))

#define AUTHORSHIP		{ fprintf(stderr, "\n\t\t\t\t2018, Daniel Garcia-Castellanos\n");}

#define	LENGTHVERS	40

extern float g;
extern float pi;
extern float sqrt2;
extern float CGU;
extern float number_e;
extern float densice;
extern float denswater;
extern float E;
extern float nu;
extern float Rearth;
extern float viscwater;

extern int Nx, Ny, Nz, switch_ps, verbose_level;

extern float dx, dy, dz, dxy;

extern char version[LENGTHVERS];
extern char version_input[LENGTHVERS];

extern bool switch_geograph_coor;
extern bool switch_write_file;

#define	ET2RIG(x)	(E*pow(x,3.)/(12*(1.-nu*nu)))		/*Converts Elastic Thickness into Rigidity*/
#define	RIG2ET(x)	(pow((x)/E*12.*(1.-nu*nu), 1./3.))	/*Converts Rigidity into Elastic Thickness*/

extern int grav_anom_type, isost_model, water_load;

extern float Te_default, crust_thick_default, upper_crust_thick_default;
extern float densasthen, densmantle, denscrust, densinfill, denssedim, densenv;
extern float sea_level;
extern float temp_sea_level;
extern float Time, Timeini, Timefinal;
extern float dt, dt_eros, tau;

extern char projectname[MAXLONFICH];
extern char title[MAXLONLINE];

#endif /* GEOMODEL_H */
