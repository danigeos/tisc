/*
	include file for most of my programs.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
//#include <malloc.h>
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
/*TAKE_LINE_N reads N numbers from the first line with at least N numerical values. Put it inside a 'for' loop, it will break at the end of the file*/
#define TAKE_LINE_4(x, y, z, t)	{ char auxstr[MAXLENLINE], *lin; int nfields=0; while (nfields<4) {lin=fgets(auxstr, MAXLENLINE-1, file); if (lin==NULL) break; nfields=sscanf(lin, "%f %f %f %f", &x, &y, &z, &t);};	if (lin==NULL) break;}
#define TAKE_LINE_3(x, y, z)	{ char auxstr[MAXLENLINE], *lin; int nfields=0; while (nfields<3) {lin=fgets(auxstr, MAXLENLINE-1, file); if (lin==NULL) break; nfields=sscanf(lin, "%f %f %f",    &x, &y, &z);};   	if (lin==NULL) break;}
#define TAKE_LINE_2(x, y)	{ char auxstr[MAXLENLINE], *lin; int nfields=0; while (nfields<2) {lin=fgets(auxstr, MAXLENLINE-1, file); if (lin==NULL) break; nfields=sscanf(lin, "%f %f",       &x, &y);};      	if (lin==NULL) break;}
#define TAKE_LINE_1(x)		{ char auxstr[MAXLENLINE], *lin; int nfields=0; while (nfields<1) {lin=fgets(auxstr, MAXLENLINE-1, file); if (lin==NULL) break; nfields=sscanf(lin, "%f",          &x);};          	if (lin==NULL) break;}
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
#define GAMMA_AIR		.008 	/*"lapse rate constant" in K/m*/
#define GAMMA			.0065 	/*K/m*/
#define TEMPERATURE_GROUND(topo)	(temp_sea_level - GAMMA_GROUND*(topo)) /*[K] air temperature at ground level*/
#define TEMPERATURE_AIR(topo,height)	(TEMPERATURE_GROUND(topo) - GAMMA_AIR*(height)) /*[K] air temperature*/
#define TEMPERATURE(altitude)	(temp_sea_level - GAMMA*(altitude)) /*[K]; "lapse rate constant" in K/m*/
#define TEMPERATURE_ICE(altitude)	(temp_sea_level + 5*cos(Time/(100e3*secsperyr)*2*3.1415927) - GAMMA*(altitude)) /*in [K]; "lapse rate constant" in K/m*/
#define IF_LAKE_IS_SEA(il) 	if (il) if (Lake[il].n_sd) if (topo[Lake[il].row_sd[0]][Lake[il].col_sd[0]]<sea_level && AT_BORDER(Lake[il].row_sd[0], Lake[il].col_sd[0]))




#define AUTHORSHIP		{ fprintf(stderr, "\n\t\t\t\t2008, Daniel Garcia-Castellanos\n");}

/*YES=1;  NO=0. Defines boolean type for switches with values true or false (YES or NO)*/
typedef enum BOOLEAN	{NO, YES} BOOL;


struct BLOCK { 		/*for TISC*/
	float	**thick;		/*Present thickness at each x point*/
	float	**detr_ratio;		/*Only used for sediment Blocks: % of detrital sediment (non carbonatic)*/
	float	**detr_grsize;		/*Only used for sediment Blocks: grain size of the detrital sediment*/
	float	age;			/*Age of initial file reading*/
	float	density;		/*Density*/
	float	erodibility;		/*erosion parameter*/
	float	last_shift_x;		/*Previous x shift of Block*/
	float	last_shift_y;		/*Previous y shift of Block*/
	float	last_vel_time; 		/*Last time in which velocity changed*/
	float	shift_x;		/*Total x shift of Block*/
	float	shift_y;		/*Total y shift of Block*/
	float	time_stop;		/*Time to stop*/
	char 	type;			/*'T' means thin_sheet*/
	float	**vel_x;		/*Velocity in x direction*/
	float	**vel_y;		/*Velocity in y direction*/
	float	**visc;			/*Viscosity (only for thin sheet calculations)*/
	float	**viscTer;		/*Viscosity thermal term (only used for thin sheet calculations, in the first step)*/
};

struct BLOCK_1D {	/*for tAo*/
	float	*thick;			/*Present thickness at each x point*/
	float	*detr_ratio;		/*Only used for sediment Blocks: % of detrital sediment (non carbonatic)*/
	float	*detr_grsize;		/*Only used for sediment Blocks: grain size of the detrital sediment*/
	float	age;			/*Age of initial file read*/
	float	density;		/*Density*/
	float	erodibility;		/*erosion parameter*/
	float	last_shift;		/*Espected shift (not affected by finite differences discretization)*/
	float	last_vel_time; 		/*Last time in which velocity changed*/
	float	shift;			/*Total horizontal shift of Block*/
	float	time_stop;		/*Time in wich Block will stop*/
	char 	type;			/*'T' means thin_sheet*/
	float	vel;			/*Velocity at wich Block moves*/
};

struct GRIDNODE {
	int row;
	int col;
};

struct DRAINAGE {
	int dr_row;		/*row of the node to where drains*/
	int dr_col;		/*column of the node to where drains*/
	float discharge;	/*water flow through the node [m3/s]*/
	float masstr;		/*sediment load: mass exiting the cell [kg/s]*/
	float grainsize;	/*average grain size of the sedload[m]*/
	char type;		/*type (lake, river, sea, etc)*/
	int lake;		/*number of the lake: > 0 means is well defined; < 0 means is not still defined; 0 means it is not a lake*/
};

struct LAKE_INFO {		/*For lakes*/
	int n;			/*number of nodes INCLUDING SADDLES*/
	int *row;
	int *col;
	int n_sd;		/*number of saddles and transferring nodes*/
	int *row_sd;
	int *col_sd;
	float alt;		/*Altitude of the lake water level*/
	float vol;		/*Volume of the lake water body*/
};

struct CS2D {
	float *horiz;
	float x;
	float y;
	float l;
};

struct DRAINAGE_1D {
	int dr;		/*row of the node to where drains*/
	float discharge;	/*water flow through the node [m3/s]*/
	float masstr;		/*sediment load: mass exiting the cell [kg/s]*/
	float grainsize;	/*average grain size of the sedload[m]*/
	char type;		/*type (lake, river, sea, etc)*/
	int lake;		/*number of the lake: > 0 means is well defined; < 0 means is not still defined; 0 means it is not a lake*/
};

struct LAKE_INFO_1D {		/*For lakes*/
	int n;			/*number of nodes including saddles*/
	int *cell;
	int n_sd;		/*number of saddles and transferring nodes*/
	int *sd;
	float alt;		/*Altitude of the lake water level*/
	float vol;		/*Volume of the lake water body*/
};

