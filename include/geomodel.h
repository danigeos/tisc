/*
		GENERAL INCLUDE FILE FOR GEOPHYSICAL MODELS
*/

#include "universal.h"			/*Most general definitions and types*/


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

