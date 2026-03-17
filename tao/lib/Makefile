################################################################################
########                       LIBRARIES MAKE FILE                      ########
################################################################################

#------------------------------------------------------------------------
#You should not modify these directories
LIB	= .
OBJ	= .
INCL	= ../include
CC	= gcc #gcc cc
FC	= gfortran #/opt/intel/fce/10.1.008/bin/ifort #f90 -O2 -YEXT_NAMES=LCS -s -B108 -YCFRL=1
LIBS	= -lm -lc 
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#You may need to modify these variables
#Options depending on the compiler:
OPTS_linux = -g -O3 -std=c99 #-g
OPTS_AIX_RS6000 = -g -O3 -Q #-qsrcmsg -v
OPTS_SUN = -xO5 #-fast #-g -xO5 -fast -O -O3
#Chose your compiler:
OPTS	= $(OPTS_linux)
#------------------------------------------------------------------------



libreria.o : libreria.c
	$(CC) -c -o libreria.o $(OPTS) -I$(INCL) libreria.c

thin_sheet : thin_sheet.f sistbanda.f
	$(FC) -c thin_sheet.f sistbanda.f

