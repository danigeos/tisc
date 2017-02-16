#------------------------------------------------------------------------
#You may need to modify these variables
#C compiler:
CC	= gcc #gcc #cc 
#Fortran compiler (only needed if THIN_SHEET is defined below):
FC	= gfortran #fort77 /opt/intel/fce/10.1.008/bin/ifort -O3 -nofor-main -static #f90 #g77 #f90 
LIBS	= -lm -lc 
#Program parts to include:
DEFS =   -DSURFACE_TRANSPORT  #-DTHIN_SHEET 
# Things you might add to DEFS:
# -DTHIN_SHEET 
# -DSURFACE_TRANSPORT 

#Options depending on the compiler/system:
OPTS_mac = -g -O3 -w -std=c99 #-static #-O3 -Wuninitialized #-Wno-unused-result #-Wparentheses 
OPTS_linux = $(OPTS_mac)
OPTS_AIX_RS6000 = -g -O3 #-Q -qsrcmsg -qmaxmem=4000 #-v
OPTS_SUN = -xO5 #-fast #-O3 -xO5 -g -fast -W
#CHOSE YOUR SYSTEM:
OPTS	= $(OPTS_mac)
#------------------------------------------------------------------------
