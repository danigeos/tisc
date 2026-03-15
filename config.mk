#------------------------------------------------------------------------
VERSION = TISC_2026-03-14_factorized  #Automatically updated in template.PRM !!
#You may need to modify these variables
#C compiler:
CC	= gcc #gcc #cc 
#Fortran compiler (only needed if THIN_SHEET is defined below):
FC	= gfortran #fort77 /opt/intel/fce/10.1.008/bin/ifort -O3 -nofor-main -static #f90 #g77 #f90 
#Program parts to include:
DEFS =   -DSURFACE_TRANSPORT  #-DTHIN_SHEET 
# Things you might add to DEFS:
# -DTHIN_SHEET 
# -DSURFACE_TRANSPORT 

# --- OpenMP Parallelization ---
# To turn OFF parallelization later, simply comment out the 4 OMP_ lines below:
OMP_PREFIX := $(shell /opt/homebrew/bin/brew --prefix libomp 2>/dev/null || brew --prefix libomp 2>/dev/null || echo /usr/local/opt/libomp)
OMP_OPTS_mac = -Xpreprocessor -fopenmp -I$(OMP_PREFIX)/include
OMP_LIBS_mac = -L$(OMP_PREFIX)/lib -lomp
OMP_OPTS_linux = -fopenmp
OMP_LIBS_linux = -fopenmp

#Options depending on the compiler/system:
OPTS_mac = -g -O3 -w -std=c99 $(OMP_OPTS_mac)
LIBS_mac = -lm -lc $(OMP_LIBS_mac)

OPTS_linux = -g -O3 -w -std=c99 $(OMP_OPTS_linux) 
LIBS_linux = -lm -lc $(OMP_LIBS_linux)

OPTS_AIX_RS6000 = -g -O3 #-Q -qsrcmsg -qmaxmem=4000 #-v
OPTS_SUN = -xO5 #-fast #-O3 -xO5 -g -fast -W
#CHOSE YOUR SYSTEM:
OPTS	= $(OPTS_mac)
LIBS	= $(LIBS_mac)
#------------------------------------------------------------------------
