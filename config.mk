# --- Shared Configuration for tAo and TISC ---

TAO_VERSION  = tAo_2026-03-22_factorized_stable
TISC_VERSION = TISC_2026-03-22_factorized_stable

# Compilers
CC = gcc
FC = gfortran

# Specific TISC definitions
TISC_DEFS = -DSURFACE_TRANSPORT #-DTHIN_SHEET 

# --- OpenMP Parallelization (Shared) ---
OMP_PREFIX := $(shell /opt/homebrew/bin/brew --prefix libomp 2>/dev/null || brew --prefix libomp 2>/dev/null || echo /usr/local/opt/libomp)
OMP_OPTS_mac = -Xpreprocessor -fopenmp -I$(OMP_PREFIX)/include
OMP_LIBS_mac = -L$(OMP_PREFIX)/lib -lomp
OMP_OPTS_linux = -fopenmp
OMP_LIBS_linux = -fopenmp

# Base Options
OPTS_mac = -g -O3 -w -std=c99 $(OMP_OPTS_mac)
LIBS_mac = -lm -lc $(OMP_LIBS_mac)

OPTS_linux = -g -O3 -w -std=c99 $(OMP_OPTS_linux) 
LIBS_linux = -lm -lc $(OMP_LIBS_linux)

# ACTIVE CONFIGURATION: Auto-detect OS
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    OPTS = $(OPTS_mac)
    LIBS = $(LIBS_mac)
else
    OPTS = $(OPTS_linux)
    LIBS = $(LIBS_linux)
endif