#---------------------------- tisc makefile ----------------------------
#
#Read and modify options in ../config.mk
#
#Type  'make'  in this directory to compile.
#
#------------------------------------------------------------------------

include ../config.mk


#------------------------------------------------------------------------
#You should not modify these directories
tisc_DIR = $(abspath ..)
HOME_DIR = $(HOME)
OS 	= $(shell uname -s)
ARCH 	= $(shell uname -m)
BIN	= ../bin
OBJ	= ../build
LIB	= ../lib
INCL	= ../include 
MACOSX_DEPLOYMENT_TARGET = 10.9

#Other things:
BACKUPS	= $(HOME_DIR)/var/backups
ifneq ($(findstring THIN_SHEET,$(DEFS)),) 
	TS_LIBS = $(LIB)/thin_sheet.o $(LIB)/sistbanda.o 
	LC = $(FC)
else
	LC = $(CC)
endif

ifneq ($(findstring SURFACE_TRANSPORT,$(DEFS)),) 
	SP_LIBS = $(OBJ)/surf_proc.o 
endif

.PHONY: all backup clean libreria thin_sheet notify help

all: ## Compile all binaries
	$(MAKE) tisc cuthrz lonlat2xy outin blockxyz $(OBJ)/call_tisc.o

backup: ## Copy source files to backup directory
	$(CP) -f tisc* $(BACKUPS)

clean: ## Remove generated object files and binaries
	$(RM) -r $(OBJ)/*.o $(OBJ)/*.d $(LIB)/*.o $(LIB)/*.a $(BIN)/.dSYM 

DEPFLAGS = -MMD -MP
CFLAGS = $(OPTS) -I$(INCL) $(DEFS) $(DEPFLAGS)

TISC_OBJS = $(OBJ)/tisclib.o $(OBJ)/tiscio.o $(OBJ)/globals.o

$(OBJ):
	mkdir -p $(OBJ)

$(OBJ)/%.o: %.c | $(OBJ)
	$(CC) $(CFLAGS) -D'VERSION="$(VERSION)"' -D'TISCDIR="$(tisc_DIR)"' -c $< -o $@

blockxyz: $(OBJ)/blockxyz.o $(OBJ)/globals.o libreria | $(OBJ)
	$(CC) $(OBJ)/blockxyz.o $(OBJ)/globals.o $(LIB)/libreria.o $(LIBS) -o $(BIN)/blockxyz 
	chmod a+x $(BIN)/blockxyz

cuthrz: $(OBJ)/cuthrz.o $(TISC_OBJS) libreria | $(OBJ)
	$(CC) $(OBJ)/cuthrz.o $(TISC_OBJS) $(LIB)/libreria.o $(SP_LIBS) $(LIBS) -o $(BIN)/cuthrz \
		-D'TISCDIR="$(tisc_DIR)"'

outin: $(OBJ)/outin.o $(OBJ)/globals.o libreria | $(OBJ)
	$(CC) $(OBJ)/outin.o $(OBJ)/globals.o $(LIB)/libreria.o $(LIBS) -o $(BIN)/outin 
	chmod a+x $(BIN)/outin

lonlat2xy: $(OBJ)/lonlat2xy.o $(OBJ)/globals.o libreria | $(OBJ)
	$(CC) $(OBJ)/lonlat2xy.o $(OBJ)/globals.o $(LIB)/libreria.o $(LIBS) -o $(BIN)/lonlat2xy

projbanda:
	if [ $(findstring THIN_SHEET,$(DEFS)) ]; then \
		$(FC) projbanda.f -o $(BIN)/projbanda -O3 -I$(INCL) -L../lib/ ; \
		chmod  a+x  $(BIN)/$@
	fi

$(OBJ)/surf_proc.o : surf_proc.c | $(OBJ)
	$(CC) -c -o $@ $(CFLAGS) surf_proc.c

tisc:	$(OBJ)/tisc.o $(TISC_OBJS) libreria $(SP_LIBS) thin_sheet notify
	$(LC)  $(OBJ)/tisc.o $(TISC_OBJS) $(LIB)/libreria.o $(SP_LIBS)  \
		$(TS_LIBS) $(LIBS) -o $(BIN)/tisc
	chmod a+x $(BIN)/tisc

$(OBJ)/tisc.o: tisc.c | $(OBJ)
	@printf "\ntisc HOME DIRECTORY: $(tisc_DIR)\nARCH=$(ARCH)\nOP. SYSTEM=$(OS)\nUser home: $(HOME_DIR)\nYou are currently $(shell whoami) using the shell $(SHELL)\n\n"
	$(CC) $(CFLAGS) -o $@ -D'VERSION="$(VERSION)"' -D'TISCDIR="$(tisc_DIR)"' -c tisc.c
	@if [ -x /usr/local/bin/ctags ]; then \
		/usr/local/bin/ctags *.c ../lib/thin_sheet.f ../lib/libreria.c ../include/*.h 2>/dev/null;\
		echo tags file produced.;\
	fi

$(OBJ)/call_tisc.o : call_tisc.c | $(OBJ)
	if [ $(findstring SURFACE_TRANSPORT,$(DEFS)) ]; then \
		$(CC) -c -o $@ $(CFLAGS) -DSURFACE_TRANSPORT call_tisc.c ;\
	fi

libreria :
	$(CC) -c -o $(LIB)/libreria.o $(CFLAGS) $(LIB)/libreria.c

thin_sheet :
	if [ $(findstring THIN_SHEET,$(DEFS)) ]; then \
		$(MAKE) -C $(LIB) thin_sheet; \
	fi

notify :
	@if [ -f ../doc/.first_compilation.txt ]; \
	then echo ; else echo First compilation of TISC in `hostname` by `whoami`>../doc/.first_compilation.txt; pwd>>../doc/.first_compilation.txt; \
		date>>../doc/.first_compilation.txt ; echo>> ../doc/.first_compilation.txt; finger `whoami` >> ../doc/.first_compilation.txt; mail danielgc@ictja.csic.es<../doc/.first_compilation.txt; \
	fi

help: ## Show this help message
	@awk 'BEGIN {FS = ":.*##"; printf "\nUsage:\n  make \033[36m<target>\033[0m\n\nTargets:\n"} /^[a-zA-Z_-]+:.*?##/ { printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2 }' $(MAKEFILE_LIST)

-include $(wildcard $(OBJ)/*.d)
