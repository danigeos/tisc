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
	TS_LIBS = $(OBJ)/thin_sheet.o $(OBJ)/sistbanda.o 
	LC = $(FC)
else
	LC = $(CC)
endif

ifneq ($(findstring SURFACE_TRANSPORT,$(DEFS)),) 
	SP_LIBS = $(OBJ)/surf_proc.o 
endif

# --- Fancy output & verbosity control ---
V ?= 0
ifeq ($(V),0)
  Q = @
else
  Q =
endif

BLUE  = \033[34m
GREEN = \033[32m
RESET = \033[0m

.PHONY: all backup clean notify help

all: ## Compile all binaries
	@$(MAKE) tisc cuthrz lonlat2xy outin blockxyz $(OBJ)/call_tisc.o

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
	@printf "$(BLUE)[CC]$(RESET) $<\n"
	$(Q)$(CC) $(CFLAGS) -D'VERSION="$(VERSION)"' -D'TISCDIR="$(tisc_DIR)"' -c $< -o $@

blockxyz: $(OBJ)/blockxyz.o $(OBJ)/globals.o $(OBJ)/libreria.o | $(OBJ)
	@printf "$(GREEN)[LD]$(RESET) $@\n"
	$(Q)$(CC) $(OBJ)/blockxyz.o $(OBJ)/globals.o $(OBJ)/libreria.o $(LIBS) -o $(BIN)/blockxyz 
	$(Q)chmod a+x $(BIN)/blockxyz

cuthrz: $(OBJ)/cuthrz.o $(TISC_OBJS) $(OBJ)/libreria.o | $(OBJ)
	@printf "$(GREEN)[LD]$(RESET) $@\n"
	$(Q)$(CC) $(OBJ)/cuthrz.o $(TISC_OBJS) $(OBJ)/libreria.o $(SP_LIBS) $(LIBS) -o $(BIN)/cuthrz \
		-D'TISCDIR="$(tisc_DIR)"'

outin: $(OBJ)/outin.o $(OBJ)/globals.o $(OBJ)/libreria.o | $(OBJ)
	@printf "$(GREEN)[LD]$(RESET) $@\n"
	$(Q)$(CC) $(OBJ)/outin.o $(OBJ)/globals.o $(OBJ)/libreria.o $(LIBS) -o $(BIN)/outin 
	$(Q)chmod a+x $(BIN)/outin

lonlat2xy: $(OBJ)/lonlat2xy.o $(OBJ)/globals.o $(OBJ)/libreria.o | $(OBJ)
	@printf "$(GREEN)[LD]$(RESET) $@\n"
	$(Q)$(CC) $(OBJ)/lonlat2xy.o $(OBJ)/globals.o $(OBJ)/libreria.o $(LIBS) -o $(BIN)/lonlat2xy

projbanda:
	@if [ $(findstring THIN_SHEET,$(DEFS)) ]; then \
		printf "$(GREEN)[FC]$(RESET) $@\n"; \
		$(FC) projbanda.f -o $(BIN)/projbanda -O3 -I$(INCL) -L../lib/ ; \
		chmod  a+x  $(BIN)/$@
	fi

$(OBJ)/surf_proc.o : surf_proc.c | $(OBJ)
	@printf "$(BLUE)[CC]$(RESET) $<\n"
	$(Q)$(CC) -c -o $@ $(CFLAGS) surf_proc.c

tisc:	$(OBJ)/tisc.o $(TISC_OBJS) $(OBJ)/libreria.o $(SP_LIBS) $(TS_LIBS) notify
	@printf "$(GREEN)[LD]$(RESET) $@\n"
	$(Q)$(LC)  $(OBJ)/tisc.o $(TISC_OBJS) $(OBJ)/libreria.o $(SP_LIBS)  \
		$(TS_LIBS) $(LIBS) -o $(BIN)/tisc
	$(Q)chmod a+x $(BIN)/tisc

$(OBJ)/tisc.o: tisc.c | $(OBJ)
	@printf "$(BLUE)[CC]$(RESET) $<\n"
	$(Q)$(CC) $(CFLAGS) -o $@ -D'VERSION="$(VERSION)"' -D'TISCDIR="$(tisc_DIR)"' -c tisc.c
	$(Q)if [ -x /usr/local/bin/ctags ]; then \
		/usr/local/bin/ctags *.c ../lib/thin_sheet.f ../lib/libreria.c ../include/*.h 2>/dev/null;\
		echo tags file produced.;\
	fi

$(OBJ)/call_tisc.o : call_tisc.c | $(OBJ)
	@if [ $(findstring SURFACE_TRANSPORT,$(DEFS)) ]; then \
		printf "$(BLUE)[CC]$(RESET) $<\n"; \
		$(CC) -c -o $@ $(CFLAGS) -DSURFACE_TRANSPORT call_tisc.c ; \
	fi

$(OBJ)/libreria.o : ../lib/libreria.c | $(OBJ)
	@printf "$(BLUE)[CC]$(RESET) $<\n"
	$(Q)$(CC) $(CFLAGS) -c $< -o $@

$(OBJ)/thin_sheet.o : ../lib/thin_sheet.f | $(OBJ)
	@printf "$(GREEN)[FC]$(RESET) $<\n"
	$(Q)$(FC) -O3 -c $< -o $@

$(OBJ)/sistbanda.o : ../lib/sistbanda.f | $(OBJ)
	@printf "$(GREEN)[FC]$(RESET) $<\n"
	$(Q)$(FC) -O3 -c $< -o $@

notify :
	@if [ -f ../doc/.first_compilation.txt ]; \
	then echo ; else echo First compilation of TISC in `hostname` by `whoami`>../doc/.first_compilation.txt; pwd>>../doc/.first_compilation.txt; \
		date>>../doc/.first_compilation.txt ; echo>> ../doc/.first_compilation.txt; finger `whoami` >> ../doc/.first_compilation.txt; mail danielgc@ictja.csic.es<../doc/.first_compilation.txt; \
	fi

help: ## Show this help message
	@awk 'BEGIN {FS = ":.*##"; printf "\nUsage:\n  make \033[36m<target>\033[0m\n\nTargets:\n"} /^[a-zA-Z_-]+:.*?##/ { printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2 }' $(MAKEFILE_LIST)

-include $(wildcard $(OBJ)/*.d)
