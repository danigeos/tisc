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
tisc_DIR = \"`cd ..; pwd`\"
HOME_DIR = `cd; pwd`
OS 	= `../script/match_OS`
ARCH 	= `../script/match_arch`
BIN	= ../bin
OBJ	= .
LIB	= ../lib
INCL	= ../include 

#Other things:
BACKUPS	= $(HOME_DIR)/var/backups
ifneq ($(findstring THIN_SHEET,$(DEFS)),) 
	TS_LIBS = $(LIB)/thin_sheet.o $(LIB)/sistbanda.o 
	LC = $(FC)
else
	LC = $(CC)
endif

ifneq ($(findstring SURFACE_TRANSPORT,$(DEFS)),) 
	SP_LIBS = surf_proc.o 
endif

all:	clean tisc cuthrz lonlat2xy invertfile outin blockxyz call_tisc.o #projbanda
	exit

backup:
	cp -f	tisc*	$(BACKUPS)

clean:
	rm -f $(OBJ)/*.o $(LIB)/*.o $(LIB)/*.a


blockxyz: blockxyz.c 
	$(CC) blockxyz.c $(OPTS) $(LIB)/libreria.o $(LIBS) -I$(INCL) -o $(BIN)/blockxyz 
	chmod a+x $(BIN)/blockxyz

cuthrz: cuthrz.c tisclib.c tiscio.c libreria
	$(CC) cuthrz.c -o $(BIN)/cuthrz $(OPTS) $(LIBS) $(LIB)/libreria.o $(SP_LIBS) \
		-DTISCDIR=$(tisc_DIR) -I$(INCL)


invertfile: invertfile.c 
	$(CC) invertfile.c $(OPTS) $(LIBS) -I$(INCL) -o $(BIN)/invertfile
	chmod a+x $(BIN)/invertfile

outin: outin.c 
	$(CC) outin.c $(OPTS) $(LIB)/libreria.o $(LIBS) -I$(INCL) -o $(BIN)/outin 
	chmod a+x $(BIN)/outin


lonlat2xy: lonlat2xy.c libreria
	$(CC) lonlat2xy.c -o $(BIN)/lonlat2xy $(OPTS) $(LIBS) $(LIB)/libreria.o -I$(INCL)

projbanda:
	if [ $(findstring THIN_SHEET,$(DEFS)) ]; then \
		$(FC) projbanda.f -o $(BIN)/projbanda -O3 -I$(INCL) -L../lib/
		chmod  a+x  $(BIN)/$@
	fi

%:	
	$(CC) $*.c -c $(OPTS) -I$(INCL) 
	$(CC) $@.o -o $(BIN)/$@ $(OPTS) $(LIBS) $(LIB)/libreria.o
	chmod  a+x  $(BIN)/$@

$(SP_LIBS) : surf_proc.c 
	$(CC) -c -o surf_proc.o $(OPTS) -I$(INCL) surf_proc.c

tisc:	tisc.o libreria $(SP_LIBS) thin_sheet notify
	$(LC)  tisc.o $(LIB)/libreria.o $(SP_LIBS)  \
		$(LIBS) $(TS_LIBS) -o $(BIN)/tisc
	chmod a+x $(BIN)/tisc


tisc.o: tisc.c tisclib.c tiscio.c 
	@echo; echo tisc HOME DIRECTORY: $(tisc_DIR);  echo ARCH=$(ARCH);  
	@echo OP. SYSTEM=$(OS);  echo User home: $(HOME_DIR); 
	@echo You are currently $(shell whoami)	using the shell $(SHELL) ; echo
	$(CC) $(OPTS) -o tisc.o -DTISCDIR=$(tisc_DIR) -I$(INCL) -c tisc.c $(DEFS)
	if [ -f /usr/local/bin/ctags ]; then \
		/usr/local/bin/ctags *.c ../lib/thin_sheet.f ../lib/libreria.c ../include/*.h ;\
		echo tags file produced.;\
	fi

call_tisc.o : call_tisc.c
	if [ $(findstring SURFACE_TRANSPORT,$(DEFS)) ]; then \
		$(CC) -c -o call_tisc.o $(OPTS) -I$(INCL) -DSURFACE_TRANSPORT call_tisc.c ;\
	fi

libreria :
	(cd $(LIB)  ; make libreria.o ; )

thin_sheet :
	if [ $(findstring THIN_SHEET,$(DEFS)) ]; then \
		cd $(LIB); make thin_sheet; \
	fi

notify :
	@if [ -f ../doc/.first_compilation.txt ]; \
	then echo ; else echo First compilation of TISC in `hostname` by `whoami`>../doc/.first_compilation.txt; pwd>>../doc/.first_compilation.txt; \
		date>>../doc/.first_compilation.txt ; echo>> ../doc/.first_compilation.txt; finger `whoami` >> ../doc/.first_compilation.txt; mail danielgc@ictja.csic.es<../doc/.first_compilation.txt; \
	fi

