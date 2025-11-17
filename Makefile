#---------------------------- TISC makefile ----------------------------
#First read and modify options in ./config.mk
#
#Type  'make'  in this directory to compile.
#
#TISC has been succesfully compiled with this Makefile in: 
#  macOS 11, macOS 16, linux, 
#Earlier versions were functional for:
#  IBM AIX Version 3.2 for IBM RISC 6000 workstations, Hewlett Packard Envizex. Sun Solaris OS5
#------------------------------------------------------------------------

include config.mk

all:
	@echo; echo; echo Compiling version $(VERSION)
	(cd src; make all)
	@echo; echo; echo Compilation succeeded!
	@(echo "ADD TO YOUR PATH: `pwd`/bin/  AND  `pwd`/script/")
	@(echo "ADD IN .cshrc:    setenv tisc_dir `pwd` ")
	@(echo "ADD IN .bashrc:   export tisc_dir=`pwd` ")

clean_for_tar:
	(cd src; make clean)
	(cd demo; rm -f `find . -name '*.all' -print`)
	(cd demo; rm -f `find . -name '*.bas' -print`)
	(cd demo; rm -f `find . -name '*.lak' -print`)
	(cd demo; rm -f `find . -name '*.tmp' -print`)
	(cd demo; rm -f `find . -name '*[0-9][0-9][2-9].jpg' -print`)
	(rm -f `find . -name core -print`)


vers: 	clean_for_tar
	echo "CLEANING for packing"
	rm -R -f tmp tisc_copy_for_upload
	mkdir tmp tmp/bin
	cp -R -L Makefile config.mk README demo doc include lib script src   tmp 
	rm -f tmp/doc/.first_compilation.txt
	if [ $(findstring THIN_SHEET,$(DEFS)) ]; then echo Including thin sheet stuff; else \
		echo Removing thin sheet stuff; \
		rm tmp/lib/*thin_sheet* ; \
		rm tmp/lib/sistbanda* ; \
	fi
	if [ $(findstring SURFACE_TRANSPORT,$(DEFS)) ]; then echo Including surface processes stuff; else \
		echo Removing surface processes stuff; \
		rm tmp/src/*surf_proc* ; \
	fi
	echo "PACKING"
	tar -chf tisc.tar tmp
	chmod og-r tisc.tar
	gzip -f tisc.tar
	touch tmp/bin/touch_something #needed by git add
	mv tmp tisc_copy_for_upload
	make upload


upload:
	echo "UPLOADING to github."
	cd tisc_copy_for_upload
	#For initialization:  
	#git init; git remote add tisc https://github.com/danigeos/tisc; git add Makefile README config.mk bin demo doc include lib script src; git rm --cached doc/.first_compilation.txt
	git commit -a -mTISC_newVersion
	git config http.postBuffer 524288000; git config http.maxRequestBuffer 100M; git config core.compression 0
	#add --force to pass by the remote version 
	git push -u -f tisc master

