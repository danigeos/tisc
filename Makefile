#---------------------------- TISC makefile ----------------------------
##First read and modify options in ./config.mk
#
#Type  'make'  in this directory to compile.
#
#tisc has been succesfully compiled with this Makefile in: 
#  iOS 11, linux, 
#Earlier versions were functional in:
#  IBM AIX Version 3.2 for IBM RISC 6000 workstations, Hewlett Packard Envizex. Sun Solaris OS5
#------------------------------------------------------------------------

include config.mk

all:
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
	rm -R -f tisc tisc_version
	mkdir tisc tisc/bin
	cp -R -L Makefile config.mk README demo doc include lib script src    tisc 
	rm -f tisc/doc/.first_compilation.txt
	if [ $(findstring THIN_SHEET,$(DEFS)) ]; then echo Including thin sheet stuff; else \
		echo Removing thin sheet stuff; \
		rm tisc/lib/*thin_sheet* ; \
		rm tisc/lib/sistbanda* ; \
	fi
	if [ $(findstring SURFACE_TRANSPORT,$(DEFS)) ]; then echo Including surface processes stuff; else \
		echo Removing surface processes stuff; \
		rm tisc/src/*surf_proc* ; \
	fi
	tar -chf tisc.tar tisc
	chmod og-r tisc.tar
	gzip -f tisc.tar
	echo "UPLOADING to github."
	touch tisc/bin/touch_something #needed by git add
	mv tisc tisc_copy_for_upload
	make upload


upload:
	cd tisc_copy_for_upload
	#For initialization:  
	#git init; git remote add tisc https://github.com/danigeos/tisc; git add Makefile README config.mk bin demo doc include lib script src; git rm --cached doc/.first_compilation.txt
	git commit -a -mTISC_newVersion
	git config http.postBuffer 524288000; git config http.maxRequestBuffer 100M; git config core.compression 0
	#add --force to pass by the remote version 
	git push -u -f tisc master

