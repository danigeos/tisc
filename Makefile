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

.PHONY: all clean clean_for_tar help

.DEFAULT_GOAL := all

all: ## Compile the project
	@printf "\n\nCompiling version $(VERSION)\n"
	$(MAKE) -C src all
	@printf "\n\nCompilation succeeded!\n"
	@echo "Updating version in doc/template.PRM"
	LANG=C sed -E -i.bak 's/^version[[:space:]]+.*/version\t\t$(VERSION)/' doc/template.PRM
	@echo "ADD TO YOUR PATH: $(CURDIR)/bin/  AND  $(CURDIR)/script/"
	@echo "ADD IN .cshrc:    setenv tisc_dir $(CURDIR) "
	@echo "ADD IN .bashrc:   export tisc_dir=$(CURDIR) "

clean: ## Clean generated files for packaging
	$(MAKE) -C src clean

clean_for_tar:
	$(MAKE) -C src clean
	find demo -type f \( -name '*.all' -o -name '*.bas' -o -name '*.lak' -o -name '*.tmp' -o -name '*[0-9][0-9][2-9].jpg' \) -delete
	find . -type f -name 'core' -delete

help: ## Show this help message
	@awk 'BEGIN {FS = ":.*##"; printf "\nUsage:\n  make \033[36m<target>\033[0m\n\nTargets:\n"} /^[a-zA-Z_-]+:.*?##/ { printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2 }' $(MAKEFILE_LIST)
