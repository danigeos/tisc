#---------------------------- TISC makefile ----------------------------
#First read and modify options in ./config.mk
#Then type  'make'  in this directory to compile.
#------------------------------------------------------------------------

include config.mk

# --- Colors ---
BLUE  := \033[34m
GREEN := \033[32m
RESET := \033[0m

.PHONY: all clean clean_for_tar help

.DEFAULT_GOAL := all

all: ## Compile the project
	@printf "\nTISC HOME DIRECTORY: $(CURDIR)\nARCH=$(shell uname -m)\nOP. SYSTEM=$(shell uname -s)\nYou are currently $(shell whoami) using the shell $(SHELL)\n"
	@printf "\n$(BLUE)=== Compiling version $(VERSION) ===$(RESET)\n"
	@$(MAKE) -C src all
	@printf "\n$(GREEN)=== Compilation succeeded! ===$(RESET)\n"
	@echo "Updating version in doc/template.PRM"
	LANG=C sed -E -i.bak 's/^version[[:space:]]+.*/version\t\t$(VERSION)/' doc/template.PRM
	@echo "ADD TO YOUR PATH: $(CURDIR)/bin/ AND $(CURDIR)/script/"
	@echo "ADD IN .bashrc:   export tisc_dir=$(CURDIR) "

clean: ## Clean generated files for packaging
	$(MAKE) -C src clean

clean_for_tar:
	$(MAKE) -C src clean
	find demo -type f \( -name '*.all' -o -name '*.bas' -o -name '*.lak' -o -name '*.tmp' -o -name '*[0-9][0-9][2-9].jpg' \) -delete
	find . -type f -name 'core' -delete

help: ## Show this help message
	@awk 'BEGIN {FS = ":.*##"; printf "\nUsage:\n  make $(BLUE)<target>$(RESET)\n\nTargets:\n"} /^[a-zA-Z_-]+:.*?##/ { printf "  $(BLUE)%-15s$(RESET) %s\n", $$1, $$2 }' $(MAKEFILE_LIST)
