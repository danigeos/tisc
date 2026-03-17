#----------------------------- tao makefile -----------------------------
#
#First read and modify options in ./config.mk
#
#Type  'make'  in this directory to compile tao 
#
#tAo has been succesfully compiled with this Makefile in: 
#  macOS 11, macOS 16, linux, 
#Earlier versions were functional for:
#  IBM AIX Version 3.2 for IBM RISC 6000 workstations, Hewlett Packard Envizex. Sun Solaris OS5
#------------------------------------------------------------------------

include config.mk

# --- Colors ---
BLUE  := \033[34m
GREEN := \033[32m
RESET := \033[0m

.PHONY: all clean clean_for_tar help

.DEFAULT_GOAL := all

all: ## Compile the project
	@printf "\ntAo HOME DIRECTORY: $(CURDIR)\nARCH=$(shell uname -m)\nOP. SYSTEM=$(shell uname -s)\nYou are currently $(shell whoami) using the shell $(SHELL)\n"
	@printf "\n$(BLUE)=== Compiling version $(VERSION) ===$(RESET)\n"
	@$(MAKE) -C src all
	@printf "\n$(GREEN)=== Compilation succeeded! ===$(RESET)\n"
	@echo "ADD TO YOUR PATH: $(CURDIR)/bin/ AND $(CURDIR)/script/"
	@echo "ADD IN .bashrc:   export tao_dir=$(CURDIR) "

clean: ## Clean generated files for packaging
	$(MAKE) -C src clean

clean_for_tar:
	$(MAKE) -C src clean
	find demo -type f \( -name '*[0-9][0-9][2-9].jpg' -o -name '*[0-9][0-9][2-9].png' \) -delete

help: ## Show this help message
	@awk 'BEGIN {FS = ":.*##"; printf "\nUsage:\n  make $(BLUE)<target>$(RESET)\n\nTargets:\n"} /^[a-zA-Z_-]+:.*?##/ { printf "  $(BLUE)%-15s$(RESET) %s\n", $$1, $$2 }' $(MAKEFILE_LIST)
