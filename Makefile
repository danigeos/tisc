# --- Root Makefile for tAo and TISC Monorepo ---

tisc_DIR := $(CURDIR)
export tisc_DIR

BIN := $(tisc_DIR)/bin
BUILD := $(tisc_DIR)/tao+tisc_commons/build

include config.mk

OS        := $(shell uname -s)
ARCH      := $(shell uname -m)
PROCESSOR := $(shell uname -p)
USERNAME  := $(shell whoami)

BOLD    := \033[1m
CYAN    := \033[36m
BLUE    := \033[34m
GREEN   := \033[32m
YELLOW  := \033[33m
MAGENTA := \033[35m
RESET   := \033[0m

.PHONY: all clean help tao tisc dirs info update_prm_versions

.DEFAULT_GOAL := all

all: info dirs update_prm_versions tao tisc ## Compile both tAo and TISC
	@printf "\n"
	@printf "$(MAGENTA)$(BOLD)============================================================$(RESET)\n"
	@printf "$(GREEN)$(BOLD)      🎉 SUCCESS! TISC and tAo are ready to use! 🎉$(RESET)\n"
	@printf "$(MAGENTA)$(BOLD)============================================================$(RESET)\n"
	@echo ""
	@echo "======================================================================="
	@echo "To run the executables, add the following to your shell profile (e.g., .bashrc or .zshrc):"
	@echo "  export tisc_dir=\"$(tisc_DIR)/\""
	@echo "  export PATH=\"$(BIN):$(tisc_DIR)/tisc/script/:$(tisc_DIR)/tao/script/:\$$PATH:\""
	@echo "======================================================================="

info:
	@printf "\n"
	@printf "$(MAGENTA)$(BOLD)============================================================$(RESET)\n"
	@printf "$(YELLOW)$(BOLD)      🌟 COMPILING TISC (Tectonics, Isostasy, Surface) 🌟$(RESET)\n"
	@printf "$(MAGENTA)$(BOLD)============================================================$(RESET)\n"
	@printf "$(CYAN)  User:$(RESET)        $(USERNAME)\n"
	@printf "$(CYAN)  OS:$(RESET)          $(OS)\n"
	@printf "$(CYAN)  Arch:$(RESET)        $(ARCH)\n"
	@printf "$(CYAN)  Processor:$(RESET)   $(PROCESSOR)\n"
	@printf "$(MAGENTA)$(BOLD)============================================================$(RESET)\n\n"

update_prm_versions:
	@if ! grep -q "^version[[:space:]]*$(TISC_VERSION)" tisc/doc/template.PRM; then \
		echo "Updating TISC version tag in template.PRM..."; \
		LC_ALL=C sed -i.bak -e "s/^version[[:space:]]*[^[:space:]]*/version		$(TISC_VERSION)/" tisc/doc/template.PRM; \
		rm -f tisc/doc/template.PRM.bak; \
	fi
	@if ! grep -q "^version[[:space:]]*$(TAO_VERSION)" tao/doc/template.PRM; then \
		echo "Updating tAo version tag in template.PRM..."; \
		LC_ALL=C sed -i.bak -e "s/^version[[:space:]]*[^[:space:]]*/version		$(TAO_VERSION)/" tao/doc/template.PRM; \
		rm -f tao/doc/template.PRM.bak; \
	fi

tao: update_prm_versions ## Compile tAo
	@$(MAKE) -C tao/src all

tisc: update_prm_versions ## Compile TISC
	@$(MAKE) -C tisc/src all

clean: ## Clean generated files for both projects
	@$(MAKE) -C tao/src clean
	@$(MAKE) -C tisc/src clean
	@rm -rf $(BIN) $(BUILD)
	@echo "Cleaned both tAo and TISC."

dirs:
	@mkdir -p $(BIN) $(BUILD)

help: ## Show this help message
	@awk 'BEGIN {FS = ":.*##"; printf "\nUsage:\n  make \033[34m<target>\033[0m\n\nTargets:\n"} /^[a-zA-Z_-]+:.*?##/ { printf "  \033[34m%-15s\033[0m %s\n", $$1, $$2 }' $(MAKEFILE_LIST)
	@echo ""
