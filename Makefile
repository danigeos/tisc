# --- Root Makefile for tAo and TISC Monorepo ---

tisc_DIR := $(CURDIR)
export tisc_DIR

BIN := $(tisc_DIR)/bin
BUILD := $(tisc_DIR)/tao+tisc_commons/build

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

.PHONY: all clean help tao tisc dirs info

.DEFAULT_GOAL := all

all: info dirs tao tisc ## Compile both tAo and TISC
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

tao: ## Compile tAo
	@$(MAKE) -C tao/src all

tisc: ## Compile TISC
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
