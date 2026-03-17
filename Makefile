# --- Root Makefile for tAo and TISC Monorepo ---

BIN := $(CURDIR)/bin
BUILD := $(CURDIR)/build

.PHONY: all clean help tao tisc dirs

.DEFAULT_GOAL := all

all: dirs tao tisc ## Compile both tAo and TISC
	@echo ""
	@echo "======================================================================="
	@echo "To run the executables, add the following to your shell profile (e.g., .bashrc or .zshrc):"
	@echo "  export tisc_dir=\"$(CURDIR)\""
	@echo "  export PATH=\"$$PATH:$(BIN)\""
	@echo "======================================================================="

tao: ## Compile tAo
	@$(MAKE) -C tao

tisc: ## Compile TISC
	@$(MAKE) -C tisc

clean: ## Clean generated files for both projects
	@$(MAKE) -C tao clean
	@$(MAKE) -C tisc clean
	@rm -rf $(BIN) $(BUILD)
	@echo "Cleaned both tAo and TISC."

dirs:
	@mkdir -p $(BIN) $(BUILD)

help: ## Show this help message
	@awk 'BEGIN {FS = ":.*##"; printf "\nUsage:\n  make \033[34m<target>\033[0m\n\nTargets:\n"} /^[a-zA-Z_-]+:.*?##/ { printf "  \033[34m%-15s\033[0m %s\n", $$1, $$2 }' $(MAKEFILE_LIST)
	@echo ""