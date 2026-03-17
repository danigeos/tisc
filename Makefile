# --- Root Makefile for tAo and TISC Monorepo ---

.PHONY: all clean help tao tisc

.DEFAULT_GOAL := all

all: tao tisc ## Compile both tAo and TISC

tao: ## Compile tAo
	@$(MAKE) -C tao

tisc: ## Compile TISC
	@$(MAKE) -C tisc

clean: ## Clean generated files for both projects
	@$(MAKE) -C tao clean
	@$(MAKE) -C tisc clean
	@echo "Cleaned both tAo and TISC."

help: ## Show this help message
	@awk 'BEGIN {FS = ":.*##"; printf "\nUsage:\n  make \033[34m<target>\033[0m\n\nTargets:\n"} /^[a-zA-Z_-]+:.*?##/ { printf "  \033[34m%-15s\033[0m %s\n", $$1, $$2 }' $(MAKEFILE_LIST)
	@echo ""