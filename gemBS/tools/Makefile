#============================================================================
# PROJECT: gemBS
# FILE: Makefile
# DATE: 01/01/2015
# AUTHOR(S): Marcos Fernandez Callejo <mfernandez@cnag.crg.eu>
# DESCRIPTION: Top level makefile
#============================================================================

# Definitions
ROOT_PATH=$(CURDIR)

FOLDER_BIN=./bin/

all: release

static:	setup
	$(MAKE) --directory=src static
	cd loki; ./configure;
	$(MAKE) --directory=loki static
	$(MAKE) --directory=filter_vcf static
	$(MAKE) --directory=bs_call static
	cd gem3-mapper; ./configure;
	$(MAKE) --directory=gem3-mapper static

release: setup	
	$(MAKE) --directory=src	
	cd loki; ./configure;      
	$(MAKE) --directory=loki
	$(MAKE) --directory=filter_vcf
	$(MAKE) --directory=bs_call
	cd gem3-mapper; ./configure;
	$(MAKE) --directory=gem3-mapper

debug: setup
	$(MAKE) --directory=src debug
	cd loki; ./configure
	$(MAKE) --directory=loki debug
	$(MAKE) --directory=filter_vcf debug
	$(MAKE) --directory=bs_call debug
	cd gem3-mapper; ./configure;
	$(MAKE) --directory=gem3-mapper debug

setup: 
	@mkdir -p $(FOLDER_BIN) 

clean:
	@rm -f $(FOLDER_BIN)/*
	$(MAKE) --directory=gem3-mapper clean
	$(MAKE) --directory=bs_call clean

