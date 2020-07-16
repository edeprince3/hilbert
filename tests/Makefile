# hilbert test suite

SHELL := /bin/bash

# add new tests here
subdirs := doci pp2rdm pccd v2rdm_doci v2rdm2 v2rdm3 v2rdm4 v2rdm5 v2rdm6 v2rdm8

all-tests := $(addsuffix .test, $(subdirs))

quick-tests := $(addsuffix .test, pccd)

.PHONY : test all %.test 

test: $(all-tests)

quick: $(quick-tests)

%.test : 
	@echo ""
	@echo "    $(basename $@):"
	@echo ""
	@cd $(basename $@); psi4 #-d
	@echo ""

