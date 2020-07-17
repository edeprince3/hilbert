# hilbert test suite

SHELL := /bin/bash

# add new tests here
subdirs := doci pp2rdm pccd v2rdm_doci v2rdm2 v2rdm3 v2rdm4 v2rdm5 v2rdm6 

all-tests := $(addsuffix .test, $(subdirs))

quick-tests := $(addsuffix .test, pccd)

v2rdm-tests := $(addsuffix .test, v2rdm2 v2rdm3 v2rdm4 v2rdm5 v2rdm6 v2rdm8)

.PHONY : test all %.test 

test: $(all-tests)

quick: $(quick-tests)

v2rdm: $(v2rdm-tests)

%.test : 
	@echo ""
	@echo "    $(basename $@):"
	@echo ""
	@cd $(basename $@); psi4 #-d
	@echo ""

