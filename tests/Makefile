# hilbert test suite

SHELL := /bin/bash

# add new tests here
subdirs := doci pp2rdm pccd

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

