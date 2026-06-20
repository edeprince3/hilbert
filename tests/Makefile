# hilbert test suite

SHELL := /bin/bash

# add new tests here
subdirs := qed_response-1 qed_response-2 qed_response-3 qed_ccsd qed_dft qed_dft_gradient qed_tddft-1 qed_tddft-2 doci pp2rdm pccd v2rdm_doci mcpdft v2rdm2 v2rdm3 v2rdm4 v2rdm5 v2rdm6 jellium_scf

all-tests := $(addsuffix .test, $(subdirs))

quick-tests := $(addsuffix .test, pccd)

v2rdm-tests := $(addsuffix .test, v2rdm2 v2rdm3 v2rdm4 v2rdm5 v2rdm6 v2rdm8)

qed-tests := $(addsuffix .test, qed_ccsd qed_tddft-1 qed_tddft-2 qed_dft qed_dft_gradient qed_response-1 qed_response-2 qed_response-3)

qed-ta-tests := $(addsuffix .test, qed_ccsd_ta-1 qed_ccsd_ta-2)

.PHONY : test all %.test 

test: $(all-tests)

quick: $(quick-tests)

v2rdm: $(v2rdm-tests)

qed: $(qed-tests)

qed_ta: $(qed-ta-tests)

%.test : 
	@echo ""
	@echo "    $(basename $@):"
	@echo ""
	@cd $(basename $@); psi4 #-d
	@echo ""

