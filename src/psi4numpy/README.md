polaritonic_scf
==============================

A Psi4Numpy implementation of ab inition - cavity quantum electrodynamics methods 

### Copyright
Copyright (c) 2021, Foley Lab, William Paterson University

### Requirements
- python3.7 or greater
- numpy
- scipy
- pytest
- psi4

### Contents
- polaritonic_scf/helper_cqed_rhf.py -> helper function to compute CQED-RHF energy, wavefunction, and associated quantities in the coherent-state basis
- polaritonic_scf/helper_cs_cqed_cis.py -> helper function to compute CQED-CIS energy, wavefunction, and associated quantities in the coherent-state basis
- polaritonic_scf/test_cqed_rhf.py -> unit test of helper_cqed_rhf.py   
- polaritonic_scf/test_cqed_cis.py -> unit test of helper_cs_cqed_cis.py

To run unit tests: From this directory, 

`pytest -v`  



