import psi4
import os
from ..helper_cqed_rhf import cqed_rhf
import numpy as np
import pytest


numpy_memory = 2


def test_cqed_rhf():
    # dictionary for psi4 options for test calculations
    psi4_options_dict = {
                  'basis':        'def2-tzvppd',
                  'scf_type':     'pk',
                  'reference':    'rhf',
                  'mp2_type':     'conv',
                  'save_jk': True,
                  'e_convergence': 1e-10,
                  'd_convergence': 1e-10
    }
    psi4.set_options(psi4_options_dict)
    
    # strings to define test molecules
    NaF_string = """
    
    0 1
        NA           0.000000000000     0.000000000000    -0.875819904077
        F            0.000000000000     0.000000000000     1.059820520433
    no_reorient
    symmetry c1
    """

    NaCl_string = """
    
    0 1
        NA           0.000000000000     0.000000000000    -1.429419641344
        CL           0.000000000000     0.000000000000     0.939751385626
    no_reorient
    symmetry c1
    """
    
    expected_NaF =  -261.371070718358
    expected_NaCl = -621.438985539266
    # electric field
    lambda_vector = np.array([0, 0, 0.01])

    naf_dict = cqed_rhf(lambda_vector, NaF_string, psi4_options_dict)
    nacl_dict = cqed_rhf(lambda_vector, NaCl_string, psi4_options_dict)

    assert np.isclose(naf_dict['CQED-RHF ENERGY'], expected_NaF,5e-5)
    assert np.isclose(nacl_dict['CQED-RHF ENERGY'], expected_NaCl,5e-5)


