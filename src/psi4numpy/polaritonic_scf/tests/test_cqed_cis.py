import psi4
import os
from ..helper_cs_cqed_cis import cs_cqed_cis
from psi4.driver.procrouting.response.scf_response import tdscf_excitations
import numpy as np
import pytest


numpy_memory = 2


def test_cqed_cis():
    # dictionary for psi4 options for test calculations
    psi4_options_dict = {
        'basis': 'cc-pVDZ',
        'save_jk': True, 
        'scf_type': 'pk'
    }
    mol_str = """
    Mg
    H 1 1.3
    symmetry c1
    1 1
    """
    psi4.set_options(psi4_options_dict)
    mol = psi4.geometry(mol_str)
    
    # electric field
    lambda_vector = np.array([0, 0, 0])
    omega = 0.
    
    n_states = 5
    # run psi4 SCF
    psi4_rhf_e, wfn = psi4.energy("scf/cc-pVDZ", return_wfn=True, molecule=mol)
    
    # calculate the excited-state energies and save them to a dictionary called 'res'
    res = tdscf_excitations(wfn, states=n_states, triplets = "NONE", tda=True)
    # parse res for excitation energies
    psi4_excitation_e = [r["EXCITATION ENERGY"] for r in res]

    mgh_dict = cs_cqed_cis(lambda_vector, omega, mol_str, psi4_options_dict)

    assert np.isclose(mgh_dict['RHF ENERGY'], psi4_rhf_e)
    assert np.isclose(mgh_dict['CQED-RHF ENERGY'], psi4_rhf_e)
    cqed_cis_e = mgh_dict['CQED-CIS ENERGY']
    print("printing comparison")
    assert np.isclose(cqed_cis_e[2], psi4_excitation_e[0])
    assert np.isclose(cqed_cis_e[4], psi4_excitation_e[1])
    assert np.isclose(cqed_cis_e[6], psi4_excitation_e[2])
    assert np.isclose(cqed_cis_e[8], psi4_excitation_e[3])
    assert np.isclose(cqed_cis_e[10], psi4_excitation_e[4])







