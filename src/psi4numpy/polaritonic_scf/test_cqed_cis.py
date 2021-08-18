import psi4
from helper_cs_cqed_cis import cs_cqed_cis
from psi4.driver.procrouting.response.scf_response import tdscf_excitations
import numpy as np


numpy_memory = 2


def test_cqed_cis_no_field():
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

def test_cs_cqed_cis_with_field():

    # options dict
    options_dict = {'basis': 'cc-pVDZ',
        'save_jk': True,
        'scf_type' : 'pk',
        'e_convergence' : 1e-12,
        'd_convergence' : 1e-12
    }

    molstr = """
    0 1
        O      0.000000000000   0.000000000000  -0.068516219320
        H      0.000000000000  -0.790689573744   0.543701060715
        H      0.000000000000   0.790689573744   0.543701060715
    no_reorient
    symmetry c1
    """

    lam = np.array([0., 0., 0.05])
    om = 2./27.21138602

    expected_eigenvalues = np.array([-76.016613491776,-75.943171858458,-75.696248394443,-75.634194182018,-75.611459823919])
    cqed_dict = cs_cqed_cis(lam, om, molstr, options_dict)

    computed_eigenvalues = cqed_dict['CQED-CIS ENERGY']+cqed_dict['CQED-RHF ENERGY']
    assert np.allclose(expected_eigenvalues, computed_eigenvalues[:5])


#if test_cqed_cis_no_field():
#    print("CQED_CIS WITH NO FIELD PASSED!\n")
#else:
#    print("CQED_CIS WITH NO FIELD FAILED!\n")
#
#if test_cs_cqed_cis_with_field():
#    print("CQED_CIS WITH FIELD PASSED!\n")
#else:
#    print("CQED_CIS WITH FIELD FAILED!\n")





