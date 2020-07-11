import pytest


@pytest.mark.quick

def test_doci():
    """hilbert/tests/doci"""
    #! 6-31g H2O DOCI Test 

    import psi4

    import sys
    sys.path.insert(0, '../..')
    import hilbert

    h2o = psi4.geometry("""
    0 1
    o
    h 1 1.0
    h 1 1.0 2 104.5
    symmetry c1
    """)

    interloper = psi4.geometry("""
    0 1
    O
    H 1 1.0
    H 1 1.0 2 90.0
    """)

    psi4.set_options({
      'basis': '6-31g',
      'scf_type': 'df',
      'e_convergence': 1e-10,
      'r_convergence': 1e-8,
      'orbopt_gradient_convergence': 1e-6,
      'orbopt_energy_convergence': 1e-8,
      'maxiter': 500,
    })
    psi4.set_module_options('hilbert', {
      'orbopt_maxiter': 1,
      'localize_orbitals': True,
      'noisy_orbitals': True,
      'optimize_orbitals': True,
    })

    psi4.activate(h2o)

    ref_scf  = -75.98014193580194
    ref_doci = -76.052851232676

    # be sure to save three-index integrals after scf
    psi4.core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    # get scf wfn
    scf_energy,ref_wfn = psi4.energy('scf',return_wfn=True)

    # grab options object
    options = psi4.core.get_options()
    options.set_current_module('HILBERT')

    # evaluate doci energy
    doci = hilbert.DOCIHelper(ref_wfn,options)
    current_energy = doci.compute_energy()

    assert psi4.compare_values(ref_scf, scf_energy, 8, "SCF total energy")
    assert psi4.compare_values(ref_doci, current_energy, 5, "DOCI total energy")

def test_pp2rdm():
    """hilbert/tests/pp2rdm"""
    #! 6-31g H2O pp2RDM Test 

    import psi4

    import sys
    sys.path.insert(0, '../..')
    import hilbert

    h2o = psi4.geometry("""
    0 1
    o
    h 1 1.0
    h 1 1.0 2 104.5
    symmetry c1
    """)

    interloper = psi4.geometry("""
    0 1
    O
    H 1 1.0
    H 1 1.0 2 90.0
    """)

    psi4.set_options({
      'basis': '6-31g',
      'scf_type': 'df',
      'e_convergence': 1e-10,
      'r_convergence': 1e-8,
      'orbopt_gradient_convergence': 1e-6,
      'orbopt_energy_convergence': 1e-8,
      'maxiter': 500,
    })
    psi4.set_module_options('hilbert', {
      'orbopt_maxiter': 1,
      'localize_orbitals': True,
      'noisy_orbitals': True,
      'optimize_orbitals': True,
      'p2rdm_type': 'k',
    })

    psi4.activate(h2o)

    ref_scf  = -75.98014193580194
    ref_pp2rdm = -76.052804542359

    # be sure to save three-index integrals after scf
    psi4.core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    # get scf wfn
    scf_energy,ref_wfn = psi4.energy('scf',return_wfn=True)

    # grab options object
    options = psi4.core.get_options()
    options.set_current_module('HILBERT')

    # evaluate pp2RDM energy
    pp2rdm = hilbert.pp2RDMHelper(ref_wfn,options)
    current_energy = pp2rdm.compute_energy()

    assert psi4.compare_values(ref_scf, scf_energy, 8, "SCF total energy")
    assert psi4.compare_values(ref_pp2rdm, current_energy, 5, "pp2RDM total energy")

def test_pccd():
    """hilbert/tests/pccd"""
    #! 6-31g H2O pCCD Test 

    import psi4

    import sys
    sys.path.insert(0, '../..')
    import hilbert

    h2o = psi4.geometry("""
    0 1
    o
    h 1 1.0
    h 1 1.0 2 104.5
    symmetry c1
    """)

    interloper = psi4.geometry("""
    0 1
    O
    H 1 1.0
    H 1 1.0 2 90.0
    """)

    psi4.set_options({
      'basis': '6-31g',
      'scf_type': 'df',
      'e_convergence': 1e-10,
      'r_convergence': 1e-8,
      'orbopt_gradient_convergence': 1e-6,
      'orbopt_energy_convergence': 1e-8,
      'maxiter': 500,
    })
    psi4.set_module_options('hilbert', {
      'orbopt_maxiter': 1,
      'localize_orbitals': True,
      'noisy_orbitals': True,
      'optimize_orbitals': True,
      'p2rdm_type': 'ccd',
    })

    psi4.activate(h2o)

    ref_scf  = -75.98014193580194
    ref_pccd = -76.052839394329

    # be sure to save three-index integrals after scf
    psi4.core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    # get scf wfn
    scf_energy,ref_wfn = psi4.energy('scf',return_wfn=True)

    # grab options object
    options = psi4.core.get_options()
    options.set_current_module('HILBERT')

    # evaluate pp2RDM energy
    pccd = hilbert.pp2RDMHelper(ref_wfn,options)
    current_energy = pccd.compute_energy()

    assert psi4.compare_values(ref_scf, scf_energy, 8, "SCF total energy")
    assert psi4.compare_values(ref_pccd, current_energy, 5, "pCCD total energy")

def test_v2rdm_doci():
    """hilbert/tests/v2rdm_doci"""
    #! 6-31g H2O v2RDM-DOCI Test 

    import psi4

    import sys
    sys.path.insert(0, '../..')
    import hilbert

    h2o = psi4.geometry("""
    0 1
    o
    h 1 1.0
    h 1 1.0 2 104.5
    symmetry c1
    """)

    interloper = psi4.geometry("""
    0 1
    O
    H 1 1.0
    H 1 1.0 2 90.0
    """)

    psi4.set_options({
      'basis': '6-31g',
      'scf_type': 'df',
      'e_convergence': 1e-4,
      'r_convergence': 1e-6,
      'maxiter': 100000,
    })
    psi4.set_module_options('hilbert', {
      'orbopt_maxiter': 20,
      'localize_orbitals': True,
      'noisy_orbitals': True,
      'optimize_orbitals': True,
      'orbopt_gradient_convergence': 1e-6,
      'orbopt_energy_convergence': 1e-8,
      'orbopt_frequency': 1000,
    })

    psi4.activate(h2o)

    ref_scf        = -75.98014193580194
    ref_v2rdm_doci = -76.053681341564

    # be sure to save three-index integrals after scf
    psi4.core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    # get scf wfn
    scf_energy,ref_wfn = psi4.energy('scf',return_wfn=True)

    # grab options object
    options = psi4.core.get_options()
    options.set_current_module('HILBERT')

    # evaluate v2RDM-DOCI energy
    v2rdm_doci = hilbert.v2RDMHelper(ref_wfn,options)
    current_energy = v2rdm_doci.compute_energy()

    assert psi4.compare_values(ref_scf, scf_energy, 8, "SCF total energy")
    assert psi4.compare_values(ref_v2rdm_doci, current_energy, 5, "v2RDM-DOCI total energy")

#test_doci()
#test_pp2rdm()
#test_pccd()
test_v2rdm_doci()
