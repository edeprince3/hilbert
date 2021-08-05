
# psi4 / hilbert
import psi4
import sys
sys.path.insert(0, '../../..')
import hilbert

# pyscf stuff
import pyscf
import openfermion as of
from openfermion.chem.molecular_data import spinorb_from_spatial
from openfermionpyscf import run_pyscf
from pyscf.cc.addons import spatial2spin
from pyscf import cc
import numpy as np
from functools import reduce
from pyscf import ao2mo

def v2rdm_psi4():

    mol = psi4.geometry("""
    0 1
    b
    h 1 1.1
    """)

    psi4.set_options({
      'basis':           'sto-3g',
      'scf_type':        'out_of_core',
      'd_convergence':   1e-10,
      'maxiter':         500,
    })
    psi4.set_module_options('hilbert', {
      'positivity':      'dqg',
      'r_convergence':   1e-5,
      'e_convergence':   1e-4,
      'maxiter':         20000,
    })

    psi4.activate(mol)

    # get scf wfn
    scf_energy,ref_wfn = psi4.energy('scf',return_wfn=True)

    # grab options object
    options = psi4.core.get_options()
    options.set_current_module('HILBERT')

    v2rdm_psi4 = hilbert.v2RDMHelper(ref_wfn,options)
    current_energy = v2rdm_psi4.compute_energy()

    return current_energy, mol.nuclear_repulsion_energy()

def v2rdm_pyscf():

    # pyscf reference
    np.set_printoptions(linewidth=500)
    basis = 'STO-3G'
    mol = pyscf.M(
        atom='B 0 0 0; H 0 0 {}'.format(1.1),
        basis=basis)

    mf = mol.RHF()
    mf.verbose = 3
    mf.run()

    # make h1 and spatial integrals in MO basis
    eri = ao2mo.kernel(mol, mf.mo_coeff)
    eri = ao2mo.restore(1, eri, mf.mo_coeff.shape[1])

    # this produces spatial MO h1 integrals
    h1 = reduce(np.dot, (mf.mo_coeff.T, mf.get_hcore(), mf.mo_coeff))

    # these parameters are hard-coded for now
    nalpha = 3
    nbeta  = 3
    nmo    = 6

    # hilbert options
    psi4.set_module_options('hilbert', {
      'positivity':      'dqg',
      'r_convergence':   1e-5,
      'e_convergence':   1e-4,
      'maxiter':         20000,
    })

    # grab options object
    options = psi4.core.get_options()
    options.set_current_module('HILBERT')

    v2rdm_pyscf = hilbert.v2RDMHelper(nalpha,nbeta,nmo,h1.flatten(),eri.flatten(),options)
    current_energy = v2rdm_pyscf.compute_energy()

    return current_energy


if __name__ == "__main__":

    energy_psi4, enuc  = v2rdm_psi4()
    energy_pyscf       = v2rdm_pyscf()

    print()
    print("  ==> Psi4 vs PySCF <== ")
    print()
    print("    v2RDM with Psi4 integrals:  %20.12f" % (energy_psi4))
    print("    v2RDM with PySCF integrals: %20.12f" % (energy_pyscf + enuc))
    print()

