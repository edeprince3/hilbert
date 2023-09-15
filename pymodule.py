#
# @BEGIN LICENSE
#
# hilbert by Psi4 Developer, a plugin to:
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import numpy as np

import psi4
import psi4.driver.p4util as p4util
from psi4.driver.procrouting import proc_util
from psi4.driver.procrouting import proc
#from psi4.driver.qcdb import molecule

# to build a fake molecule
import qcelemental as qcel
from psi4.driver import qcdb

def run_polaritonic_scf(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    polaritonic scf can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('polaritonic-rhf')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    optstash = p4util.OptionsState(
        ['SCF', 'DF_INTS_IO'])

    psi4.core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    if ( lowername == 'polaritonic-rhf' ):
        psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_RHF')
    elif ( lowername == 'polaritonic-uhf' ):
        psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_UHF')
    elif ( lowername == 'polaritonic-rohf' ):
        psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_ROHF')
    elif ( lowername == 'polaritonic-uks' ):
        psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_UKS')
    elif ( lowername == 'polaritonic-rks' ):
        psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_RKS')
    elif ( lowername == 'polaritonic-rcis' ):
        psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_RCIS')
    elif ( lowername == 'polaritonic-uccsd' ):
        psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_UCCSD')
    elif ( lowername == 'polaritonic-tddft' ):
        psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_TDDFT')
    elif ( lowername == 'polaritonic-rpa' ):
        psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_RPA')

    # Compute a SCF reference, a wavefunction is return which holds the molecule used, orbitals
    # Fock matrices, and more
    #print('Attention! This SCF may be density-fitted.')
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = psi4.driver.scf_helper(name, **kwargs)

    scf_aux_basis = psi4.core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_SCF",
                                        psi4.core.get_option("SCF", "DF_BASIS_SCF"),
                                        "JKFIT", psi4.core.get_global_option('BASIS'),
                                        puream=ref_wfn.basisset().has_puream())
    ref_wfn.set_basisset("DF_BASIS_SCF", scf_aux_basis)

    aux_basis = psi4.core.BasisSet.build(ref_wfn.molecule(), "DF_BASIS_CC",
                                        psi4.core.get_global_option("DF_BASIS_CC"),
                                        "RIFIT", psi4.core.get_global_option("BASIS"))
    ref_wfn.set_basisset("DF_BASIS_CC", aux_basis)

    # Ensure IWL files have been written when not using DF/CD
    scf_type = psi4.core.get_option('SCF', 'SCF_TYPE')
    if ( scf_type == 'PK' or scf_type == 'DIRECT' ):
        proc_util.check_iwl_file_from_scf_type(psi4.core.get_option('SCF', 'SCF_TYPE'), ref_wfn)

    # Call the Psi4 plugin
    # Please note that setting the reference wavefunction in this way is ONLY for plugins
    rhf_wfn = psi4.core.plugin('hilbert.so', ref_wfn)

    optstash.restore()

    return rhf_wfn

def run_doci(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    doci can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('doci')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    optstash = p4util.OptionsState(
        ['SCF', 'DF_INTS_IO'])

    psi4.core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'DOCI')

    # Compute a SCF reference, a wavefunction is return which holds the molecule used, orbitals
    # Fock matrices, and more
    #print('Attention! This SCF may be density-fitted.')
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = psi4.driver.scf_helper(name, **kwargs)

    # Ensure IWL files have been written when not using DF/CD
    scf_type = psi4.core.get_option('SCF', 'SCF_TYPE')
    if ( scf_type == 'PK' or scf_type == 'DIRECT' ):
        proc_util.check_iwl_file_from_scf_type(psi4.core.get_option('SCF', 'SCF_TYPE'), ref_wfn)

    # Call the Psi4 plugin
    # Please note that setting the reference wavefunction in this way is ONLY for plugins
    doci_wfn = psi4.core.plugin('hilbert.so', ref_wfn)

    optstash.restore()

    return doci_wfn

def run_pp2rdm(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    pp2rdm can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('pp2rdm')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    optstash = p4util.OptionsState(
        ['SCF', 'DF_INTS_IO'])

    psi4.core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'PP2RDM')

    if lowername == 'pp2rdm':
        psi4.core.set_local_option('HILBERT', 'P2RDM_TYPE', 'K')
    elif lowername == 'pccd':
        psi4.core.set_local_option('HILBERT', 'P2RDM_TYPE', 'CCD')
    elif lowername == 'pcid':
        psi4.core.set_local_option('HILBERT', 'P2RDM_TYPE', 'CID')
    elif lowername == 'pcepa(1)':
        psi4.core.set_local_option('HILBERT', 'P2RDM_TYPE', 'CEPA(1)')
    elif lowername == 'pcepa(0)':
        psi4.core.set_local_option('HILBERT', 'P2RDM_TYPE', 'CEPA(0)')
    elif lowername == 'pacpf':
        psi4.core.set_local_option('HILBERT', 'P2RDM_TYPE', 'ACPF')
    elif lowername == 'paqcc':
        psi4.core.set_local_option('HILBERT', 'P2RDM_TYPE', 'AQCC')

    # Compute a SCF reference, a wavefunction is return which holds the molecule used, orbitals
    # Fock matrices, and more
    #print('Attention! This SCF may be density-fitted.')
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = psi4.driver.scf_helper(name, **kwargs)

    # Ensure IWL files have been written when not using DF/CD
    scf_type = psi4.core.get_option('SCF', 'SCF_TYPE')
    if ( scf_type == 'PK' or scf_type == 'DIRECT' ):
        proc_util.check_iwl_file_from_scf_type(psi4.core.get_option('SCF', 'SCF_TYPE'), ref_wfn)

    # Call the Psi4 plugin
    # Please note that setting the reference wavefunction in this way is ONLY for plugins
    pp2rdm_wfn = psi4.core.plugin('hilbert.so', ref_wfn)

    optstash.restore()

    return pp2rdm_wfn

def run_v2rdm_doci(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    v2rdm_doci can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('v2rdm_doci')

    """

    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    optstash = p4util.OptionsState(
        ['SCF', 'DF_INTS_IO'])

    psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'V2RDM_DOCI')

    psi4.core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    # Your plugin's psi4 run sequence goes here
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = psi4.driver.scf_helper(name, **kwargs)

    # if restarting from a checkpoint file, this file
    # needs to be in scratch with the correct name
    filename = psi4.core.get_option("HILBERT","RESTART_FROM_CHECKPOINT_FILE")

    # todo PSIF_V2RDM_CHECKPOINT should be definied in psifiles.h
    #if ( filename != "" and psi4.core.get_global_option('DERTYPE') != 'FIRST' ):
    #    molname = psi4.wavefunction().molecule().name()
    #    p4util.copy_file_to_scratch(filename,'psi',molname,269,False)

    # Ensure IWL files have been written when not using DF/CD
    scf_type = psi4.core.get_option('SCF', 'SCF_TYPE')
    if ( scf_type == 'PK' or scf_type == 'DIRECT' ):
        proc_util.check_iwl_file_from_scf_type(psi4.core.get_option('SCF', 'SCF_TYPE'), ref_wfn)

    v2rdm_doci_wfn = psi4.core.plugin('hilbert.so', ref_wfn)

    optstash.restore()

    return v2rdm_doci_wfn

def run_v2rdm_casscf(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    v2rdm_casscf can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('v2rdm_casscf')

    """

    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    optstash = p4util.OptionsState(
        ['SCF', 'DF_INTS_IO'])

    psi4.core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'V2RDM_CASSCF')

    # Your plugin's psi4 run sequence goes here
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = psi4.driver.scf_helper(name, **kwargs)

    # if restarting from a checkpoint file, this file
    # needs to be in scratch with the correct name
    filename = psi4.core.get_option("HILBERT","RESTART_FROM_CHECKPOINT_FILE")

    # Ensure IWL files have been written when not using DF/CD
    scf_type = psi4.core.get_option('SCF', 'SCF_TYPE')
    if ( scf_type == 'PK' or scf_type == 'DIRECT' ):
        proc_util.check_iwl_file_from_scf_type(psi4.core.get_option('SCF', 'SCF_TYPE'), ref_wfn)

    # reorder wavefuntions based on user input
    # apply a list of 2x2 rotation matrices to the orbitals in the form of [irrep, orbital1, orbital2, theta]
    # where an angle of 0 would do nothing and an angle of 90 would switch the two orbitals.
    # the indices of irreps and orbitals start from 0
    reorder_orbitals = psi4.core.get_option("HILBERT","MCSCF_ROTATE")
    for orbord in reorder_orbitals:
        if type(orbord) != list :
            raise psi4.p4util.PsiException("Each element of the orbtial rotate list requires 4 arguements (irrep, orb1, orb2, theta).")
        if len(orbord) != 4:
            raise psi4.p4util.PsiException("Each element of the orbtial rotate list requires 4 arguements (irrep, orb1, orb2, theta).")

        irrep, orb1, orb2, theta = orbord

        if irrep > ref_wfn.Ca().nirrep():
            raise psi4.p4util.PsiException("REORDER_ORBITALS: Expression %s irrep number is larger than the number of irreps" %
                                    (str(orbord)))

        if max(orb1, orb2) > ref_wfn.Ca().coldim()[irrep]:
            raise psi4.p4util.PsiException("REORDER_ORBITALS: Expression %s orbital number exceeds number of orbitals in irrep" %
                                    (str(orbord)))

        theta = np.deg2rad(theta)

        x_a = ref_wfn.Ca().nph[irrep][:, orb1].copy()
        y_a = ref_wfn.Ca().nph[irrep][:, orb2].copy()

        xp_a = np.cos(theta) * x_a - np.sin(theta) * y_a
        yp_a = np.sin(theta) * x_a + np.cos(theta) * y_a

        ref_wfn.Ca().nph[irrep][:, orb1] = xp_a
        ref_wfn.Ca().nph[irrep][:, orb2] = yp_a

        # beta orbitals are not changed because v2rdm does not use beta orbitals
        # if/when beta orbitals are used, you should change beta orbitals too,
        # but keep in mind that in RHF wavefunctions, beta orbital pointer may
        # point to the same location as alpha orbital pointer

    returnvalue = psi4.core.plugin('hilbert.so', ref_wfn)

    optstash.restore()

    return returnvalue

def run_v2rdm_casscf_gradient(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    v2rdm_casscf can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> gradient('v2rdm_casscf')

    """

    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    optstash = p4util.OptionsState(
        ['GLOBALS', 'DERTYPE'],
        ['HILBERT', 'OPTIMIZE_ORBITALS'],
        ['HILBERT', 'SEMICANONICALIZE_ORBITALS'],
        ['HILBERT', 'ORBOPT_ACTIVE_ACTIVE_ROTATIONS'],
        ['HILBERT', 'RESTART_FROM_CHECKPOINT_FILE'],
        ['HILBERT', 'WRITE_CHECKPOINT_FILE'])

    psi4.core.set_global_option('DERTYPE', 'FIRST')
    psi4.core.set_local_option("HILBERT","OPTIMIZE_ORBITALS",True)
    psi4.core.set_local_option("HILBERT","ORBOPT_ACTIVE_ACTIVE_ROTATIONS",True)
    psi4.core.set_local_option("HILBERT","SEMICANONICALIZE_ORBITALS",False)
    psi4.core.set_local_option("HILBERT","RESTART_FROM_CHECKPOINT_FILE","DUMMY")
    psi4.core.set_local_option("HILBERT","WRITE_CHECKPOINT_FILE",True)

    # analytic derivatives do not work with scf_type df/cd
    scf_type = psi4.core.get_option('SCF', 'SCF_TYPE')
    if ( scf_type == 'CD' or scf_type == 'DF' ):
        raise ValidationError("""Error: analytic v2RDM-CASSCF gradients not implemented for scf_type %s.""" % scf_type)

    v2rdm_wfn = run_v2rdm_casscf(name,**kwargs)
    derivobj = psi4.core.Deriv(v2rdm_wfn)
    derivobj.set_deriv_density_backtransformed(True)
    derivobj.set_ignore_reference(True)
    grad = derivobj.compute()

    v2rdm_wfn.set_gradient(grad)

    optstash.restore()

    return v2rdm_wfn

def run_p2rdm(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    p2rdm can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('p2rdm')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    optstash = p4util.OptionsState(
        ['SCF', 'DF_INTS_IO'])

    psi4.core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'P2RDM')

    #if lowername == 'p2rdm':
    #    psi4.core.set_local_option('HILBERT', 'P2RDM_TYPE', 'K')
    #elif lowername == 'cid':
    #    psi4.core.set_local_option('HILBERT', 'P2RDM_TYPE', 'CID')
    #elif lowername == 'cepa(0)':
    #    psi4.core.set_local_option('HILBERT', 'P2RDM_TYPE', 'CEPA(0)')

    # Compute a SCF reference, a wavefunction is return which holds the molecule used, orbitals
    # Fock matrices, and more
    #print('Attention! This SCF may be density-fitted.')
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = psi4.driver.scf_helper(name, **kwargs)

    # Ensure IWL files have been written when not using DF/CD
    scf_type = psi4.core.get_option('SCF', 'SCF_TYPE')
    if ( scf_type == 'PK' or scf_type == 'DIRECT' ):
        proc_util.check_iwl_file_from_scf_type(psi4.core.get_option('SCF', 'SCF_TYPE'), ref_wfn)

    # Call the Psi4 plugin
    # Please note that setting the reference wavefunction in this way is ONLY for plugins
    p2rdm_wfn = psi4.core.plugin('hilbert.so', ref_wfn)

    optstash.restore()

    return p2rdm_wfn

def run_jellium_scf(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    jellium_scf can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('jellium-scf')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # build empty reference wavefunction to pass into plugin

    #ref_molecule = kwargs.get('molecule', psi4.core.get_active_molecule())
    mol = """H 0 0 0
    H 0 0 1"""
    ref_molecule = psi4.core.Molecule.from_string(mol)
    base_wfn = psi4.core.Wavefunction.build(ref_molecule, 'STO-3G')
    ref_wfn = proc.scf_wavefunction_factory('HF', base_wfn, psi4.core.get_global_option('REFERENCE'))

    psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'JELLIUM_SCF')

    jellium_scf_wfn = psi4.core.plugin('hilbert.so', ref_wfn)

    return jellium_scf_wfn

def density_analysis(**kwargs):
    r"""Function to evaluate real-space density"""

    kwargs = p4util.kwargs_lower(kwargs)

    psi4.core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        raise ValidationError("""Error: density_analysis requires a reference wave function.""" )

    func = 'M06-2X'
    ref_molecule = kwargs.get('molecule', psi4.core.get_active_molecule())
    base_wfn = psi4.core.Wavefunction.build(ref_molecule, psi4.core.get_global_option('BASIS'))
    new_wfn = proc.scf_wavefunction_factory('M06-2X', base_wfn, 'UKS')

    # push reference orbitals onto new wave function 
    for irrep in range (0,ref_wfn.Ca().nirrep()):
        new_wfn.Ca().nph[irrep][:,:] = ref_wfn.Ca().nph[irrep][:,:]
        new_wfn.Cb().nph[irrep][:,:] = ref_wfn.Cb().nph[irrep][:,:]

    # push reference energies onto new wave function
    for irrep in range (0,ref_wfn.epsilon_a().nirrep()):
        new_wfn.epsilon_a().nph[irrep][:] = ref_wfn.epsilon_a().nph[irrep][:]
        new_wfn.epsilon_b().nph[irrep][:] = ref_wfn.epsilon_b().nph[irrep][:]

    # grab options object
    options = psi4.core.get_options()
    options.set_current_module('HILBERT')

    # evaluate doci energy
    # build real-space density
    import hilbert
    real_space_density = hilbert.RealSpaceDensityHelper(new_wfn,options)

    # Call the Psi4 plugin
    # Please note that setting the reference wavefunction in this way is ONLY for plugins
    #mcpdft_wfn = psi4.core.plugin('mcpdft.so', ref_wfn)

    return real_space_density

# Integration with driver routines

# jellium-scf
psi4.driver.procedures['energy']['jellium-scf'] = run_jellium_scf

# p2rdm
psi4.driver.procedures['energy']['p2rdm'] = run_p2rdm
#psi4.driver.procedures['energy']['cid']   = run_p2rdm

# pair methods:
psi4.driver.procedures['energy']['pp2rdm']   = run_pp2rdm
psi4.driver.procedures['energy']['pcid']     = run_pp2rdm
psi4.driver.procedures['energy']['pccd']     = run_pp2rdm
psi4.driver.procedures['energy']['pcepa(0)'] = run_pp2rdm
psi4.driver.procedures['energy']['pacpf']    = run_pp2rdm
psi4.driver.procedures['energy']['paqcc']    = run_pp2rdm

# doci
psi4.driver.procedures['energy']['doci'] = run_doci

# v2rdm-doci
psi4.driver.procedures['energy']['v2rdm-doci'] = run_v2rdm_doci

# v2rdm-casscf
psi4.driver.procedures['energy']['v2rdm-casscf'] = run_v2rdm_casscf
psi4.driver.procedures['gradient']['v2rdm-casscf'] = run_v2rdm_casscf_gradient

# polaritonic scf and cc
psi4.driver.procedures['energy']['polaritonic-rhf'] = run_polaritonic_scf
psi4.driver.procedures['energy']['polaritonic-uhf'] = run_polaritonic_scf
psi4.driver.procedures['energy']['polaritonic-rohf'] = run_polaritonic_scf
psi4.driver.procedures['energy']['polaritonic-uks'] = run_polaritonic_scf
psi4.driver.procedures['energy']['polaritonic-rks'] = run_polaritonic_scf
psi4.driver.procedures['energy']['polaritonic-rcis'] = run_polaritonic_scf
psi4.driver.procedures['energy']['polaritonic-tddft'] = run_polaritonic_scf
psi4.driver.procedures['energy']['polaritonic-uccsd'] = run_polaritonic_scf
psi4.driver.procedures['energy']['polaritonic-rpa'] = run_polaritonic_scf

