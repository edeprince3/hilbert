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

def run_qed_scf(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    qed-scf can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('qed-scf')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    optstash = p4util.OptionsState(
        ['SCF', 'DF_INTS_IO'])

    psi4.core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    reference = psi4.core.get_global_option('REFERENCE').lower()

    if ( lowername == 'qed-scf'):
        if ( reference == 'rhf'):
            psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_RHF')
        elif ( reference == 'rohf'):
            psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_ROHF')
        else:
            psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_UHF')
    elif ( lowername == 'qed-dft' ):
        if ( reference == 'rks' or reference == 'rhf'):
            psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_RKS')
        else:
            psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_UKS')
    elif ( lowername == 'qed-cis' ):
        psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_RCIS')
    elif ( lowername == 'qed-ccsd' ):
        psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_UCCSD')
    elif ( lowername == 'qed-tddft' ):
        if ( reference == 'rks' or reference == 'rhf'):
            psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_RTDDFT')
        else:
            psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_UTDDFT')

    # Compute a SCF reference, a wavefunction is return which holds the molecule used, orbitals
    # Fock matrices, and more
    #print('Attention! This SCF may be density-fitted.')
    ref_wfn = kwargs.get('ref_wfn', None)
    #if ref_wfn is None:
    #    ref_wfn = psi4.driver.scf_helper(name, **kwargs)
    if ref_wfn is None:
        if ( lowername == 'qed-dft' or lowername == 'qed-tddft'):
            func = psi4.core.get_option('HILBERT','QED_DFT_FUNCTIONAL')
            en, ref_wfn = psi4.driver.energy(func, **kwargs, return_wfn=True)
        else :
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

def run_qed_scf_gradient(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    qed-scf can be called via :py:func:`~driver.gradient`. For post-scf plugins.

    >>> energy('qed-scf')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    optstash = p4util.OptionsState(
        ['SCF', 'DF_INTS_IO'])

    psi4.core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    reference = psi4.core.get_global_option('REFERENCE').lower()

    if ( lowername == 'qed-scf'):
        if ( reference == 'rhf'):
            psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_RHF')
        elif ( reference == 'rohf'):
            psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_ROHF')
        else:
            psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_UHF')
    elif ( lowername == 'qed-dft' ):
        if ( reference == 'rks' or reference == 'rhf'):
            psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_RKS')
        else:
            psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_UKS')
    elif ( lowername == 'qed-cis' ):
        psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_RCIS')
    elif ( lowername == 'qed-ccsd' ):
        psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_UCCSD')
    elif ( lowername == 'qed-tddft' ):
        if ( reference == 'rks' or reference == 'rhf'):
            psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_RTDDFT')
        else:
            psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_UTDDFT')

    # Compute a SCF reference, a wavefunction is return which holds the molecule used, orbitals
    # Fock matrices, and more
    #print('Attention! This SCF may be density-fitted.')
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        if ( lowername == 'qed-dft' or lowername == 'qed-tddft' ):

            # get functional from options
            func = psi4.core.get_option('HILBERT','QED_DFT_FUNCTIONAL')

            # check if dertype is present in kwargs and handle accordingly
            try:
                dertype = kwargs.pop('dertype') # must remove dertype for energy call if present
            except:
                 dertype = "gradient" # set default dertype for analytic gradient if not present

            # call energy and grab wfn
            en, ref_wfn = psi4.driver.energy(func, **kwargs, return_wfn=True)
            kwargs['dertype'] = dertype # restore dertype for gradient call
        else :
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

    # check if reference wave function is restricted
    if ("rks" in lowername or "rhf" in lowername or "rohf" in lowername):
        # copy alpha quantities to beta quantities in polaritonic wave function
        for irrep in range (0,ref_wfn.Cb().nirrep()):
            rhf_wfn.Cb().nph[irrep][:,:] = rhf_wfn.Ca().nph[irrep][:,:]
            rhf_wfn.Db().nph[irrep][:,:] = rhf_wfn.Da().nph[irrep][:,:]
            rhf_wfn.epsilon_b().nph[irrep][:] = rhf_wfn.epsilon_a().nph[irrep][:]

    # gradient of photon-free hamiltonian

    # some quantities aren't set correctly in hilbert's wave functions, so we can't call
    # scfgrad directly. to get the photon-free part of the gradient, just push 
    # (i)   polaritonic-scf orbitals 
    # (ii)  polaritonic-scf orbital energies
    # (iii) polaritonic-scf densities
    # onto reference wave function 

    # set alpha orbitals, densities, and energies
    for irrep in range (0,ref_wfn.Ca().nirrep()):
        ref_wfn.Ca().nph[irrep][:,:] = rhf_wfn.Ca().nph[irrep][:,:]
        ref_wfn.Cb().nph[irrep][:,:] = rhf_wfn.Cb().nph[irrep][:,:]
        ref_wfn.Da().nph[irrep][:,:] = rhf_wfn.Da().nph[irrep][:,:]
        ref_wfn.Db().nph[irrep][:,:] = rhf_wfn.Db().nph[irrep][:,:]
        ref_wfn.epsilon_a().nph[irrep][:] = rhf_wfn.epsilon_a().nph[irrep][:]
        ref_wfn.epsilon_b().nph[irrep][:] = rhf_wfn.epsilon_b().nph[irrep][:]

    #### call scfgrad for electron-only part of gradient ####
    gradient = psi4.core.scfgrad(ref_wfn)

    #### dipole self energy portion of gradient ####

    # OPDM
    Da = np.asarray(rhf_wfn.Da())
    Db = np.asarray(rhf_wfn.Db())

    # dipole integrals
    mints = psi4.core.MintsHelper(ref_wfn.basisset())
    dipole = mints.so_dipole()

    mu_z = np.asarray(dipole[2])
    if ( psi4.core.get_option("HILBERT","ROTATE_POLARIZATION_AXIS") == "YZX" ):
        mu_z = np.asarray(dipole[0])
    if ( psi4.core.get_option("HILBERT","ROTATE_POLARIZATION_AXIS") == "ZXY" ):
        mu_z = np.asarray(dipole[1])

    # exchange contribution to dipole self energy 

    #### D(p,q) = - mu(r,s) [ Da(p,r)Da(s,q) + Db(p,r) Da(s,q) ] ####

    tmpa = -np.einsum('rs,pr,sq->pq',mu_z, Da, Da) 
    tmpb = -np.einsum('rs,pr,sq->pq',mu_z, Db, Db)

    # test exchange energy from dressed RDM
    g = psi4.core.get_option("HILBERT","CAVITY_COUPLING_STRENGTH")
    w = psi4.core.get_option("HILBERT","CAVITY_FREQUENCY")
    lambda_z = g[2] * np.sqrt(2.0 * w[2])

    en  = 0.5 * lambda_z * lambda_z * np.einsum('pq,pq',tmpa,mu_z)
    en += 0.5 * lambda_z * lambda_z * np.einsum('pq,pq',tmpb,mu_z)

    D = tmpa + tmpb

    # symmetrize D because dipole_grad only uses 1/2 the elements
    D = 0.5 * ( D + np.einsum('rs->sr',D) )

    D = psi4.core.Matrix.from_array(D)

    # number of atoms
    mol = psi4.core.get_active_molecule()
    natom = mol.natom()

    tmp = mints.dipole_grad(D)
    dse_gradient = np.asarray(tmp)

    # unpack z-component 3N x 3 matrix (the third column)
    dse_gradient_z = np.zeros((natom,3))
    zdir = 2
    if ( psi4.core.get_option("HILBERT","ROTATE_POLARIZATION_AXIS") == "YZX" ):
        zdir = 0
    if ( psi4.core.get_option("HILBERT","ROTATE_POLARIZATION_AXIS") == "ZXY" ):
        zdir = 1
    for atom in range (0,natom):
        for cart in range (0,3):
            dse_gradient_z[atom,cart] = dse_gradient[atom*3+cart,zdir] 

    # scale by lambda^2
    dse_gradient_z_scaled = psi4.core.Matrix.from_array(dse_gradient_z)
    dse_gradient_z_scaled.scale(lambda_z*lambda_z)

    #### quadrupole integral gradient ####
      
    C = [0.0, 0.0, 0.0] # origin
    maxorder = 2 # quadrupole
    D = Da + Db # OPDM
    
    # symmetrize D because dipole_grad only uses 1/2 the elements
    D = 0.5 * ( D + np.einsum('rs->sr',D) )
    D = psi4.core.Matrix.from_array(D)
    
    # 3N x 9 matrix of quadrupole derivatives
    quad_grad = np.asarray(mints.multipole_grad(D, maxorder, C))
    
    # get requested component of quadrupole gradient
    zzdir = 8 # zz component
    if ( psi4.core.get_option("HILBERT","ROTATE_POLARIZATION_AXIS") == "YZX" ):
        zzdir = 3 # xx component
    if ( psi4.core.get_option("HILBERT","ROTATE_POLARIZATION_AXIS") == "ZXY" ):
        zzdir = 6 # yy component

    # unpack zz-component 3N x 3 matrix (the 9th column)
    dse_gradient_zz = np.zeros((natom,3))
    for atom in range (0,natom):
        for cart in range (0,3):
            dse_gradient_zz[atom,cart] = quad_grad[atom*3+cart,zzdir]
    
    dse_gradient_z_scaled_2 = psi4.core.Matrix.from_array(dse_gradient_zz)
    dse_gradient_z_scaled_2.scale(-0.5 * lambda_z*lambda_z)
    dse_gradient_z_scaled.add(dse_gradient_z_scaled_2)

    #### print out gradients ####

    # electronic gradient
    eg_norm = np.linalg.norm(gradient)
    eg_norm_xyz = np.linalg.norm(gradient, axis=0)
    psi4.core.print_out(f"\nElectronic Gradient: norm = {eg_norm:-20.12f}\n\n") # total norm
    gradient.print_out()
    psi4.core.Vector.from_array(eg_norm_xyz, name="Electronic Gradient |xyz|").print_out() # norm along each axis

    # total polaritonic gradient
    gradient.add(dse_gradient_z_scaled)
    pg_norm = np.linalg.norm(gradient)
    pg_norm_xyz = np.linalg.norm(gradient, axis=0)
    psi4.core.print_out(f"\nPolaritonic Gradient: norm = {pg_norm:-20.12f}\n\n") # total norm
    gradient.print_out() 
    psi4.core.Vector.from_array(pg_norm_xyz, name="Polaritonic Gradient |xyz|").print_out() # norm along each axis

    # difference between polaritonic and electronic gradients
    psi4.core.print_out(f"Gradient Difference: norm = {pg_norm-eg_norm:-20.12f}\n") # total norm
    pg_norm_xyz -= eg_norm_xyz
    psi4.core.Vector.from_array(pg_norm_xyz, name="Gradient Difference |xyz|").print_out() # norm along each axis

    optstash.restore()

    # set the gradient and return the wavefunction
    rhf_wfn.set_gradient(gradient)
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

def run_mcpdft(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    mcpdft can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('mcpdft')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    psi4.core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        raise ValidationError("""Error: mcpdft requires a reference wave function.""" )

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

    import hilbert
    rho_helper = hilbert.RealSpaceDensityHelper(new_wfn,options)

    # need Da and Db to build T+V+J
    Da = rho_helper.Da() 
    Db = rho_helper.Db() 

    # T + V
    mints = psi4.core.MintsHelper(new_wfn.basisset())

    V = np.asarray(mints.so_potential())
    T = np.asarray(mints.so_kinetic())

    Ta = psi4.core.Matrix.from_array(T)
    Tb = psi4.core.Matrix.from_array(T)

    Ta.transform(new_wfn.Ca())
    Tb.transform(new_wfn.Cb())

    kinetic_energy = Da.vector_dot(Ta)
    kinetic_energy += Db.vector_dot(Tb)

    Va = psi4.core.Matrix.from_array(V)
    Vb = psi4.core.Matrix.from_array(V)

    Va.transform(new_wfn.Ca())
    Vb.transform(new_wfn.Cb())

    en_potential_energy = Da.vector_dot(Va)
    en_potential_energy += Db.vector_dot(Vb)


    # classical coulomb energy
    jk = psi4.core.JK.build(new_wfn.get_basisset("ORBITAL"),
                           aux=new_wfn.get_basisset("DF_BASIS_SCF"))

    jk.set_do_K(False)
    jk.set_do_wK(False)
    jk.initialize()

    Ca = np.asarray(new_wfn.Ca())
    Cb = np.asarray(new_wfn.Cb())

    Cra = psi4.core.Matrix.from_array(Ca)
    Crb = psi4.core.Matrix.from_array(Cb)

    Cla = psi4.core.Matrix.from_array(Ca)
    Clb = psi4.core.Matrix.from_array(Cb)

    Cla.zero();
    Cla.gemm(False, True, 1.0, Cra, Da, 0.0);
    jk.C_left_add(Cla)
    jk.C_right_add(Cra)

    Clb.zero();
    Clb.gemm(False, True, 1.0, Crb, Db, 0.0);
    jk.C_left_add(Clb)
    jk.C_right_add(Crb)

    jk.compute()

    Ja = psi4.core.Matrix.from_array(np.asarray(jk.J()[0]))
    Jb = psi4.core.Matrix.from_array(np.asarray(jk.J()[1]))

    Ja.transform(new_wfn.Ca())
    Jb.transform(new_wfn.Cb())

    coulomb_energy = Da.vector_dot(Ja)
    coulomb_energy += Db.vector_dot(Jb)

    # xc contribution to the energy

    # density in real space
    rho_a = np.asarray(rho_helper.rho_a())
    rho_b = np.asarray(rho_helper.rho_b())
    rho = rho_a + rho_b

    # on-top pair density in real space
    pi = np.asarray(rho_helper.pi())

    # on-top ratio
    R = 4.0 * np.divide(pi, rho * rho)

    # translated densities
    translated_rho_a = np.zeros_like(rho)
    translated_rho_b = np.zeros_like(rho)

    # rhoa = [1 + zeta] * rho / 2
    # rhob = [1 - zeta] * rho / 2
    # zeta = sqrt(1-R), where 1-R > 0, 0 otherwise
    zeta = np.sqrt( 1.0 - R, out = np.zeros_like(R), where = 1.0 - R > 0 )

    translated_rho_a =  0.5 * rho * (1.0 + zeta)
    translated_rho_b =  0.5 * rho * (1.0 - zeta)

    # with translated rho_a, rho_b, evaluate xc contribution to the energy

    # LSDA:
    grid_w = rho_helper.grid_w()

    ex = 0.0
    ec = 0.0

    ex += np.sum(translated_rho_a**(4.0/3.0) * grid_w)
    ex += np.sum(translated_rho_b**(4.0/3.0) * grid_w)
    ex *= -2.0 ** (1.0 / 3.0 ) * 2.0 / 3.0 * 9.0 / 8.0 * (3.0 / np.pi)**(1.0 / 3.0)

    #na = np.sum(translated_rho_a * grid_w)
    #nb = np.sum(translated_rho_b * grid_w)

    #print('integrated number of electrons', na)
    #print('integrated number of electrons', nb)
    #print('kinetic energy', kinetic_energy)
    #print('electron-nucleus potential energy', en_potential_energy)
    #print('classical coulomb energy', coulomb_energy)
    #print('exchange energy', ex)
    #print('corrrelation energy', ec)
    #exit()

    optstash = p4util.OptionsState(
        ['SCF', 'DF_INTS_IO'])

    psi4.core.set_local_option('SCF', 'DF_INTS_IO', 'SAVE')

    v2rdm_wfn = kwargs.get('ref_wfn', None)
    if v2rdm_wfn is None:
        raise ValidationError("""Error: %s requires a reference wave function (v2rdm-casscf).""" % name)

    psi4.core.set_variable("V2RDM TOTAL ENERGY",v2rdm_wfn.energy())
   
    if ( (psi4.core.get_option('HILBERT', 'MCPDFT_METHOD') == '1DH_MCPDFT')
    or (psi4.core.get_option('HILBERT', 'MCPDFT_METHOD') == 'LS1DH_MCPDFT') ): 
        proc.run_dfmp2('mp2',**kwargs)

    if ('WBLYP' == psi4.core.get_option('HILBERT','MCPDFT_FUNCTIONAL')):
       func = 'BLYP'
    else:
       func = psi4.core.get_option('HILBERT','MCPDFT_FUNCTIONAL')
    ref_molecule = kwargs.get('molecule', psi4.core.get_active_molecule())
    base_wfn = psi4.core.Wavefunction.build(ref_molecule, psi4.core.get_global_option('BASIS'))
    ref_wfn = proc.scf_wavefunction_factory(func, base_wfn, psi4.core.get_global_option('REFERENCE'))

    # push v2rdm-casscf orbitals onto reference 
    for irrep in range (0,v2rdm_wfn.Ca().nirrep()): 
        ref_wfn.Ca().nph[irrep][:,:] = v2rdm_wfn.Ca().nph[irrep][:,:]
        ref_wfn.Cb().nph[irrep][:,:] = v2rdm_wfn.Cb().nph[irrep][:,:]

    # push v2rdm-casscf energies onto reference
    for irrep in range (0,v2rdm_wfn.epsilon_a().nirrep()): 
        ref_wfn.epsilon_a().nph[irrep][:] = v2rdm_wfn.epsilon_a().nph[irrep][:]
        ref_wfn.epsilon_b().nph[irrep][:] = v2rdm_wfn.epsilon_b().nph[irrep][:]

    # set the hilbert method
    psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'MCPDFT')

    # Call the Psi4 plugin
    # Please note that setting the reference wavefunction in this way is ONLY for plugins
    mcpdft_wfn = psi4.core.plugin('hilbert.so', ref_wfn)

    optstash.restore()

    return mcpdft_wfn


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

# qed-scf,dft,cc,tddft
psi4.driver.procedures['energy']['qed-scf'] = run_qed_scf
psi4.driver.procedures['energy']['qed-dft'] = run_qed_scf
psi4.driver.procedures['energy']['qed-tddft'] = run_qed_scf
psi4.driver.procedures['energy']['qed-ccsd'] = run_qed_scf

psi4.driver.procedures['gradient']['qed-scf'] = run_qed_scf_gradient
psi4.driver.procedures['gradient']['qed-dft'] = run_qed_scf_gradient

# mcpdft
psi4.driver.procedures['energy']['mcpdft'] = run_mcpdft

