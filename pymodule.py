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

###############################################################
#
# begin andy's functions for quadrupole derivative integrals
#
###############################################################

from scipy.special import hyp1f1
from scipy.special import factorial2 as doublefactorial

np.set_printoptions(suppress=True, linewidth=200)

# These integrals are implemented using the McMurchie-Davidson scheme, described in Helgaker's
# purple book (Molecular Electronic Structure Theory).  The various equation numbers scattered
# throughout the code correspond to equation numbers in that book
#
# Andy Simmonett (01/22)

class GaussianShell(object):
    def __init__(self, origin, L, exps, coefs, atom_number):
        assert len(exps) == len(coefs)
        assert len(origin) == 3
        self.origin = np.array(origin)
        self.L = L
        self.exps = np.array(exps)
        self.coefs = np.array(coefs)
        self.atom_number = atom_number

def generate_quanta(am):
    """ Generates angular momentum compontents in CCA's lexicographic order http://dx.doi.org/10.1002/jcc.20815 """
    index = 0
    for l in range(am, -1, -1):
        for n in range(am - l + 1):
            m = am - l - n
            yield l, m, n, index
            index += 1

def cumulative_cart_dim(L):
    """ The number of Cartesian components in all shell with angular momentum L and lower """
    return ((L + 1) * (L + 2) * (L + 3)) // 6

def cart_dim(L):
    """ The number of Cartesian components in a shell with angular momentum L """
    return ((L + 1) * (L + 2)) // 2

def generate_E_matrix(maxam1, maxam2, P, A, B, a, b):
    """ Makes the Hermite->Cartesian conversion factors by recursion """
    Ex = np.zeros((maxam1 + 1, maxam2 + 1, maxam1 + maxam2 + 2))
    Ey = np.zeros((maxam1 + 1, maxam2 + 1, maxam1 + maxam2 + 2))
    Ez = np.zeros((maxam1 + 1, maxam2 + 1, maxam1 + maxam2 + 2))

    # 9.2.11
    p = a + b
    # 9.2.14
    AB = A - B
    PA = P - A
    PB = P - B
    # 9.2.12
    mu = a * b / p

    # 9.2.15
    Ex[0, 0, 0] = np.exp(-mu * AB[0] * AB[0])
    Ey[0, 0, 0] = np.exp(-mu * AB[1] * AB[1])
    Ez[0, 0, 0] = np.exp(-mu * AB[2] * AB[2])
    oo2p = 1 / (2 * p)
    for i in range(maxam1 + 1):
        for j in range(maxam2 + 1):
            uppert = i + j + 1
            # handle t = 0 case
            if i > 0:
                Ex[i, j, 0] += PA[0] * Ex[i - 1, j, 0] + Ex[i - 1, j, 1]
                Ey[i, j, 0] += PA[1] * Ey[i - 1, j, 0] + Ey[i - 1, j, 1]
                Ez[i, j, 0] += PA[2] * Ez[i - 1, j, 0] + Ez[i - 1, j, 1]
            elif j > 0:
                Ex[i, j, 0] += PB[0] * Ex[i, j - 1, 0] + Ex[i, j - 1, 1]
                Ey[i, j, 0] += PB[1] * Ey[i, j - 1, 0] + Ey[i, j - 1, 1]
                Ez[i, j, 0] += PB[2] * Ez[i, j - 1, 0] + Ez[i, j - 1, 1]
            for t in range(1, uppert):
                if i > 0:
                    # 9.5.6
                    Ex[i, j, t] += PA[0] * Ex[i - 1, j, t] + (t + 1) * Ex[i - 1, j, t + 1] + oo2p * Ex[i - 1, j, t - 1]
                    Ey[i, j, t] += PA[1] * Ey[i - 1, j, t] + (t + 1) * Ey[i - 1, j, t + 1] + oo2p * Ey[i - 1, j, t - 1]
                    Ez[i, j, t] += PA[2] * Ez[i - 1, j, t] + (t + 1) * Ez[i - 1, j, t + 1] + oo2p * Ez[i - 1, j, t - 1]
                elif j > 0:
                    # 9.5.7
                    Ex[i, j, t] += PB[0] * Ex[i, j - 1, t] + (t + 1) * Ex[i, j - 1, t + 1] + oo2p * Ex[i, j - 1, t - 1]
                    Ey[i, j, t] += PB[1] * Ey[i, j - 1, t] + (t + 1) * Ey[i, j - 1, t + 1] + oo2p * Ey[i, j - 1, t - 1]
                    Ez[i, j, t] += PB[2] * Ez[i, j - 1, t] + (t + 1) * Ez[i, j - 1, t + 1] + oo2p * Ez[i, j - 1, t - 1]
    return Ex, Ey, Ez

def generate_M_matrix(maxam, maxpow, PC, exp_a, exp_b):
    """ Generate multipole intermediates using equations
        (9.5.31) to (9.5.36) from Helgaker's book """
    p = exp_a + exp_b
    prefac = 1 / (2 * p)
    Mx = np.zeros((maxpow + 1, maxam + 3))
    My = np.zeros((maxpow + 1, maxam + 3))
    Mz = np.zeros((maxpow + 1, maxam + 3))

    Mx[0, 0] = My[0, 0] = Mz[0, 0] = np.sqrt(np.pi / p)
    for e in range(1, maxpow + 1):
        # t > e is zero!
        uppert = min(e + 1, maxam + 2)
        Mx[e, 0] += PC[0] * Mx[e - 1, 0] + prefac * Mx[e - 1, 1]
        My[e, 0] += PC[1] * My[e - 1, 0] + prefac * My[e - 1, 1]
        Mz[e, 0] += PC[2] * Mz[e - 1, 0] + prefac * Mz[e - 1, 1]
        for t in range(1, uppert):
            Mx[e, t] += PC[0] * Mx[e - 1, t] + prefac * Mx[e - 1, t + 1] + t * Mx[e - 1, t - 1]
            My[e, t] += PC[1] * My[e - 1, t] + prefac * My[e - 1, t + 1] + t * My[e - 1, t - 1]
            Mz[e, t] += PC[2] * Mz[e - 1, t] + prefac * Mz[e - 1, t + 1] + t * Mz[e - 1, t - 1]
    return Mx, My, Mz

def generate_M_pair_deriv(shell1, shell2, C, max_order):
    """ Computes derivatives of multipole moment integrals for a pair of shells """
    am1 = shell1.L
    am2 = shell2.L
    A = shell1.origin
    B = shell2.origin
    m_mats = [ np.zeros((cart_dim(am1), cart_dim(am2))) for i in range(6*cumulative_cart_dim(max_order)) ]
    for ca, a in zip(shell1.coefs, shell1.exps):
        for cb, b in zip(shell2.coefs, shell2.exps):
            p = a + b
            P = (a * A + b * B) / p
            PC = P - C
            prefac = ca * cb
            Ex, Ey, Ez = generate_E_matrix(am1 + 1, am2 + 1, P, A, B, a, b)
            Mx, My, Mz = generate_M_matrix(am1 + am2, max_order, PC, a, b)
            Sx = np.zeros((am1 + 2, am2 + 2, max_order + 1))
            Sy = np.zeros((am1 + 2, am2 + 2, max_order + 1))
            Sz = np.zeros((am1 + 2, am2 + 2, max_order + 1))
            for i in range(am1 + 2):
                for j in range(am2 + 2):
                    for e in range(max_order + 1):
                        for t in range(min(i + j, e) + 1):
                            # 9.5.39
                            Sx[i, j, e] += Ex[i, j, t] * Mx[e, t]
                            Sy[i, j, e] += Ey[i, j, t] * My[e, t]
                            Sz[i, j, e] += Ez[i, j, t] * Mz[e, t]

            m_count = 0
            for order in range(max_order + 1):
                for ex, ey, ez, _ in generate_quanta(order):
                    for l1, m1, n1, index1 in generate_quanta(am1):
                        for l2, m2, n2, index2 in generate_quanta(am2):
                            sx = Sx[l1, l2, ex]
                            sy = Sy[m1, m2, ey]
                            sz = Sz[n1, n2, ez]
                            # Use eqn. 9.3.30 to define derivatives
                            # Ax
                            DAx = -2 * a * Sx[l1+1, l2, ex]
                            if l1 > 0:
                                DAx += l1 * Sx[l1-1, l2, ex]
                            m_mats[m_count+0][index1, index2] += prefac * DAx * sy * sz
                            # Ay
                            DAy = -2 * a * Sy[m1+1, m2, ey]
                            if m1 > 0:
                                DAy += m1 * Sy[m1-1, m2, ey]
                            m_mats[m_count+1][index1, index2] += prefac * sx * DAy * sz
                            # Az
                            DAz = -2 * a * Sz[n1+1, n2, ez]
                            if n1 > 0:
                                DAz += n1 * Sz[n1-1, n2, ez]
                            m_mats[m_count+2][index1, index2] += prefac * sx * sy * DAz
                            # Bx
                            DBx = -2 * b * Sx[l1, l2+1, ex]
                            if l2 > 0:
                                DBx += l2 * Sx[l1, l2-1, ex]
                            m_mats[m_count+3][index1, index2] += prefac * DBx * sy * sz
                            # By
                            DBy = -2 * b * Sy[m1, m2+1, ey]
                            if m2 > 0:
                                DBy += m2 * Sy[m1, m2-1, ey]
                            m_mats[m_count+4][index1, index2] += prefac * sx * DBy * sz
                            # Bz
                            DBz = -2 * b * Sz[n1, n2+1, ez]
                            if n2 > 0:
                                DBz += n2 * Sz[n1, n2-1, ez]
                            m_mats[m_count+5][index1, index2] += prefac * sx * sy * DBz
                    m_count += 6
    return m_mats

def build_derivative_integral_matrices(shells, maxorder, D, natoms, callback, antisymmetric=False):
    nbf = np.sum([cart_dim(s.L) for s in shells])
    mats = [np.zeros((natoms, 3)) for i in range(cumulative_cart_dim(maxorder))]
    Poff = 0
    for P, shellP in enumerate(shells):
        Pdim = cart_dim(shellP.L)
        Patom = shellP.atom_number
        Qoff = 0
        for Q, shellQ in enumerate(shells[:P+1]):
            Qdim = cart_dim(shellQ.L)
            Qatom = shellQ.atom_number
            # make sure we don't double count diagonal blocks when (anti)symmetrizing
            scale = 1.0 if P == Q else 2.0
            Dblock = D[Poff : Poff + Pdim, Qoff : Qoff + Qdim]
            contributions = callback(shellP, shellQ)
            assert len(contributions) % 6 == 0
            nchunks = len(contributions) // 6
            for chunk in range(nchunks):
                # Ax
                mats[chunk][Patom, 0] += scale * np.tensordot(contributions[6*chunk+0], Dblock)
                # Ay
                mats[chunk][Patom, 1] += scale * np.tensordot(contributions[6*chunk+1], Dblock)
                # Az
                mats[chunk][Patom, 2] += scale * np.tensordot(contributions[6*chunk+2], Dblock)
                # Bx
                mats[chunk][Qatom, 0] += scale * np.tensordot(contributions[6*chunk+3], Dblock)
                # By
                mats[chunk][Qatom, 1] += scale * np.tensordot(contributions[6*chunk+4], Dblock)
                # Bz
                mats[chunk][Qatom, 2] += scale * np.tensordot(contributions[6*chunk+5], Dblock)
            Qoff += Qdim
        Poff += Pdim
    return mats

###############################################################
#
# end andy's functions for quadrupole derivative integrals
#
###############################################################


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
    elif ( lowername == 'polaritonic-rcis' ):
        psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_RCIS')
    elif ( lowername == 'polaritonic-uccsd' ):
        psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_UCCSD')

    # Compute a SCF reference, a wavefunction is return which holds the molecule used, orbitals
    # Fock matrices, and more
    #print('Attention! This SCF may be density-fitted.')
    ref_wfn = kwargs.get('ref_wfn', None)
    #if ref_wfn is None:
    #    ref_wfn = psi4.driver.scf_helper(name, **kwargs)
    if ref_wfn is None:
        if ( lowername == 'polaritonic-uks' ):
            func = psi4.core.get_option('HILBERT','CAVITY_QED_DFT_FUNCTIONAL')
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

def run_polaritonic_scf_gradient(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    polaritonic scf can be called via :py:func:`~driver.gradient`. For post-scf plugins.

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
    elif ( lowername == 'polaritonic-rcis' ):
        psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_RCIS')
    elif ( lowername == 'polaritonic-uccsd' ):
        psi4.core.set_local_option('HILBERT', 'HILBERT_METHOD', 'POLARITONIC_UCCSD')

    # Compute a SCF reference, a wavefunction is return which holds the molecule used, orbitals
    # Fock matrices, and more
    #print('Attention! This SCF may be density-fitted.')
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        if ( lowername == 'polaritonic-uks' ):
            func = psi4.core.get_option('HILBERT','CAVITY_QED_DFT_FUNCTIONAL')
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

    # gradient of photon-free hamiltonian

    # some quantities aren't set correctly in hilbert's wave functions, so we can't call
    # scfgrad directly. to get the photon-free part of the gradient, just push 
    # (i)   polaritonic-scf orbitals 
    # (ii)  polaritonic-scf orbital energies
    # (iii) polaritonic-scf densities
    # onto reference wave function 

    for irrep in range (0,ref_wfn.Ca().nirrep()):
        ref_wfn.Ca().nph[irrep][:,:] = rhf_wfn.Ca().nph[irrep][:,:]
        ref_wfn.Cb().nph[irrep][:,:] = rhf_wfn.Cb().nph[irrep][:,:]
        ref_wfn.Da().nph[irrep][:,:] = rhf_wfn.Da().nph[irrep][:,:]
        ref_wfn.Db().nph[irrep][:,:] = rhf_wfn.Db().nph[irrep][:,:]
        ref_wfn.epsilon_a().nph[irrep][:] = rhf_wfn.epsilon_a().nph[irrep][:]
        ref_wfn.epsilon_b().nph[irrep][:] = rhf_wfn.epsilon_b().nph[irrep][:]

    gradient = psi4.core.scfgrad(ref_wfn)

    # dipole self energy portion of gradient

    # OPDM
    Da = np.asarray(rhf_wfn.Da())
    Db = np.asarray(rhf_wfn.Db())

    # dipole integrals
    mints = psi4.core.MintsHelper(ref_wfn.basisset())
    dipole = mints.so_dipole()

    mu_z = np.asarray(dipole[2])

    # exchange contribution to dipole self energy 

    # D(p,q) = - mu(r,s) [ Da(p,r)Da(s,q) + Db(p,r) Da(s,q) ]

    tmpa = -np.einsum('rs,pr,sq->pq',mu_z, Da, Da) 
    tmpb = -np.einsum('rs,pr,sq->pq',mu_z, Db, Db)

    # test exchange energy from dressed RDM
    g = psi4.core.get_option("HILBERT","CAVITY_COUPLING_STRENGTH")
    w = psi4.core.get_option("HILBERT","CAVITY_FREQUENCY")
    lambda_z = g[2] * np.sqrt(2.0 * w[2])

    en  = 0.5 * lambda_z * lambda_z * np.einsum('pq,pq',tmpa,mu_z)
    en += 0.5 * lambda_z * lambda_z * np.einsum('pq,pq',tmpb,mu_z)
    #print('exchange dse:',en)

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
    for atom in range (0,natom):
        for cart in range (0,3):
            dse_gradient_z[atom,cart] = dse_gradient[atom*3+cart,2] 

    # scale by lambda^2
    dse_gradient_z_scaled = psi4.core.Matrix.from_array(dse_gradient_z)
    dse_gradient_z_scaled.scale(lambda_z*lambda_z)

    # quadrupole integral gradient

    # Make a lightweight list of shells
    basis = ref_wfn.basisset()
    shells = []
    for shellnum in range(basis.nshell()):
        shell = basis.shell(shellnum)
        xyz = mol.xyz(shell.ncenter)
        coefs = [shell.coef(primitive) for primitive in range(shell.nprimitive)]
        exps = [shell.exp(primitive) for primitive in range(shell.nprimitive)]
        shells.append(GaussianShell([xyz[0], xyz[1], xyz[2]], shell.am, exps, coefs, shell.ncenter))
    natoms = basis.molecule().natom()
    
    C = [0.0, 0.0, 0.0]
    maxorder = 2

    # now get densities in cartesian ao basis
    Da = np.asarray(rhf_wfn.Da_subset("CartAO"))
    Db = np.asarray(rhf_wfn.Db_subset("CartAO"))
    D = Da + Db
    derivmats = build_derivative_integral_matrices(shells, maxorder, D, natoms, lambda P, Q: generate_M_pair_deriv(P, Q, C, maxorder))
    #count = 0
    #for order in range(maxorder + 1):
    #    for ex, ey, ez, _ in generate_quanta(order):
    #        print(count,"X"*ex + "Y"*ey + "Z"*ez + " multipole derivatives")
    #        print(derivmats[count])
    #        print()
    #        count += 1    
    dse_gradient_z_scaled_2 = psi4.core.Matrix.from_array(derivmats[9])
    dse_gradient_z_scaled_2.scale(-0.5 * lambda_z*lambda_z)

    # print
    #dse_gradient_z_scaled.print_out()
    #dse_gradient_z_scaled_2.print_out()

    dse_gradient_z_scaled.add(dse_gradient_z_scaled_2)
    #dse_gradient_z_scaled.print_out()

    gradient.print_out()

    # total gradient
    gradient.add(dse_gradient_z_scaled)
    gradient.print_out()

    optstash.restore()

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
psi4.driver.procedures['energy']['polaritonic-rcis'] = run_polaritonic_scf
psi4.driver.procedures['energy']['polaritonic-uccsd'] = run_polaritonic_scf

psi4.driver.procedures['gradient']['polaritonic-uks'] = run_polaritonic_scf_gradient
psi4.driver.procedures['gradient']['polaritonic-uhf'] = run_polaritonic_scf_gradient

#psi4.driver.procedures['energy']['polaritonic-uhf'] = run_polaritonic_scf

