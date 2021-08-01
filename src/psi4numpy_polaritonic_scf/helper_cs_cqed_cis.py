"""
Helper function for CQED_CIS in the coherent state basis optimized by cqed_rhf

"""

__authors__   = ["Jon McTague", "Jonathan Foley"]
__credits__   = ["Jon McTague", "Jonathan Foley"]

__copyright_amp__ = "(c) 2014-2018, The Psi4NumPy Developers"
__license__   = "BSD-3-Clause"
__date__      = "2021-01-15"

# ==> Import Psi4, NumPy, & SciPy <==
import psi4
import numpy as np
import scipy.linalg as la
import time
from helper_cqed_rhf import *

def cs_cqed_cis(lam, molecule_string, psi4_options_dict, omega_val):
    """ Computes the QED-CIS energy and wavefunction

        Arguments
        ---------
        lambda_vector : 1 x 3 array of floats
            the electric field vector
        
        molecule_string : string
            specifies the molecular geometry

        omega_val : float
            the photon energy

        Returns
        -------
        psi4_scf_energy : float
            Ground state energy from canonical RHF wavefunction from

        cqed_cis_energy : float
            Array of excitation energies computed from cs_cqed_cis Hamiltonian

        cqed_cis_wavefunction : 1 x nocc * nvirt * nphoton array of floats
            cs_cqed_cis eigenvectors

        Example
        -------
        >>> psi4_cis_energy, cqed_cis_energy, cqed_cis_vec = cqed_cis([0., 0., 1e-2], '''\nMg\nH 1 1.7\nsymmetry c1\n1 1\n''', 0.02)
        
    """
    # define geometry using the molecule_string
    mol = psi4.geometry(molecule_string)
    # define options for the calculation
    psi4.set_options(psi4_options_dict)
    # run psi4 to get ordinary scf energy and wavefunction object
    #scf_e, wfn = psi4.energy('scf', return_wfn=True)

    # run cqed_rhf method
    cqed_rhf_dict = cqed_rhf(lam, molecule_string)
    
    # grab necessary quantities from cqed_rhf_dict
    scf_e  = cqed_rhf_dict['rhf_energy']
    cqed_scf_e = cqed_rhf_dict['cqed_rhf_energy']
    wfn = cqed_rhf_dict['psi4_wfn']
    C = cqed_rhf_dict['cqed_rhf_transformation_vectors']
    eps = cqed_rhf_dict['cqed_rhf_orbital_energies']
    cqed_rhf_dipole_moment = cqed_rhf_dict['cqed_rhf_dipole_moment']

    nmo = wfn.nmo()

    # Create instance of MintsHelper class
    mints = psi4.core.MintsHelper(wfn.basisset())
    
    # Grab data from wavfunction
    
    # number of doubly occupied orbitals
    ndocc   = wfn.nalpha()
    
    # total number of orbitals
    nmo     = wfn.nmo()
    
    # number of virtual orbitals
    nvirt   = nmo - ndocc

    # need to update the Co and Cv core matrix objects so we can
    # utlize psi4s fast integral transformation!

    # first collect rhf wfn object as dictionary
    wfn_dict = psi4.core.Wavefunction.to_file(wfn)
    
    # Get orbitals from CQED
    wfn_dict['matrix']['Ca'] = C
    wfn_dict['matrix']['Cb'] = C
    # update wfn object
    wfn = psi4.core.Wavefunction.from_file(wfn_dict) 

    # occupied orbitals as psi4 objects but they correspond to CQED-RHF orbitals
    Co = wfn.Ca_subset("AO", "OCC")
    
    # virtual orbitals same way
    Cv = wfn.Ca_subset("AO", "VIR")
    
    # 2 electron integrals in CQED-RHF basis
    ovov = np.asarray(mints.mo_eri(Co, Cv, Co, Cv))
    
    # build the (oo|vv) integrals:
    oovv = np.asarray(mints.mo_eri(Co, Co, Cv, Cv))

    # strip out occupied orbital energies, eps_o spans 0..ndocc-1
    eps_o = eps[:ndocc]
    
    # strip out virtual orbital energies, eps_v spans 0..nvirt-1
    eps_v = eps[ndocc:]
    
    # Extra terms for Pauli-Fierz Hamiltonian
    # nuclear dipole
    mu_nuc_x = mol.nuclear_dipole()[0]
    mu_nuc_y = mol.nuclear_dipole()[1]
    mu_nuc_z = mol.nuclear_dipole()[2]
    
    # dipole arrays in AO basis
    mu_ao_x = np.asarray(mints.ao_dipole()[0])
    mu_ao_y = np.asarray(mints.ao_dipole()[1])
    mu_ao_z = np.asarray(mints.ao_dipole()[2])
    
    # transform dipole array to canonical MO basis from ordinary RHF (no photon)
    mu_cmo_x = np.dot(C.T, mu_ao_x).dot(C)
    mu_cmo_y = np.dot(C.T, mu_ao_y).dot(C)
    mu_cmo_z = np.dot(C.T, mu_ao_z).dot(C)
    
    # compute components of electronic dipole moment <mu> from ordinary RHF (no photon)
    mu_exp_x = 0.0
    mu_exp_y = 0.0
    mu_exp_z = 0.0
    
    for i in range(0, ndocc):
        # double because this is only alpha terms!
        mu_exp_x += 2 * mu_cmo_x[i, i]
        mu_exp_y += 2 * mu_cmo_y[i, i]
        mu_exp_z += 2 * mu_cmo_z[i, i]
        
    # need to add the nuclear term to the expectation values above whic
    # only included the electronic term!
    mu_exp_x += mu_nuc_x
    mu_exp_y += mu_nuc_y
    mu_exp_z += mu_nuc_z
    
    # also have the dipole moment expectation value from the cqed_rhf_dict...
    # check to see if it matches!
    assert(np.isclose(mu_exp_x, cqed_rhf_dipole_moment[0]))
    assert(np.isclose(mu_exp_y, cqed_rhf_dipole_moment[1]))
    assert(np.isclose(mu_exp_z, cqed_rhf_dipole_moment[2]))

    # We need to carry around the electric field dotted into the nuclear dipole moment
    # and the electric field dotted into the RHF electronic dipole expectation value...
    # so let's compute them here!
    
    # \lambda \cdot \mu_{nuc}
    l_dot_mu_nuc = lam[0] * mu_nuc_x + lam[1] * mu_nuc_y + lam[2] * mu_nuc_z

    # \lambda \cdot < \mu > where <\mu> contains electronic and nuclear contributions
    l_dot_mu_exp = lam[0] * mu_exp_x + lam[1] * mu_exp_y + lam[2] * mu_exp_z

    # \lambda \cdot \mu_{el}
    l_dot_mu_el =  lam[0] * mu_cmo_x 
    l_dot_mu_el += lam[1] * mu_cmo_y
    l_dot_mu_el += lam[2] * mu_cmo_z
    
    # dipole constants to add to E_RHF
    #  0.5 * (\lambda \cdot \mu_{nuc})** 2 
    #      - (\lambda \cdot <\mu> ) ( \lambda \cdot \mu_{nuc})
    # +0.5 * (\lambda \cdot <\mu>) ** 2
    d_c = 0.5 * l_dot_mu_nuc **2 - l_dot_mu_nuc * l_dot_mu_exp + 0.5 * l_dot_mu_exp ** 2

    # again check to see if this thing is close to what we have from CQED-RHF calculation!
    assert np.isclose(d_c, cqed_rhf_dict['Nuclear Dipolar Energy'])
    
    # quadrupole arrays
    # Q_ao_xx[0,0] would correspond to the integral
    # of \phi_1 * e^2 * x * x * \phi_1 
    # where e is the charge of the electron, x is the
    # x-corrdinate, phi_1 is atomic orbital 1
    
    # Q_ao_xy[0,0] would correspond to the integral
    # of \phi_1 * e^2 * x * y * \phi_1
    Q_ao_xx = np.asarray(mints.ao_quadrupole()[0])
    Q_ao_xy = np.asarray(mints.ao_quadrupole()[1])
    Q_ao_xz = np.asarray(mints.ao_quadrupole()[2])
    Q_ao_yy = np.asarray(mints.ao_quadrupole()[3])
    Q_ao_yz = np.asarray(mints.ao_quadrupole()[4])
    Q_ao_zz = np.asarray(mints.ao_quadrupole()[5])
    
    # transform quadrupole array to canonical MO basis from ordinary RHF (no photon)
    Q_cmo_xx = np.dot(C.T, Q_ao_xx).dot(C)
    Q_cmo_xy = np.dot(C.T, Q_ao_xy).dot(C)
    Q_cmo_xz = np.dot(C.T, Q_ao_xz).dot(C)
    Q_cmo_yy = np.dot(C.T, Q_ao_yy).dot(C)
    Q_cmo_yz = np.dot(C.T, Q_ao_yz).dot(C)
    Q_cmo_zz = np.dot(C.T, Q_ao_zz).dot(C)

    # sum all terms of the quadrupole terms together
    Q_PF =  lam[0] * lam[0] * Q_cmo_xx
    Q_PF += lam[1] * lam[1] * Q_cmo_yy
    Q_PF += lam[2] * lam[2] * Q_cmo_zz
    Q_PF += 2 * lam[0] * lam[1] * Q_cmo_xy
    Q_PF += 2 * lam[0] * lam[2] * Q_cmo_xz
    Q_PF += 2 * lam[1] * lam[2] * Q_cmo_yz



    # create Hamiltonian for elements H[ias, jbt]
    HCIS = np.zeros((ndocc * nvirt * 2 + 2, ndocc * nvirt * 2 + 2))

    HCIS[0,0] = d_c
    HCIS[1,1] = np.sqrt(1) * omega_val + d_c

    
    # (\lambda \cdot \mu_nuc - \lambda \cdot <\mu>) term
    dc_offset = l_dot_mu_nuc - l_dot_mu_exp

    # elements corresponding to <s|<\Phi_0 | H | \Phi_i^a>|t>
    for s in range(0,2):
        for i in range(0,ndocc):
            for a in range(0,nvirt):
                A = a + ndocc
                for t in range(0,2):
                    iat = 2*(i*nvirt + a) + t + 2
                    HCIS[s,iat] = -np.sqrt(omega_val/2) *  l_dot_mu_el[i,A] * (s==t+1)
                    HCIS[s,iat] -= np.sqrt(omega_val/2) *  l_dot_mu_el[i,A] * (s+1==t)
                    HCIS[iat,s] = -np.sqrt(omega_val/2) *  l_dot_mu_el[i,A] * (s==t+1)
                    HCIS[iat,s] -= np.sqrt(omega_val/2) *  l_dot_mu_el[i,A] * (s+1==t)
    

    # elements corresponding to <s|<\Phi_i^a| H | \Phi_j^b|t>
    for i in range(0, ndocc):
        for a in range(0, nvirt):
            A = a+ndocc
            for s in range(0,2):
                ias = 2*(i*nvirt + a) + s + 2
                
                for j in range(0, ndocc):
                    for b in range(0, nvirt):
                        B = b+ndocc
                        for t in range(0,2):
                            jbt = 2*(j*nvirt + b) + t + 2
                            # ERIs
                            HCIS[ias,jbt] =  (2.0 * ovov[i, a, j, b] - oovv[i, j, a, b]) * (s==t)
                            # 2-electron dipole terms
                            # ordinary
                            HCIS[ias,jbt] += (2.0 * l_dot_mu_el[i,A] * l_dot_mu_el[j,B]) * (s==t)
                            # exchange
                            HCIS[ias,jbt] -= l_dot_mu_el[i,j] * l_dot_mu_el[A,B] * (s==t)
                            # 1-electron orbital energies
                            HCIS[ias,jbt] += eps_v[a] * (s==t) * (a==b) * (i==j)
                            HCIS[ias,jbt] -= eps_o[i] * (s==t) * (a==b) * (i==j)
                            # photonic and dipole energy term
                            HCIS[ias,jbt] += (omega_val * np.sqrt(t) + d_c) * (s==t) * (i==j) * (a==b)
                            # 1-electron dipole and quadrupole terms
                            HCIS[ias,jbt] += (dc_offset * l_dot_mu_el[A,B] - 0.5 * Q_PF[A,B]) * (s==t) * (i==j)
                            HCIS[ias,jbt] -= (dc_offset * l_dot_mu_el[i,j] - 0.5 * Q_PF[i,j]) * (s==t) * (a==b)
                            # bilinear coupling
                            HCIS[ias,jbt] += np.sqrt(omega_val/2) * l_dot_mu_exp * (i==j) * (a==b) * (s==t+1)
                            HCIS[ias,jbt] += np.sqrt(omega_val/2) * l_dot_mu_exp * (i==j) * (a==b) * (s+1==t)
                            HCIS[ias,jbt] += np.sqrt(omega_val/2) * l_dot_mu_el[i,j] * (a==b) * (s==t+1)
                            HCIS[ias,jbt] += np.sqrt(omega_val/2) * l_dot_mu_el[i,j] * (a==b) * (s+1==t)
                            HCIS[ias,jbt] -= np.sqrt(omega_val/2) * l_dot_mu_el[A,B] * (i==j) * (s==t+1)
                            HCIS[ias,jbt] -= np.sqrt(omega_val/2) * l_dot_mu_el[A,B] * (i==j) * (s+1==t)

    #print("now formed")                        
    #print(HCIS)
    # now diagonalize H
    # use eigh if Hermitian
    if np.isclose(np.imag(omega_val),0,1e-6):
        ECIS, CCIS = np.linalg.eigh(HCIS)
    # use eig if not-Hermitian... note will not
    # return both the left and right eigenvectos
    # but the equivalent scipy routine will...
    # need to do a bit more searching to figure out 
    # what we want to do with L and R vecs!
    else:
        ECIS, CCIS = np.linalg.eig(HCIS)

    # just return the energies for now!
    return scf_e, ECIS, CCIS
