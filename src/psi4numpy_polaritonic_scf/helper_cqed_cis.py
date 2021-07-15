"""
Helper function for CQED_RHF

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

def cqed_cis(lam, molecule_string, psi4_options_dict, omega_val, include_dse=True):
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
        psi4_cis_energy : float
            Ground state energy from canonical RHF wavefunction from

        cqed_cis_energy : float
            Ground state energy of the CQED_RHF wavefunction

        cqed_cis_wavefunction : 1 x nocc * nvirt * nphoton array of floats
            Transformation vectors corresponding to the CQED_RHF orbitals

        Example
        -------
        >>> psi4_cis_energy, cqed_cis_energy, cqed_cis_vec = cqed_cis([0., 0., 1e-2], '''\nMg\nH 1 1.7\nsymmetry c1\n1 1\n''', 0.02)
        
    """
    # define geometry using the molecule_string
    mol = psi4.geometry(molecule_string)
    # define options for the calculation
    psi4.set_options(psi4_options_dict)
    # run psi4 to get ordinary scf energy and wavefunction object
    scf_e, wfn = psi4.energy('scf', return_wfn=True)
    
    # ==> Nuclear Repulsion Energy <==
    E_nuc = mol.nuclear_repulsion_energy()
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
    
    # grab all transformation vectors and store to a numpy array!
    C = np.asarray(wfn.Ca())
    
    # occupied orbitals:
    Co = wfn.Ca_subset("AO", "OCC")
    
    # virtual orbitals:
    Cv = wfn.Ca_subset("AO", "VIR")
    
    # grab all transformation vectors and store to a numpy array!
    C = np.asarray(wfn.Ca())
    
    # orbital energies
    eps     = np.asarray(wfn.epsilon_a())
    
    # ==> Nuclear Repulsion Energy <==
    E_nuc = mol.nuclear_repulsion_energy()
    
    print("\nNumber of occupied orbitals: %d" % ndocc)
    
    # 2 electron integrals in ao basis
    #I = np.asarray(mints.ao_eri())

    # 2 electron integrals in mo basis
    ovov = np.asarray(mints.mo_eri(Co, Cv, Co, Cv))
    
    # build the (oo|vv) integrals:
    oovv = np.asarray(mints.mo_eri(Co, Co, Cv, Cv))
    
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
        
    # We need to carry around the electric field dotted into the nuclear dipole moment
    # and the electric field dotted into the RHF electronic dipole expectation value...
    # so let's compute them here!
    
    # \lambda \cdot \mu_{nuc}
    l_dot_mu_nuc = lam[0] * mu_nuc_x + lam[1] * mu_nuc_y + lam[2] * mu_nuc_z
    # \lambda \cdot < \mu > where <\mu> contains electronic and nuclear contributions
    l_dot_mu_exp = lam[0] * mu_exp_x + lam[1] * mu_exp_y + lam[2] * mu_exp_z
    
    # dipole constants to add to E_RHF
    #  0.5 * (\lambda \cdot \mu_{nuc})** 2 
    #      - (\lambda \cdot <\mu> ) ( \lambda \cdot \mu_{nuc})
    # +0.5 * (\lambda \cdot <\mu>) ** 2
    d_c = 0.5 * l_dot_mu_nuc **2 - l_dot_mu_nuc * l_dot_mu_exp + 0.5 * l_dot_mu_exp ** 2
    
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
    
    # build the (ov|ov) integrals:
    ovov = np.asarray(mints.mo_eri(Co, Cv, Co, Cv))
    
    # build the (oo|vv) integrals:
    oovv = np.asarray(mints.mo_eri(Co, Co, Cv, Cv))
    # strip out occupied orbital energies, eps_o spans 0..ndocc-1
    eps_o = eps[:ndocc]
    
    # strip out virtual orbital energies, eps_v spans 0..nvirt-1
    eps_v = eps[ndocc:]

    # create Hamiltonian
    HCIS = np.zeros((2 + ndocc * nvirt * 2, 2 + ndocc * nvirt * 2))
    
    # (\lambda \cdot \mu_nuc - \lambda \cdot <\mu>) term
    dc_offset = l_dot_mu_nuc - l_dot_mu_exp
    
    # elements corresponding to <s|<\Phi_0 | H | \Phi_0>|t> go here
    # <0|\Phi_0| H |\Phi_0>0>
    
    # constant terms in the diagonals
    H_constants = scf_e
    H_dse = d_c
    H_blc = 0.0
    
    
    # add 1-electron contributions to diagonals
    for i in range(0,ndocc):
        # dipole terms scaled by dc_offset term
        H_dse += dc_offset * lam[0] * mu_cmo_x[i,i]       
        H_dse += dc_offset * lam[1] * mu_cmo_y[i,i]         
        H_dse += dc_offset * lam[2] * mu_cmo_z[i,i]
        
        # quadrupole terms
        H_dse += 0.5 * lam[0] * lam[0] * Q_cmo_xx[i,i]     
        H_dse += 0.5 * lam[1] * lam[1] * Q_cmo_yy[i,i]  
        H_dse += 0.5 * lam[2] * lam[2] * Q_cmo_zz[i,i] 
        H_dse += lam[0] * lam[1] * Q_cmo_xy[i,i]
        H_dse += lam[0] * lam[2] * Q_cmo_xz[i,i]
        H_dse += lam[1] * lam[2] * Q_cmo_yz[i,i]
        
    # add 2-electron contributions to diagonals
    for i in range(0,ndocc):
        for j in range(i+1, ndocc):
            # diagonal terms (xx, yy, zz)
            # xx
            H_dse += lam[0] * lam[0] * mu_cmo_x[i,i] * mu_cmo_x[j,j]
            H_dse -= 0.5 * lam[0] * lam[0] * mu_cmo_x[i,j] * mu_cmo_x[j,i]
            
            # yy
            H_dse += lam[1] * lam[1] * mu_cmo_y[i,i] * mu_cmo_y[j,j]
            H_dse -= 0.5 * lam[1] * lam[1] * mu_cmo_y[i,j] * mu_cmo_y[j,i]
            
            # zz
            H_dse += lam[2] * lam[2] * mu_cmo_z[i,i] * mu_cmo_z[j,j]
            H_dse -= 0.5 * lam[2] * lam[2] * mu_cmo_z[i,j] * mu_cmo_z[j,i]

            # off-diagonal terms (xy, xz, yz)
            # xy 
            H_dse += lam[0] * lam[1] * mu_cmo_x[i,i] * mu_cmo_y[j,j]
            H_dse -= 0.5 * lam[0] * lam[1] * mu_cmo_x[i,j] * mu_cmo_y[j,i]

            # yx
            H_dse += lam[1] * lam[0] * mu_cmo_y[i,i] * mu_cmo_x[j,j]
            H_dse -= 0.5 * lam[1] * lam[0] * mu_cmo_y[i,j] * mu_cmo_x[j,i]
            
            # xz
            H_dse += lam[0] * lam[2] * mu_cmo_x[i,i] * mu_cmo_z[j,j]
            H_dse -= 0.5 * lam[0] * lam[2] * mu_cmo_x[i,j] * mu_cmo_z[j,i]

            # zx
            H_dse += lam[2] * lam[0] * mu_cmo_z[i,i] * mu_cmo_x[j,j]
            H_dse -= 0.5 * lam[2] * lam[0] * mu_cmo_z[i,j] * mu_cmo_x[j,i]
            
            # yz
            H_dse += lam[1] * lam[2] * mu_cmo_y[i,i] * mu_cmo_z[j,j]
            H_dse -= 0.5 * lam[1] * lam[2] * mu_cmo_y[i,j] * mu_cmo_z[j,i]

            # zy
            H_dse += lam[2] * lam[1] * mu_cmo_z[i,i] * mu_cmo_y[j,j]
            H_dse -= 0.5 * lam[2] * lam[1] * mu_cmo_z[i,j] * mu_cmo_y[j,i]
        
    HCIS[0,0] = H_constants  
    HCIS[1,1] = H_constants + omega_val
    if include_dse:
        HCIS[0,0] += H_dse
        HCIS[1,1] += H_dse
    
    H_blc = l_dot_mu_exp
    # now sum over occupied orbitals
    for i in range (0, ndocc):
        H_blc -= lam[0] * mu_cmo_x[i,i]
        H_blc -= lam[1] * mu_cmo_y[i,i]
        H_blc -= lam[1] * mu_cmo_z[i,i]
        
    H_blc *= np.sqrt(omega_val / 2)
    
    ### off-diagonals for this block are the same!
    HCIS[0,1] = H_blc
    HCIS[1,0] = H_blc
    
    
    # elements corresponding to <s|<\Phi_i^a| H | \Phi_0|t> and <s|<\Phi_0| H | \Phi_i^a|t> go here!
    for i in range(0, ndocc):
        for a in range(0, nvirt):
            for s in range(0,2):
                # offset by 2 to account for the <s|<\Phi_0| H|\Phi_0>|t> block
                ias = 2*(i*nvirt + a) + s + 2
                
                
                for t in range(0,2):
                    H_dse = 0.
                    if s==t:
                        # quadrupole terms 
                        # xx
                        H_dse +=  0.5 * lam[0] * lam[0] * Q_cmo_xx[i,a]
                        # yy
                        H_dse += 0.5 * lam[1] * lam[1] * Q_cmo_yy[i,a]
                        # zz
                        H_dse += 0.5 * lam[2] * lam[2] * Q_cmo_zz[i,a]
                        # xy
                        H_dse += lam[0] * lam[1] * Q_cmo_xy[i,a]
                        # xz 
                        H_dse += lam[0] * lam[2] * Q_cmo_xz[i,a]
                        # yz 
                        H_dse += lam[1] * lam[2] * Q_cmo_yz[i,a]
                        
                        # 1e dipole terms scaled by dipole-offset 
                        # x
                        H_dse += dc_offset * lam[0] * mu_cmo_x[i,a]
                        # y
                        H_dse += dc_offset * lam[1] * mu_cmo_y[i,a]
                        # z
                        H_dse += dc_offset * lam[2] * mu_cmo_z[i,a]
                        
                        # 2e dipole terms
                        for j in range(0, ndocc):
                            # xx 
                            H_dse += lam[0] * lam[0] * mu_cmo_x[i,a] * mu_cmo_x[j,j]
                            H_dse -= 0.5 * lam[0] * lam[0] * mu_cmo_x[i,j] * mu_cmo_x[j,a]
                            
                            # yy
                            H_dse += lam[1] * lam[1] * mu_cmo_y[i,a] * mu_cmo_y[j,j]
                            H_dse -= 0.5 * lam[1] * lam[1] * mu_cmo_y[i,j] * mu_cmo_y[j,a]
                            
                            # zz
                            H_dse += lam[2] * lam[2] * mu_cmo_z[i,a] * mu_cmo_z[j,j]
                            H_dse -= 0.5 * lam[2] * lam[2] * mu_cmo_z[i,j] * mu_cmo_y[j,a]

                            # xy
                            H_dse += lam[0] * lam[1] * mu_cmo_x[i,a] * mu_cmo_y[j,j]
                            H_dse -= 0.5 * lam[0] * lam[1] * mu_cmo_x[i,j] * mu_cmo_y[j,a]
                            
                            # yx
                            H_dse += lam[1] * lam[0] * mu_cmo_y[i,a] * mu_cmo_x[j,j]
                            H_dse -= 0.5 * lam[1] * lam[0] * mu_cmo_y[i,j] * mu_cmo_x[j,a]
                            
                            # xz
                            H_dse += lam[0] * lam[2] * mu_cmo_x[i,a] * mu_cmo_z[j,j]
                            H_dse -= 0.5 * lam[0] * lam[0] * mu_cmo_z[i,j] * mu_cmo_z[j,a]
                            
                            # zx
                            H_dse += lam[2] * lam[0] * mu_cmo_z[i,a] * mu_cmo_x[j,j]
                            H_dse -= 0.5 * lam[2] * lam[0] * mu_cmo_z[i,j] * mu_cmo_x[j,a]
                            
                            # yz
                            H_dse += lam[1] * lam[2] * mu_cmo_y[i,a] * mu_cmo_z[j,j]
                            H_dse -= 0.5 * lam[1] * lam[2] * mu_cmo_y[i,j] * mu_cmo_z[j,a]
                            
                            # zy
                            H_dse += lam[2] * lam[1] * mu_cmo_z[i,a] * mu_cmo_y[j,j]
                            H_dse -= 0.5 * lam[2] * lam[1] * mu_cmo_z[i,j] * mu_cmo_y[j,a]

                    else:
                        H_blc = lam[0] * mu_cmo_x[i,a]
                        H_blc += lam[1] * mu_cmo_y[i,a]
                        H_blc += lam[2] * mu_cmo_z[i,a]
                        H_blc *= np.sqrt(omega_val/2)

                if include_dse:
                    HCIS[ias,t] = H_blc + H_dse
                    HCIS[t,ias] = H_blc + H_dse
                else:
                    HCIS[ias,t] = H_blc
                    HCIS[t,ias] = H_blc
          
                
    # elements corresponding to <s|<\Phi_i^a| H | \Phi_j^b|t>
    for i in range(0, ndocc):
        for a in range(0, nvirt):
            for s in range(0,2):
                ias = 2*(i*nvirt + a) + s + 2
                
                for j in range(0, ndocc):
                    for b in range(0, nvirt):
                        for t in range(0,2):
                            jbt = 2*(j*nvirt + b) + t + 2
                            H_dse = 0.
                            H_blc = 0.
                            H_elec = 0.

                            # most restrictive constraint
                            if s==t and i==j and a==b:
                                H_dse +=  d_c
                                H_elec += eps_v[a]
                                H_elec -= eps_o[i]

                            if s==t and i==j:
                                # quadrupole terms
                                # xx
                                H_dse += 0.5 * lam[0] * lam[0] * Q_cmo_xx[a,b]
                                # yy
                                H_dse += 0.5 * lam[1] * lam[1] * Q_cmo_yy[a,b]
                                # zz
                                H_dse += 0.5 * lam[2] * lam[2] * Q_cmo_zz[a,b]
                                # xy
                                H_dse += lam[0] * lam[1] * Q_cmo_xy[a,b]
                                # xz
                                H_dse += lam[0] * lam[2] * Q_cmo_xz[a,b]
                                # yz
                                H_dse += lam[1] * lam[2] * Q_cmo_xz[a,b]

                                # scaled dipole terms
                                # x
                                H_dse += dc_offset * lam[0] * mu_cmo_x[a,b]
                                # y
                                H_dse += dc_offset * lam[1] * mu_cmo_y[a,b]
                                # z 
                                H_dse += dc_offset * lam[2] * mu_cmo_z[a,b]

                            if s==t and a==b:
                                # quadrupole terms
                                # xx
                                H_dse -= 0.5 * lam[0] * lam[0] * Q_cmo_xx[i,j]
                                # yy
                                H_dse -= 0.5 * lam[1] * lam[1] * Q_cmo_yy[i,j]
                                # zz
                                H_dse -= 0.5 * lam[2] * lam[2] * Q_cmo_zz[i,j]
                                # xy
                                H_dse -= lam[0] * lam[1] * Q_cmo_xy[i,j]
                                # xz
                                H_dse -= lam[0] * lam[2] * Q_cmo_xz[i,j]
                                # yz
                                H_dse -= lam[1] * lam[2] * Q_cmo_xz[i,j]

                                # scaled dipole terms
                                # x
                                H_dse -= dc_offset * lam[0] * mu_cmo_x[i,j]
                                # y
                                H_dse -= dc_offset * lam[1] * mu_cmo_y[i,j]
                                # z 
                                H_dse -= dc_offset * lam[2] * mu_cmo_z[i,j]


                            if s==t:
                                # 2-e dipole terms
                                # xx
                                H_dse += lam[0] * lam[0] * mu_cmo_x[i,a] * mu_cmo_x[j,b]
                                H_dse -= 0.5 * lam[0] * lam[0] * mu_cmo_x[i,j] * mu_cmo_x[a,b]

                                # yy
                                H_dse += lam[1] * lam[1] * mu_cmo_y[i,a] * mu_cmo_y[j,b]
                                H_dse -= 0.5 * lam[1] * lam[1] * mu_cmo_y[i,j] * mu_cmo_y[a,b]
                                # zz
                                H_dse += lam[2] * lam[2] * mu_cmo_z[i,a] * mu_cmo_z[j,b]
                                H_dse -= 0.5 * lam[2] * lam[2] * mu_cmo_z[i,j] * mu_cmo_z[a,b]
                                
                                # xy
                                H_dse += lam[0] * lam[1] * mu_cmo_x[i,a] * mu_cmo_y[j,b]
                                H_dse -= 0.5 * lam[0] * lam[1] * mu_cmo_x[i,j] * mu_cmo_y[a,b]

                                # yx
                                H_dse += lam[1] * lam[0] * mu_cmo_y[i,a] * mu_cmo_x[j,b]
                                H_dse -= 0.5 * lam[1] * lam[0] * mu_cmo_y[i,j] * mu_cmo_x[a,b]

                                # xz
                                H_dse += lam[0] * lam[2] * mu_cmo_x[i,a] * mu_cmo_z[j,b]
                                H_dse -= 0.5 * lam[0] * lam[2] * mu_cmo_x[i,j] * mu_cmo_z[a,b]

                                # zx
                                H_dse += lam[2] * lam[0] * mu_cmo_z[i,a] * mu_cmo_x[j,b]
                                H_dse -= 0.5 * lam[2] * lam[0] * mu_cmo_z[i,j] * mu_cmo_x[a,b]

                                # yz
                                H_dse += lam[1] * lam[2] * mu_cmo_y[i,a] * mu_cmo_z[j,b]
                                H_dse -= 0.5 * lam[1] * lam[2] * mu_cmo_y[i,j] * mu_cmo_z[a,b]

                                # zy
                                H_dse += lam[2] * lam[1] * mu_cmo_z[i,a] * mu_cmo_y[j,b]
                                H_dse -= 0.5 * lam[2] * lam[1] * mu_cmo_z[i,j] * mu_cmo_y[a,b]

                                # 2e integral terms
                                H_elec += 2 * ovov[i, a, j, b] - oovv[i,j,a,b]

                            if (s==t+1 or s+1==t) and i==j and a==b:
                                # constant l dot mu term
                                H_blc += np.sqrt(omega_val/2) * l_dot_mu_exp
                            if (s==t+1 or s+1==t) and a==b:
                                # dipole coupling terms
                                # x
                                H_blc += np.sqrt(omega_val/2) * lam[0] * mu_cmo_x[i,j]
                                # y
                                H_blc += np.sqrt(omega_val/2) * lam[1] * mu_cmo_y[i,j]
                                # z
                                H_blc += np.sqrt(omega_val/2) * lam[2] * mu_cmo_z[i,j]
                            if (s==t+1 or s+1==t) and i==j:
                                # dipole coupling terms
                                # x
                                H_blc -= np.sqrt(omega_val/2) * lam[0] * mu_cmo_x[a,b]
                                # y
                                H_blc -= np.sqrt(omega_val/2) * lam[1] * mu_cmo_y[a,b]
                                # z
                                H_blc -= np.sqrt(omega_val/2) * lam[2] * mu_cmo_z[a,b]

                            if include_dse:
                                HCIS[ias,jbt] = H_elec + H_blc + H_dse
                            else:
                                HCIS[ias,jbt] = H_elec + H_blc
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