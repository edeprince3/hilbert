from __future__ import print_function

"""
A reference implementation of cavity quantum electrodynamics 
configuration interactions singles.
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

# Set Psi4 & NumPy Memory Options
psi4.set_memory('2 GB')
psi4.core.set_output_file('output.dat', False)

numpy_memory = 2

# basis set etc
psi4.set_options({'basis':        'cc-pVDZ',
                  'scf_type':     'pk',
                  'reference':    'rhf',
                  'mp2_type':     'conv',
                  'save_jk': True,
                  'e_convergence': 1e-8,
                  'd_convergence': 1e-8})

# MgH+ string
mol_string = """Mg
H 1 2.0
symmetry c1
1 1
"""


# electric field
Ex = 1e-3
Ey = 1e-3
Ez = 5e-2
lam = np.array([Ex, Ey, Ez])


# run cqed_rhf and collect the results!
# don't update the dipole moment at each scf step:
rhf_e, no_ud_cqed_rhf_e, cqed_vecs = cqed_rhf(lam, mol_string, self_consistent_dipole=False)
# do update the dipole moment at each scf step:
rhf_e, cqed_rhf_e, cqed_vecs = cqed_rhf(lam, mol_string, self_consistent_dipole=True)

print(f'Updating the dipole in the SCF steps changed the energy by {cqed_rhf_e-no_ud_cqed_rhf_e} Hartrees')

