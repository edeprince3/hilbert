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

cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(hilbert CXX Fortran)

find_package(psi4 1.1 REQUIRED)

add_psi4_plugin(hilbert 

    plugin.cc

    doci/david.cc
    doci/doci_solver.cc
    doci/threeindexintegrals.cc
    doci/update_transformation_matrix.cc
    doci/write_tpdm.cc

    focas/focas_data.F90
    focas/focas_exponential.F90
    focas/focas_hessian.F90
    focas/focas_transform_driver.F90
    focas/focas_diis.F90
    focas/focas_genfock.F90
    focas/focas_interface.F90
    focas/focas_transform_oeints.F90
    focas/focas_driver.F90
    focas/focas_gradient.F90
    focas/focas_redundant.F90
    focas/focas_transform_teints.F90
    focas/focas_energy.F90
    focas/focas_gradient_hessian.F90
    focas/focas_semicanonical.F90
    focas/focas_full_hessian.F90
)
