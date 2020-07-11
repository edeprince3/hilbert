#
# @BEGIN LICENSE
#
# hilbert by DePrince Group, a plugin to:
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

    src/plugin.cc

    src/python_api/python_helpers.cc

    src/doci/david.cc
    src/doci/doci_solver.cc
    src/doci/write_tpdm.cc

    src/pp2rdm/diis.cc
    src/pp2rdm/minimization_routines.cc
    src/pp2rdm/pp2rdm_solver.cc
    src/pp2rdm/write_tpdm.cc

    src/misc/threeindexintegrals.cc

    src/focas/update_transformation_matrix.cc

    src/focas/focas_data.F90
    src/focas/focas_exponential.F90
    src/focas/focas_hessian.F90
    src/focas/focas_transform_driver.F90
    src/focas/focas_diis.F90
    src/focas/focas_genfock.F90
    src/focas/focas_interface.F90
    src/focas/focas_transform_oeints.F90
    src/focas/focas_driver.F90
    src/focas/focas_gradient.F90
    src/focas/focas_redundant.F90
    src/focas/focas_transform_teints.F90
    src/focas/focas_energy.F90
    src/focas/focas_gradient_hessian.F90
    src/focas/focas_semicanonical.F90
    src/focas/focas_full_hessian.F90

    src/v2rdm_doci/backtransform_tpdm.cc
    src/v2rdm_doci/basis.cc
    src/v2rdm_doci/cg_solver.cc
    src/v2rdm_doci/checkpoint.cc
    src/v2rdm_doci/d2.cc
    src/v2rdm_doci/d3.cc
    src/v2rdm_doci/g2.cc
    src/v2rdm_doci/gpu_transform_3index_teint.cc
    src/v2rdm_doci/integraltransform_sort_so_tpdm.cc
    src/v2rdm_doci/integraltransform_tpdm_unrestricted.cc
    src/v2rdm_doci/oei.cc
    src/v2rdm_doci/orbital_lagrangian.cc
    src/v2rdm_doci/q2.cc
    src/v2rdm_doci/sortintegrals.cc
    src/v2rdm_doci/t1.cc
    src/v2rdm_doci/t2.cc
    src/v2rdm_doci/tei.cc
    src/v2rdm_doci/update_primal.cc
    src/v2rdm_doci/v2rdm_solver.cc
    src/v2rdm_doci/write_tpdm.cc
    src/v2rdm_doci/write_tpdm_iwl.cc
)


install(DIRECTORY ${hilbert_BINARY_DIR}/
        DESTINATION ${CMAKE_INSTALL_PREFIX}
        USE_SOURCE_PERMISSIONS
        PATTERN "*.pyc" EXCLUDE
        PATTERN "*.mod" EXCLUDE
        PATTERN "CMakeCache.txt" EXCLUDE
        PATTERN "CMakeFiles" EXCLUDE
        PATTERN "cmake_install.cmake" EXCLUDE
        PATTERN "Makefile" EXCLUDE
        PATTERN "install_manifest.txt" EXCLUDE
)
