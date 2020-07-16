/*
 *@BEGIN LICENSE
 *
 * v2RDM-CASSCF, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * Copyright (c) 2014, The Florida State University. All rights reserved.
 * 
 *@END LICENSE
 *
 */

#ifndef _python_api2_h_
#define _python_api2_h_

#include "v2rdm_solver.h"
#include "v2rdm_helper.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <psi4/libpsi4util/process.h>
#include <psi4/libmints/wavefunction.h>

using namespace psi;

namespace py = pybind11;
using namespace pybind11::literals;

namespace hilbert{

void export_v2RDMHelper(py::module& m) {
    py::class_<v2rdm_casscf::v2RDMHelper, std::shared_ptr<v2rdm_casscf::v2RDMHelper> >(m, "v2RDMHelper")
        .def(py::init<std::shared_ptr<Wavefunction>,Options &>())
        .def("set_orbitals", &v2RDMHelper::set_orbitals)
        .def("get_orbitals", &v2RDMHelper::get_orbitals)
        .def("get_opdm", &v2RDMHelper::get_opdm)
        .def("get_tpdm", &v2RDMHelper::get_tpdm)
        .def("get_opdm_sparse", &v2RDMHelper::get_opdm_sparse)
        .def("get_tpdm_sparse", &v2RDMHelper::get_tpdm_sparse)
        .def("compute_energy", &v2RDMHelper::compute_energy);

    py::class_<v2rdm_casscf::opdm, std::shared_ptr<v2rdm_casscf::opdm> >(m, "opdm")
        .def(py::init<>())
        .def_readwrite("i", &opdm::i)
        .def_readwrite("j", &opdm::j)
        .def_readwrite("value", &opdm::value);

    py::class_<v2rdm_casscf::tpdm, std::shared_ptr<v2rdm_casscf::tpdm> >(m, "tpdm")
        .def(py::init<>())
        .def_readwrite("i", &tpdm::i)
        .def_readwrite("j", &tpdm::j)
        .def_readwrite("k", &tpdm::k)
        .def_readwrite("l", &tpdm::l)
        .def_readwrite("value", &tpdm::value);
}

PYBIND11_MODULE(v2rdm_casscf, m) {
    m.doc() = "Python API of v2rdm_casscf: a variational 2-RDM-driven CASSCF plugin to Psi4";
    export_v2RDMHelper(m);
}

v2RDMHelper::v2RDMHelper(SharedWavefunction reference_wavefunction,Options & options)
{

    v2rdm = (std::shared_ptr<v2RDMSolver>)(new v2RDMSolver(reference_wavefunction,options));

}

v2RDMHelper::~v2RDMHelper()
{

}

double v2RDMHelper::compute_energy() {

    return v2rdm->compute_energy();

}

std::shared_ptr<Matrix> v2RDMHelper::get_opdm() {

    return v2rdm->get_opdm();

}

std::vector<opdm> v2RDMHelper::get_opdm_sparse(std::string type) {

    return v2rdm->get_opdm_sparse(type);

}

std::shared_ptr<Matrix> v2RDMHelper::get_tpdm() {

    return v2rdm->get_tpdm();

}

std::vector<tpdm> v2RDMHelper::get_tpdm_sparse(std::string type) {

    return v2rdm->get_tpdm_sparse(type);

}

void v2RDMHelper::orbital_locations(const std::string& orbitals, int* start, int* end) {

    v2rdm->orbital_locations(orbitals,start,end);

}

std::shared_ptr<Matrix> v2RDMHelper::get_orbitals(const std::string& orbital_name) {

    return v2rdm->get_orbitals(orbital_name);

}

void v2RDMHelper::set_orbitals(const std::string& orbital_name, SharedMatrix orbitals) {

    v2rdm->set_orbitals(orbital_name, orbitals);

}


} // End namespaces

#endif
