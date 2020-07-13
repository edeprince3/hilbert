/* 
 *  @BEGIN LICENSE
 * 
 *  Hilbert: a space for quantum chemistry plugins to Psi4 
 * 
 *  Copyright (c) 2020 by its authors (LICENSE).
 * 
 *  The copyrights for code used from other parties are included in
 *  the corresponding files.
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 * 
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 * 
 *  @END LICENSE
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "python_helpers.h"

#include <v2rdm_doci/v2rdm_solver.h>
#include <doci/doci_solver.h>
#include <pp2rdm/pp2rdm_solver.h>

using namespace psi;

namespace py = pybind11;
using namespace pybind11::literals;

namespace hilbert{

void export_HilbertHelper(py::module& m) {
    py::class_<DOCIHelper, std::shared_ptr<DOCIHelper> >(m, "DOCIHelper")
        .def(py::init<std::shared_ptr<Wavefunction>,Options &>())
        .def("compute_energy", &DOCIHelper::compute_energy);

    py::class_<pp2RDMHelper, std::shared_ptr<pp2RDMHelper> >(m, "pp2RDMHelper")
        .def(py::init<std::shared_ptr<Wavefunction>,Options &>())
        .def("compute_energy", &pp2RDMHelper::compute_energy);

    py::class_<v2RDMHelper, std::shared_ptr<v2RDMHelper> >(m, "v2RDMHelper")
        .def(py::init<std::shared_ptr<Wavefunction>,Options &>())
        .def("compute_energy", &v2RDMHelper::compute_energy);
}

PYBIND11_MODULE(hilbert, m) {
    m.doc() = "Python API of Hilbert";
    export_HilbertHelper(m);
}

DOCIHelper::DOCIHelper(SharedWavefunction reference_wavefunction,Options & options)
{
    doci = (std::shared_ptr<DOCISolver>)(new DOCISolver(reference_wavefunction,options));
}

DOCIHelper::~DOCIHelper()
{
}

double DOCIHelper::compute_energy() {
    return doci->compute_energy();
}


pp2RDMHelper::pp2RDMHelper(SharedWavefunction reference_wavefunction,Options & options)
{
    pp2rdm = (std::shared_ptr<pp2RDMSolver>)(new pp2RDMSolver(reference_wavefunction,options));
}

pp2RDMHelper::~pp2RDMHelper()
{
}

double pp2RDMHelper::compute_energy() {
    return pp2rdm->compute_energy();
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

}

