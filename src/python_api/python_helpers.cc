#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "python_helpers.h"
#include "../v2rdm_doci/v2rdm_solver.h"
#include "../doci/doci_solver.h"
#include "../pp2rdm/pp2rdm_solver.h"

using namespace psi;

namespace py = pybind11;
using namespace pybind11::literals;

namespace psi{ 

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
    doci = (std::shared_ptr<doci::DOCISolver>)(new doci::DOCISolver(reference_wavefunction,options));
}

DOCIHelper::~DOCIHelper()
{
}

double DOCIHelper::compute_energy() {
    return doci->compute_energy();
}


pp2RDMHelper::pp2RDMHelper(SharedWavefunction reference_wavefunction,Options & options)
{
    pp2rdm = (std::shared_ptr<pp2rdm::pp2RDMSolver>)(new pp2rdm::pp2RDMSolver(reference_wavefunction,options));
}

pp2RDMHelper::~pp2RDMHelper()
{
}

double pp2RDMHelper::compute_energy() {
    return pp2rdm->compute_energy();
}

v2RDMHelper::v2RDMHelper(SharedWavefunction reference_wavefunction,Options & options)
{
    v2rdm = (std::shared_ptr<v2rdm_doci::v2RDMSolver>)(new v2rdm_doci::v2RDMSolver(reference_wavefunction,options));
}

v2RDMHelper::~v2RDMHelper()
{
}

double v2RDMHelper::compute_energy() {
    return v2rdm->compute_energy();
}

}

