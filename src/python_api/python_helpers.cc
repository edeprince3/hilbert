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

#include <v2rdm_doci/v2rdm_doci_solver.h>
#include <v2rdm_casscf/v2rdm_solver.h>
#include <doci/doci_solver.h>
#include <pp2rdm/pp2rdm_solver.h>
#include <p2rdm/p2rdm_solver.h>
#include <misc/real_space_density.h>

using namespace psi;

namespace py = pybind11;
using namespace pybind11::literals;

namespace hilbert{

void export_HilbertHelper(py::module& m) {

    // real space density
    py::class_<RealSpaceDensityHelper, std::shared_ptr<RealSpaceDensityHelper> >(m, "RealSpaceDensityHelper")
        .def(py::init<std::shared_ptr<Wavefunction>,Options &>())
        .def("grid_x", &RealSpaceDensityHelper::grid_x)
        .def("grid_y", &RealSpaceDensityHelper::grid_y)
        .def("grid_z", &RealSpaceDensityHelper::grid_z)
        .def("grid_w", &RealSpaceDensityHelper::grid_w)
        .def("pi", &RealSpaceDensityHelper::pi)
        .def("build_rho_from_disk", &RealSpaceDensityHelper::build_rho_from_disk)
        .def("rho", &RealSpaceDensityHelper::rho)
        .def("rho_a", &RealSpaceDensityHelper::rho_a)
        .def("rho_b", &RealSpaceDensityHelper::rho_b)
        .def("rho_a_x", &RealSpaceDensityHelper::rho_a_x)
        .def("rho_a_y", &RealSpaceDensityHelper::rho_a_y)
        .def("rho_a_z", &RealSpaceDensityHelper::rho_a_z)
        .def("rho_b_x", &RealSpaceDensityHelper::rho_b_x)
        .def("rho_b_y", &RealSpaceDensityHelper::rho_b_y)
        .def("rho_b_z", &RealSpaceDensityHelper::rho_b_z)
        .def("xc_hole", &RealSpaceDensityHelper::xc_hole)
        .def("Da", &RealSpaceDensityHelper::Da)
        .def("Db", &RealSpaceDensityHelper::Db);

    // doci
    py::class_<DOCIHelper, std::shared_ptr<DOCIHelper> >(m, "DOCIHelper")
        .def(py::init<std::shared_ptr<Wavefunction>,Options &>())
        .def("compute_energy", &DOCIHelper::compute_energy);

    // pp2rdm
    py::class_<pp2RDMHelper, std::shared_ptr<pp2RDMHelper> >(m, "pp2RDMHelper")
        .def(py::init<std::shared_ptr<Wavefunction>,Options &>())
        .def("compute_energy", &pp2RDMHelper::compute_energy);

    // p2rdm
    py::class_<p2RDMHelper, std::shared_ptr<p2RDMHelper> >(m, "p2RDMHelper")
        .def(py::init<std::shared_ptr<Wavefunction>,Options &>())
        .def("compute_energy", &p2RDMHelper::compute_energy);

    // v2rdm-doci
    py::class_<v2RDM_DOCIHelper, std::shared_ptr<v2RDM_DOCIHelper> >(m, "v2RDM_DOCIHelper")
        .def(py::init<std::shared_ptr<Wavefunction>,Options &>())
        .def("compute_energy", &v2RDM_DOCIHelper::compute_energy);

    // v2rdm-casscf
    py::class_<v2RDMHelper, std::shared_ptr<v2RDMHelper> >(m, "v2RDMHelper")
        .def(py::init<std::shared_ptr<Wavefunction>,Options &>())
        .def(py::init<int,int,int,std::vector<double>,std::vector<double>,Options &>())
        .def("set_orbitals", &v2RDMHelper::set_orbitals)
        .def("get_orbitals", &v2RDMHelper::get_orbitals)
        .def("get_opdm", &v2RDMHelper::get_opdm)
        .def("get_tpdm", &v2RDMHelper::get_tpdm)
        .def("get_opdm_sparse", &v2RDMHelper::get_opdm_sparse)
        .def("get_tpdm_sparse", &v2RDMHelper::get_tpdm_sparse)
        .def("compute_energy", &v2RDMHelper::compute_energy);

    py::class_<opdm, std::shared_ptr<opdm> >(m, "opdm")
        .def(py::init<>())
        .def_readwrite("i", &opdm::i)
        .def_readwrite("j", &opdm::j)
        .def_readwrite("value", &opdm::value);

    py::class_<tpdm, std::shared_ptr<tpdm> >(m, "tpdm")
        .def(py::init<>())
        .def_readwrite("i", &tpdm::i)
        .def_readwrite("j", &tpdm::j)
        .def_readwrite("k", &tpdm::k)
        .def_readwrite("l", &tpdm::l)
        .def_readwrite("value", &tpdm::value);
}

PYBIND11_MODULE(hilbert, m) {
    m.doc() = "Python API of Hilbert";
    export_HilbertHelper(m);
}

// begin real space density

RealSpaceDensityHelper::RealSpaceDensityHelper(SharedWavefunction reference_wavefunction,Options & options)
{
    real_space_density = (std::shared_ptr<RealSpaceDensity>)(new RealSpaceDensity(reference_wavefunction,options));
}

RealSpaceDensityHelper::~RealSpaceDensityHelper()
{
}

std::vector<double> RealSpaceDensityHelper::grid_x() {
    std::shared_ptr<Vector> vec = real_space_density->grid_x();
    double * vec_p = vec->pointer();
    std::vector<double> return_val(vec_p,vec_p+vec->dim(0));
    return return_val;
}
std::vector<double> RealSpaceDensityHelper::grid_y() {
    std::shared_ptr<Vector> vec = real_space_density->grid_y();
    double * vec_p = vec->pointer();
    std::vector<double> return_val(vec_p,vec_p+vec->dim(0));
    return return_val;
}
std::vector<double> RealSpaceDensityHelper::grid_z() {
    std::shared_ptr<Vector> vec = real_space_density->grid_z();
    double * vec_p = vec->pointer();
    std::vector<double> return_val(vec_p,vec_p+vec->dim(0));
    return return_val;
}
std::vector<double> RealSpaceDensityHelper::grid_w() {
    std::shared_ptr<Vector> vec = real_space_density->grid_w();
    double * vec_p = vec->pointer();
    std::vector<double> return_val(vec_p,vec_p+vec->dim(0));
    return return_val;
}
void RealSpaceDensityHelper::build_rho_from_disk() {
    real_space_density->BuildRhoFromDisk();
}
std::vector<double> RealSpaceDensityHelper::rho() {
    std::shared_ptr<Vector> vec = real_space_density->pi();
    double * vec_p = vec->pointer();
    std::vector<double> return_val(vec_p,vec_p+vec->dim(0));
    return return_val;
}
std::vector<double> RealSpaceDensityHelper::rho_a() {
    std::shared_ptr<Vector> vec = real_space_density->rho_a();
    double * vec_p = vec->pointer();
    std::vector<double> return_val(vec_p,vec_p+vec->dim(0));
    return return_val;
}
std::vector<double> RealSpaceDensityHelper::rho_b() {
    std::shared_ptr<Vector> vec = real_space_density->rho_b();
    double * vec_p = vec->pointer();
    std::vector<double> return_val(vec_p,vec_p+vec->dim(0));
    return return_val;
}
std::vector<double> RealSpaceDensityHelper::rho_a_x() {
    std::shared_ptr<Vector> vec = real_space_density->rho_a_x();
    double * vec_p = vec->pointer();
    std::vector<double> return_val(vec_p,vec_p+vec->dim(0));
    return return_val;
}
std::vector<double> RealSpaceDensityHelper::rho_a_y() {
    std::shared_ptr<Vector> vec = real_space_density->rho_a_y();
    double * vec_p = vec->pointer();
    std::vector<double> return_val(vec_p,vec_p+vec->dim(0));
    return return_val;
}
std::vector<double> RealSpaceDensityHelper::rho_a_z() {
    std::shared_ptr<Vector> vec = real_space_density->rho_a_z();
    double * vec_p = vec->pointer();
    std::vector<double> return_val(vec_p,vec_p+vec->dim(0));
    return return_val;
}
std::vector<double> RealSpaceDensityHelper::rho_b_x() {
    std::shared_ptr<Vector> vec = real_space_density->rho_b_x();
    double * vec_p = vec->pointer();
    std::vector<double> return_val(vec_p,vec_p+vec->dim(0));
    return return_val;
}
std::vector<double> RealSpaceDensityHelper::rho_b_y() {
    std::shared_ptr<Vector> vec = real_space_density->rho_b_y();
    double * vec_p = vec->pointer();
    std::vector<double> return_val(vec_p,vec_p+vec->dim(0));
    return return_val;
}
std::vector<double> RealSpaceDensityHelper::rho_b_z() {
    std::shared_ptr<Vector> vec = real_space_density->rho_b_z();
    double * vec_p = vec->pointer();
    std::vector<double> return_val(vec_p,vec_p+vec->dim(0));
    return return_val;
}
std::shared_ptr<Matrix> RealSpaceDensityHelper::Da() {
    return real_space_density->Da();
}
std::shared_ptr<Matrix> RealSpaceDensityHelper::Db() {
    return real_space_density->Db();
}
std::vector<double> RealSpaceDensityHelper::xc_hole(double x, double y, double z) {
    std::shared_ptr<Vector> vec = real_space_density->xc_hole(x,y,z);
    double * vec_p = vec->pointer();
    std::vector<double> return_val(vec_p,vec_p+vec->dim(0));
    return return_val;
}
std::vector<double> RealSpaceDensityHelper::pi() {
    std::shared_ptr<Vector> vec = real_space_density->pi();
    double * vec_p = vec->pointer();
    std::vector<double> return_val(vec_p,vec_p+vec->dim(0));
    return return_val;
}


// begin DOCI

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

// begin pair p2RDM

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

// begin p2RDM

p2RDMHelper::p2RDMHelper(SharedWavefunction reference_wavefunction,Options & options)
{
    p2rdm = (std::shared_ptr<p2RDMSolver>)(new p2RDMSolver(reference_wavefunction,options));
}

p2RDMHelper::~p2RDMHelper()
{
}

double p2RDMHelper::compute_energy() {
    return p2rdm->compute_energy();
}

// begin v2RDM-DOCI

v2RDM_DOCIHelper::v2RDM_DOCIHelper(SharedWavefunction reference_wavefunction,Options & options)
{
    v2rdm_doci = (std::shared_ptr<v2RDM_DOCISolver>)(new v2RDM_DOCISolver(reference_wavefunction,options));
}

v2RDM_DOCIHelper::~v2RDM_DOCIHelper()
{
}

double v2RDM_DOCIHelper::compute_energy() {
    return v2rdm_doci->compute_energy();
}

// begin v2RDM

// default constructor for molecular or hubbard hamiltonian
v2RDMHelper::v2RDMHelper(SharedWavefunction reference_wavefunction,Options & options){
    v2rdm = (std::shared_ptr<v2RDMSolver>)(new v2RDMSolver(reference_wavefunction,options));
}

// constructor for externally-defined hamiltonian
v2RDMHelper::v2RDMHelper(int nalpha, int nbeta, int nmo, std::vector<double> h, std::vector<double> g, Options & options) {
    v2rdm = (std::shared_ptr<v2RDMSolver>)(new v2RDMSolver(nalpha,nbeta,nmo,h,g,options));
}

v2RDMHelper::~v2RDMHelper(){
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


}

