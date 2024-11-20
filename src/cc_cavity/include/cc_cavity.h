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

#ifndef CC_CAVITY_H
#define CC_CAVITY_H


#include <tiledarray.h>
#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include "polaritonic_scf/hf.h"

#include "cc_cavity/misc/ta_helper.h"
#include "cc_cavity/misc/ta_diis.h"
#include "cc_cavity/misc/timer.h"


namespace hilbert {

    using std::string, std::unordered_map, std::unordered_set, std::vector, std::map,
          std::shared_ptr, std::to_string, std::make_shared;
    using namespace TA;
    using namespace psi;
    using namespace TA_Helper;

    // thread safe print function in scope of hilbert namespace
    static inline void Printf(const char *format, ...) {
        va_list argptr;
        va_start(argptr, format);
        char input[1024];
        vsprintf(input, format, argptr);
        va_end(argptr);
        TA::World& world_ = TA::get_default_world(); // TA world object
        world_.gop.serial_invoke(
                [=]() {
                    outfile->Printf("%s", input);
                }
        );
    }

    static inline std::string to_string(const TArrayD& array){
        std::stringstream ss;
        ss << array;
        return ss.str();
    }

    typedef map<std::string, TA::TArrayD> TArrayMap;

    class CC_Cavity : public PolaritonicHF {

    public:

        TA::World& world_ = TA::get_default_world(); // TA world object

        /**
         *  @brief CC_Cavity constructor
         *  @param[in] wfn      Wavefunction object
         *  @param[in] options  Options object
         */
        CC_Cavity(std::shared_ptr<Wavefunction> reference_wavefunction, Options & options, map<string, bool> &includes);
        ~CC_Cavity() override {
            if (epsilon_) free(epsilon_);
        }

        /**
         *  @brief initialize the CC_Cavity object
         */
        void initialize();

        /// basis set dimensions

        string cc_type_; // type of CC calculation

        size_t ns_, nQ_, // number of spin orbitals and number of basis orbitals
        o_, v_, oa_, ob_, va_, vb_; // number of occupied and virtual orbitals

        // initialize single, double, triple, and quadruple excitation dimensions
        size_t singleDim_ = 0, doubleDim_ = 0, tripleDim_ = 0, quadDim_ = 0;

        /// integrals
        TArrayMap C_blks_; // TiledArray MO transformation matrix blocks
        TArrayMap F_blks_; // Fock matrix in MO basis
        TArrayMap Dip_blks_; // dipole integrals in MO basis
        TArrayMap V_blks_; // antisymmetric 4-index integrals in MO basis

        TArrayD Qso_; // density-fitted integrals in AO basis
        TArrayD core_H_; // core Hamiltonian in AO basis
        TArrayD oe_cavity_terms_; // one-electron cavity integrals in AO basis


        /// operators
        std::map<string, string> idxs_; // map of indices for the operators
        TArrayMap amplitudes_; // amplitudes
        TArrayMap residuals_; // residuals
        TArrayMap tmps_; // temporary arrays
        std::map<string, double> scalars_; // scalar values

        map<string, double> resid_norms_; // residual norms
        TArrayMap Id_blks_; // identity matrix blocks
        
        // initialize included operator bools
        map<string, bool> &includes_;
        bool has_photon_ = includes_["t0_1"] || includes_["t1_1"] || includes_["t2_1"] || includes_["t3_1"] || includes_["t4_1"] ||
                           includes_["t0_2"] || includes_["t1_2"] || includes_["t2_2"] || includes_["t3_2"] || includes_["t4_2"] ||
                           options_.get_int("N_PHOTON_STATES") > 1;

        double lambda_[3] = {0,0,0}; // coupling strengths
        double * epsilon_; // orbital energies
        double cc_energy_ = 0.0; // CC energy


        /// timers
        Timer t_resid, t_ampUp, t_transform, t_oei, t_tei, t_ground;

        /**
         * @brief run calculation using the level of theory specified in the options
         */
        double compute_energy() override;

        /**
         * @brief get the index string for a TA object
         * @param array the TA object
         * @return the index string
         */
        static string get_index(const TArrayD& array){
            size_t rank = array.trange().rank();

            string idx = "";
            for (int i = 0; i < rank-1; i++) {
                idx += "x" + to_string(i) + ",";
            }
            idx += "x" + to_string(rank-1);
            return idx;
        }

        /**
         * @brief reset all tiles of a TA object
         * @param array the TA object
         */
        void zero_tiles(TArrayD& array) {
            array = TArrayD(world_, array.trange()); array.fill(0.0);
            world_.gop.fence();
        }

        // true if t1 amplitudes are folded into the integrals
        bool has_t1_integrals_ = false;

        /**
        * @brief transform integrals with MO Coefficients
        * @param use_t1 if true, fold the t1 amplitudes into MO Coefficients
        */
        virtual void transform_integrals(bool use_t1);


        /**
         * @brief set the density matrices from the similarity transformed density matrices
         * @param Da alpha density matrix in AO basis
         * @param Db beta density matrix in AO basis
         */
        void save_density(const SharedMatrix& Da, const SharedMatrix& Db){ Da_ = Da; Db_ = Db;}

        // DIIS object for amplitudes
        std::shared_ptr<DIISTA> diis_ta;

    protected:

        /**
         * @brief initialize the integrals and data for the calculation
         */
        virtual void init_integrals();

        /**
         * @brief print the dimensions of the calculation
         */
        virtual void print_dimensions();

        /**
         * @brief initializes the amplitudes
         */
        virtual void init_operators();

        /**
         * @brief calculate the CC energy and amplitudes
         * @return the CC energy
         */
        virtual double cc_iterations();

        /**
         * @brief use the modified MO Coefficients to build the the Fock, ERI, and Dipole integrals
         * @param C the MO Coefficients
         */
        virtual void apply_transform(TArrayMap &CL, TArrayMap &CR);

        /**
         * @brief build the OEI integrals in the MO basis.
         * @param C_new the MO coefficients to use
         */
        virtual void build_oei(TArrayMap &CL, TArrayMap &CR);

        /**
         * @brief build the TEI integrals in the MO basis.
         * @param C_new the MO coefficients to use
         */
        virtual void build_tei(TArrayMap &CL, TArrayMap &CR);

        /**
         * @brief build the antisymmetric 4-index integrals in the MO basis.
         * @param Qmo_blks the 3-index integrals in the MO basis
         */
        virtual void unpack_eris(TArrayMap &Qmo_blks);

        /**
         * @brief build the orbital energies
         */
        virtual void build_eps();

        /**
         * @brief update and extrapolate amplitudes using DIIS
         */
        virtual void update_amplitudes();

        // *** virtual functions to be implemented by derived classes ***

        /**
        * print header for CC iterations
         */
        virtual void print_iter_header() const = 0;

        /**
        * @brief print the CC energy and residual norm for current iteration
        */
        virtual void print_iteration(size_t iter, double energy, double dele, double tnorm) const = 0;

        /**
         * @brief compute the residual equations
         */
        virtual double build_residuals() = 0;

        /**
         * @brief update residual equations with the orbital energies
         */
        virtual void update_residuals() = 0;


        /**
         * @brief compute the norm of the residuals
         * @return the norm of the residuals
         */
        virtual double compute_residual_norms(bool return_tot = true) = 0;

        /**
        * @brief compute and print norms of the amplitudes
        */
        virtual void print_properties() = 0;

    };
}


#endif //CC_CAVITY_H
