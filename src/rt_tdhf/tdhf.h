#ifndef TDHF_H
#define TDHF_H
#include<psi4/libmints/wavefunction.h>
#include<psi4/libmints/vector.h>
#include"fftw3.h"

#include<polaritonic_scf/hf.h>

using namespace psi;

namespace hilbert{

class TDHF: public PolaritonicHF {
public:
    TDHF(std::shared_ptr<psi::Wavefunction> reference_wavefunction,Options & options);
    ~TDHF();

    void common_init();
    double compute_energy();
    virtual bool same_a_b_orbs() const { return true; }
    virtual bool same_a_b_dens() const { return true; }

protected:
    double total_time, time_step, laser_amp, laser_freq, laser_time, total_iter;

    double * polarization_;

    std::shared_ptr<Matrix> T_;
    std::shared_ptr<Matrix> V_;

    std::shared_ptr<Matrix> eri;

    std::shared_ptr<Matrix> Dre;
    std::shared_ptr<Matrix> Dim;
    std::shared_ptr<Matrix> Fre;
    std::shared_ptr<Matrix> Fim;

    void build_fock(double * Dre, double * Dim, double curtime);
    double get_field(double t);

    // fourier transform
    void FFTW();
    void Spectrum();
    fftw_complex * td_dipole;
    fftw_complex * td_field;
    int fftw_iter;
    int pulse_shape_;

    // -i [in1, in2]
    // -i [in1_re + i in2_im, in2_re + i in2_im]
    void commutator(double * in1_re, double * in1_im, 
                    double * in2_re, double * in2_im, 
                    double * out_re, double * out_im, int dim);

};

}

#endif
