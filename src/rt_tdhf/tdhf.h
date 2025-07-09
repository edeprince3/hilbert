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

    double * polarization;

    std::shared_ptr<Matrix> T;
    std::shared_ptr<Matrix> V;
    std::shared_ptr<Matrix> eri;
    std::shared_ptr<Matrix> Dre;
    std::shared_ptr<Matrix> Dim;
    std::shared_ptr<Matrix> Vext;
    std::shared_ptr<Matrix> Fre;
    std::shared_ptr<Matrix> Fim;
    std::vector<std::shared_ptr<Matrix> > dipole;

    void BuildFock(double * Dre, double * Dim, double curtime);
    double get_field(double t);

    // fourier transform
    void FFTW();
    void Spectrum();
    fftw_complex * td_dipole;
    fftw_complex * td_field;
    int fftw_iter;
    int pulse_shape_;

};

}

#endif
