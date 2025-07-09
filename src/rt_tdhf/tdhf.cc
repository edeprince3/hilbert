#include<psi4/liboptions/liboptions.h>
#include<psi4/libmints/mintshelper.h>
#include<psi4/libmints/wavefunction.h>
#include<psi4/libmints/vector.h>
#include<psi4/libmints/matrix.h>
#include<psi4/libpsio/psio.h>
#include<psi4/physconst.h>
#include<psi4/libqt/qt.h>
#include<psi4/psifiles.h>

#include <misc/blas.h>

#ifdef _OPENMP
    #include<omp.h>
#endif

#include"tdhf.h"

using namespace psi;
using namespace fnocc;

namespace hilbert{

TDHF::TDHF(std::shared_ptr<Wavefunction> reference_wavefunction, Options& options_):
    PolaritonicHF(reference_wavefunction,options_) {
    common_init();
}

TDHF::~TDHF(){
}
void TDHF::common_init(){

    same_a_b_orbs_ = true;
    same_a_b_dens_ = true;

    // set the wavefunction name
    name_ = "TDHF";

    // check SCF type
    if ( options_.get_str("SCF_TYPE") != "PK") {
        throw PsiException("invalid SCF_TYPE for qed-tdhf",__FILE__,__LINE__);
    }

    std::shared_ptr<MintsHelper> mints (new MintsHelper(reference_wavefunction_));

    T_ = mints->so_kinetic();
    T_->transform(Ca_);

    V_ = mints->so_potential();
    V_->transform(Ca_);

    eri = mints->mo_eri(Ca_,Ca_);

    Dre = (std::shared_ptr<Matrix>) (new Matrix(nmo_, nmo_));
    Dim = (std::shared_ptr<Matrix>) (new Matrix(nmo_, nmo_));
    Fre = (std::shared_ptr<Matrix>) (new Matrix(nmo_, nmo_));
    Fim = (std::shared_ptr<Matrix>) (new Matrix(nmo_, nmo_));

    for (int i = 0; i < nalpha_; i++) {
        Dre->pointer()[i][i] = 1.0;
    }

    // get dipole integrals:
    dipole_ = mints->so_dipole();
    for (size_t i = 0; i < 3; i++) {
        dipole_[i]->transform(Ca_);
    }

    // get polarization:
    polarization_ = (double*)malloc(sizeof(double)*3);
    if (options_["POLARIZATION"].has_changed()){
       if (options_["POLARIZATION"].size() != 3)
          throw PsiException("The POLARIZATION array has the wrong dimensions",__FILE__,__LINE__);
       for (int i = 0; i < 3; i++) polarization_[i] = options_["POLARIZATION"][i].to_double();
    }else{
       polarization_[0] = 0.0;
       polarization_[1] = 0.0;
       polarization_[2] = 1.0;
    }

    total_time = options_.get_double("TOTAL_TIME");
    time_step  = options_.get_double("TIME_STEP");
    laser_amp  = options_.get_double("LASER_AMP");
    laser_freq = options_.get_double("LASER_FREQ");
    laser_time = options_.get_double("LASER_TIME");
    total_iter = total_time / time_step + 1;

    // which pulse shape do we want?
    if (options_.get_str("LASER_SHAPE") == "SIN_SQUARED") {
        // from prl:
        pulse_shape_ = 0;
    }else if (options_.get_str("LASER_SHAPE") == "TRAPEZOID") {
        // from 2007 schlegel paper (jcp 126, 244110 (2007))
        pulse_shape_ = 1;
    }else if (options_.get_str("LASER_SHAPE") == "CONTINUOUS") {
        // continuous wave for rabi flopping
        pulse_shape_ = 2;
    }

    // correlation function or dipole acceleration (fourier transformed)
    td_dipole = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(int)(total_iter));
    td_field = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(int)(total_iter));
}

void TDHF::build_fock(double * Dre, double * Dim, double curtime) {

    Fre->zero();
    Fim->zero();
    Fre->add(V_);
    Fre->add(T_);

    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            double dumre = 0.0;
            double dumim = 0.0;
            for (int a = 0; a < nmo_; a++) {
                for (int b = 0; b < nmo_; b++) {
                    dumre += (2.0 * eri->pointer()[i*nmo_+j][a*nmo_+b] - eri->pointer()[i*nmo_+a][j*nmo_+b]) * (Dre[a*nmo_+b]);
                    dumim += (2.0 * eri->pointer()[i*nmo_+j][a*nmo_+b] - eri->pointer()[i*nmo_+a][j*nmo_+b]) * (Dim[a*nmo_+b]);
                }
            }
            Fre->pointer()[i][j] += dumre;
            Fim->pointer()[i][j] += dumim;
        }
    }

    double sigma = laser_time*0.5;

    // add external field

    double field = get_field(curtime);

    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            Fre->pointer()[i][j] -= field * dipole_[0]->pointer()[i][j] * polarization_[0];
            Fre->pointer()[i][j] -= field * dipole_[1]->pointer()[i][j] * polarization_[1];
            Fre->pointer()[i][j] -= field * dipole_[2]->pointer()[i][j] * polarization_[2];
        }
    }
}

double TDHF::get_field(double t) {

    double field = 0.0;
    if (pulse_shape_ == 0 ) {
    
        // from prl:
        if ( t < laser_time ) {
            field = sin(M_PI*t/(laser_time));
            field *= field*laser_amp*sin(laser_freq*t);
        }
    
    } else if ( pulse_shape_ == 1 ) {
    
        // from 2007 schlegel paper (jcp 126, 244110 (2007))
        if (t <= 2.0 * M_PI / laser_freq)      field = laser_freq * t / (2.0 * M_PI) * laser_amp;
        else if (t <= 4.0 * M_PI / laser_freq) field = laser_amp;
        else if (t <= 6.0 * M_PI / laser_freq) field = (3.0 - laser_freq * t / (2.0 * M_PI) ) * laser_amp;
        field *= sin(laser_freq*t);
    
    } else if ( pulse_shape_ == 2 ) {
    
        // continuous wave for rabi flopping
        field = laser_amp*sin(laser_freq*t);
    
    }

    return field;
}

// -i [in1, in2]
// -i [in1_re + i in2_im, in2_re + i in2_im]
void TDHF::commutator(double * in1_re, double * in1_im, 
                      double * in2_re, double * in2_im, 
                      double * out_re, double * out_im, int dim) {

        // real:  in1_re.in2_im + in1_im.in2_re - in2_im.in1_re - in2_re.in1_im
        F_DGEMM('n', 'n', dim, dim, dim,  1.0, in1_re, dim, in2_im, dim, 0.0, out_re, dim);
        F_DGEMM('n', 'n', dim, dim, dim,  1.0, in1_im, dim, in2_re, dim, 1.0, out_re, dim);
        F_DGEMM('n', 'n', dim, dim, dim, -1.0, in2_im, dim, in1_re, dim, 1.0, out_re, dim);
        F_DGEMM('n', 'n', dim, dim, dim, -1.0, in2_re, dim, in1_im, dim, 1.0, out_re, dim);

        // imag: -in1_re.in2_re + in1_im.in1_im + in2_re.in1_re - in2_im.in1_im
        F_DGEMM('n', 'n', dim, dim, dim, -1.0, in1_re, dim, in2_re, dim, 0.0, out_im, dim);
        F_DGEMM('n', 'n', dim, dim, dim,  1.0, in1_im, dim, in2_im, dim, 1.0, out_im, dim);
        F_DGEMM('n', 'n', dim, dim, dim,  1.0, in2_re, dim, in1_re, dim, 1.0, out_im, dim);
        F_DGEMM('n', 'n', dim, dim, dim, -1.0, in2_im, dim, in1_im, dim, 1.0, out_im, dim); 
} 

double TDHF::compute_energy() {

    double * tempre = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * tempim = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * k1re = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * k2re = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * k3re = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * k4re = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * k1im = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * k2im = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * k3im = (double*)malloc(nmo_*nmo_*sizeof(double));
    double * k4im = (double*)malloc(nmo_*nmo_*sizeof(double));
    memset((void*)tempre,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)tempim,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)k1re,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)k2re,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)k3re,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)k4re,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)k1im,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)k2im,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)k3im,'\0',nmo_*nmo_*sizeof(double));
    memset((void*)k4im,'\0',nmo_*nmo_*sizeof(double));

    // i dD / dt = [F,D]
    // real:  dDre / dt = Fre.Dim + Fim.Dre - Dim.Fre - Dre.Fim
    // imag:  dDim / dt = Fre.Dre - Fim.Dim - Dre.Fre + Dim.Fim
    fftw_iter   = 0;

    for (size_t iter = 0; iter < total_iter; iter++) {

        // RK4
        // y(n+1) = y( n ) + 1/6 h ( k1 + 2k2 + 2k3 + k4 )
        // t(n+1) = t( n ) + h

        // k1     = f( t( n )         , y( n ) )
        build_fock(&(Dre->pointer())[0][0],&(Dim->pointer())[0][0],iter * time_step);

        commutator(Dre->pointer()[0], Dim->pointer()[0], Fre->pointer()[0], Fim->pointer()[0], k1re, k1im, nmo_);

        C_DCOPY(nmo_*nmo_, Dre->pointer()[0], 1, tempre, 1);
        C_DCOPY(nmo_*nmo_, Dim->pointer()[0], 1, tempim, 1);
        C_DAXPY(nmo_*nmo_, 0.5 * time_step, k1re, 1, tempre, 1);
        C_DAXPY(nmo_*nmo_, 0.5 * time_step, k1im, 1, tempim, 1);

        // k2     = f( t( n + 1/2 h ) , y( n ) + h/2 k1 )
        build_fock(tempre,tempim,(iter+0.5)*time_step);

        commutator(tempre, tempim, Fre->pointer()[0], Fim->pointer()[0], k2re, k2im, nmo_);

        C_DCOPY(nmo_*nmo_, Dre->pointer()[0], 1, tempre, 1);
        C_DCOPY(nmo_*nmo_, Dim->pointer()[0], 1, tempim, 1);
        C_DAXPY(nmo_*nmo_, 0.5 * time_step, k2re, 1, tempre, 1);
        C_DAXPY(nmo_*nmo_, 0.5 * time_step, k2im, 1, tempim, 1);

        // k3     = f( t( n + 1/2 h ) , y( n ) + h/2 k2 )
        build_fock(tempre,tempim,(iter+0.5)*time_step);

        commutator(tempre, tempim, Fre->pointer()[0], Fim->pointer()[0], k3re, k3im, nmo_);

        C_DCOPY(nmo_*nmo_, Dre->pointer()[0], 1, tempre, 1);
        C_DCOPY(nmo_*nmo_, Dim->pointer()[0], 1, tempim, 1);
        C_DAXPY(nmo_*nmo_, time_step, k3re, 1, tempre, 1);
        C_DAXPY(nmo_*nmo_, time_step, k3im, 1, tempim, 1);

        // k4     = f( t( n + h )     , y( n ) + h k3 )
        build_fock(tempre,tempim,(iter+1)*time_step);

        commutator(tempre, tempim, Fre->pointer()[0], Fim->pointer()[0], k4re, k4im, nmo_);

        // y(n+1) = y( n ) + 1/6 h ( k1 + 2k2 + 2k3 + k4 )

        C_DAXPY(nmo_*nmo_, 1.0 / 6.0 * time_step, k1re, 1, Dre->pointer()[0], 1);
        C_DAXPY(nmo_*nmo_, 2.0 / 6.0 * time_step, k2re, 1, Dre->pointer()[0], 1);
        C_DAXPY(nmo_*nmo_, 2.0 / 6.0 * time_step, k3re, 1, Dre->pointer()[0], 1);
        C_DAXPY(nmo_*nmo_, 1.0 / 6.0 * time_step, k4re, 1, Dre->pointer()[0], 1);

        C_DAXPY(nmo_*nmo_, 1.0 / 6.0 * time_step, k1im, 1, Dim->pointer()[0], 1);
        C_DAXPY(nmo_*nmo_, 2.0 / 6.0 * time_step, k2im, 1, Dim->pointer()[0], 1);
        C_DAXPY(nmo_*nmo_, 2.0 / 6.0 * time_step, k3im, 1, Dim->pointer()[0], 1);
        C_DAXPY(nmo_*nmo_, 1.0 / 6.0 * time_step, k4im, 1, Dim->pointer()[0], 1);

        // evaluate dipole moment:
        double dpre = C_DDOT(nmo_*nmo_,&(Dre->pointer())[0][0],1,&(dipole_[0]->pointer())[0][0],1);
        dpre       += C_DDOT(nmo_*nmo_,&(Dre->pointer())[0][0],1,&(dipole_[1]->pointer())[0][0],1);
        dpre       += C_DDOT(nmo_*nmo_,&(Dre->pointer())[0][0],1,&(dipole_[2]->pointer())[0][0],1);
        double dpim = C_DDOT(nmo_*nmo_,&(Dim->pointer())[0][0],1,&(dipole_[0]->pointer())[0][0],1);
        dpim       += C_DDOT(nmo_*nmo_,&(Dim->pointer())[0][0],1,&(dipole_[1]->pointer())[0][0],1);
        dpim       += C_DDOT(nmo_*nmo_,&(Dim->pointer())[0][0],1,&(dipole_[2]->pointer())[0][0],1);

        // time-dependent dipole
        td_dipole[fftw_iter][0] = dpre;
        td_dipole[fftw_iter][1] = dpim;

        // time-dependent field
        double field = get_field(iter * time_step);
        td_field[fftw_iter][0] = field;
        td_field[fftw_iter][1] = 0.0;

        fftw_iter++;

        double en = C_DDOT(nmo_*nmo_,&Fre->pointer()[0][0],1,&Dre->pointer()[0][0],1)
                  + C_DDOT(nmo_*nmo_,&Fim->pointer()[0][0],1,&Dim->pointer()[0][0],1)
                  + C_DDOT(nmo_*nmo_,&Dre->pointer()[0][0],1,&T_->pointer()[0][0],1)
                  + C_DDOT(nmo_*nmo_,&Dre->pointer()[0][0],1,&V_->pointer()[0][0],1);
        en += enuc_;

        outfile->Printf("@TDHF TIME %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf\n",iter*time_step, dpre, dpim, field, en);

    }

    free(tempre);
    free(tempim);
    free(k1re);
    free(k2re);
    free(k3re);
    free(k4re);
    free(k1im);
    free(k2im);
    free(k3im);
    free(k4im);

    // fourier transform and spectrum
    FFTW();
    Spectrum();

    free(td_dipole);
    free(td_field);

    return 0.0;
}

void TDHF::FFTW(){
    fftw_plan p_dipole;
    p_dipole = fftw_plan_dft_1d((int)(fftw_iter),td_dipole,td_dipole,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(p_dipole);
    fftw_destroy_plan(p_dipole);

    fftw_plan p_field;
    p_field = fftw_plan_dft_1d((int)(fftw_iter),td_field,td_field,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(p_field);
    fftw_destroy_plan(p_field);
}

// output absorption spectrum
void TDHF::Spectrum(){
  outfile->Printf("\n");
  outfile->Printf("        ***********************************************************\n");
  outfile->Printf("        *                                                         *\n");
  outfile->Printf("        *             f(w) = 2/3 * w * FT( mu(t) )                *\n");
  outfile->Printf("        *                                                         *\n");
  outfile->Printf("        ***********************************************************\n");
  outfile->Printf("\n");
  outfile->Printf("                                w(eV)");
  outfile->Printf("                 f(w)\n");

  // maximum frequency to output (eV)
  double max_freq = 100.0;

  for (size_t i = 1; i < (size_t)fftw_iter; i++){
      double w = 2.0 * M_PI * i / ( fftw_iter * time_step );
      if (w * pc_hartree2ev > max_freq) break;

      double dipr = td_dipole[i][0];
      double dipi = td_dipole[i][1];
      double er = td_field[i][0];
      double ei = td_field[i][1];

      double valr  = (dipr * er + dipi * ei) / (er * er + ei * ei);
      double vali  = (dipi * er - dipr * ei) / (er * er + ei * ei);
      double val  = sqrt(valr * valr + vali * vali);

      outfile->Printf("      @Frequency %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf\n",w * pc_hartree2ev, 2.0 / 3.0 * w * valr, 2.0 / 3.0 * w * vali, 2.0 / 3.0 * w * val, td_dipole[i][0], td_dipole[i][1], td_field[i][0], td_field[i][1]);
  }
}

}
