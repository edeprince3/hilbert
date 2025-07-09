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

    T = mints->so_kinetic();
    T->transform(Ca_);

    V = mints->so_potential();
    V->transform(Ca_);

    eri = mints->mo_eri(Ca_,Ca_);

    Dre        = (std::shared_ptr<Matrix>) (new Matrix(nmo_, nmo_));
    Dim        = (std::shared_ptr<Matrix>) (new Matrix(nmo_, nmo_));
    Vext       = (std::shared_ptr<Matrix>) (new Matrix(nmo_, nmo_));
    Fre        = (std::shared_ptr<Matrix>) (new Matrix(nmo_, nmo_));
    Fim        = (std::shared_ptr<Matrix>) (new Matrix(nmo_, nmo_));

    for (int i = 0; i < nalpha_; i++) {
        Dre->pointer()[i][i] = 1.0;
    }

    // get dipole integrals:
    dipole = mints->so_dipole();
    for (size_t i = 0; i < 3; i++) {
        dipole[i]->transform(Ca_);
    }

    // get polarization:
    polarization = (double*)malloc(sizeof(double)*3);
    if (options_["POLARIZATION"].has_changed()){
       if (options_["POLARIZATION"].size() != 3)
          throw PsiException("The POLARIZATION array has the wrong dimensions",__FILE__,__LINE__);
       for (int i = 0; i < 3; i++) polarization[i] = options_["POLARIZATION"][i].to_double();
    }else{
       polarization[0] = 0.0;
       polarization[1] = 0.0;
       polarization[2] = 1.0;
    }

    BuildFock(&(Dre->pointer())[0][0],&(Dim->pointer())[0][0],0.0);

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

    // pad the correlation function with zeros just to get more output points

    // correlation function or dipole acceleration (fourier transformed)
    td_dipole = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(int)(total_iter));
    td_field = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(int)(total_iter));
}

void TDHF::BuildFock(double * Dre, double * Dim, double curtime) {

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
            dumre += T->pointer()[i][j] + V->pointer()[i][j];
            Fre->pointer()[i][j] = dumre;
            Fim->pointer()[i][j] = dumim;
        }
    }

/*
    for (int q = 0; q < nQ; q++) {
        double dumre = 0.0;
        double dumim = 0.0;
        for (int a = 0; a < nmo_; a++) {
            for (int b = 0; b < nmo_; b++) {
                dumre += 2.0 * Qmo[q*nmo_*nmo_+a*nmo_+b] * Dre[a*nmo_+b];
                dumim += 2.0 * Qmo[q*nmo_*nmo_+a*nmo_+b] * Dim[a*nmo_+b];
            }
        }
        Ire[q] = dumre;
        Iim[q] = dumim;
    }
    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            double dumre = 0.0;
            double dumim = 0.0;
            for (int q = 0; q < nQ; q++) {
                dumre += Qmo[q*nmo_*nmo_+i*nmo_+j] * Ire[q];
                dumim += Qmo[q*nmo_*nmo_+i*nmo_+j] * Iim[q];
            }
            dumre += T->pointer()[i][j] + V->pointer()[i][j];
            Fre->pointer()[i][j] = dumre;
            Fim->pointer()[i][j] = dumim;
        }
    }

    // K(re)
    F_DGEMM('t','n',nmo_,nmo_*nQ,nmo_,1.0,Dre,nmo_,Qmo,nmo_,0.0,Ire,nmo_);
    for (int q = 0; q < nQ; q++) {
        for (int a = 0; a < nmo_; a++) {
            for (int j = 0; j < nmo_; j++) {
                Iim[a*nmo_+j] = Ire[q*nmo_*nmo_+j*nmo_+a];
            }
        }
        C_DCOPY(nmo_*nmo_,Iim,1,Ire+q*nmo_*nmo_,1);
    }
    F_DGEMM('n','t',nmo_,nmo_,nQ*nmo_,-1.0,Ire,nmo_,Qmo,nmo_,1.0,&(Fre->pointer()[0][0]),nmo_);
    // K(im)
    F_DGEMM('t','n',nmo_,nmo_*nQ,nmo_,1.0,Dim,nmo_,Qmo,nmo_,0.0,Ire,nmo_);
    for (int q = 0; q < nQ; q++) {
        for (int a = 0; a < nmo_; a++) {
            for (int j = 0; j < nmo_; j++) {
                Iim[a*nmo_+j] = Ire[q*nmo_*nmo_+j*nmo_+a];
            }
        }
        C_DCOPY(nmo_*nmo_,Iim,1,Ire+q*nmo_*nmo_,1);
    }
    F_DGEMM('n','t',nmo_,nmo_,nQ*nmo_,-1.0,Ire,nmo_,Qmo,nmo_,1.0,&(Fim->pointer()[0][0]),nmo_);*/

    Vext->zero();
    double sigma = laser_time*0.5;

    // add external field

    double field = get_field(curtime);

    for (int i = 0; i < nmo_; i++) {
        for (int j = 0; j < nmo_; j++) {
            Fre->pointer()[i][j] -= field * dipole[0]->pointer()[i][j] * polarization[0];
            Fre->pointer()[i][j] -= field * dipole[1]->pointer()[i][j] * polarization[1];
            Fre->pointer()[i][j] -= field * dipole[2]->pointer()[i][j] * polarization[2];
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

    double * factorB = (double*)malloc(sizeof(double)*5);
    double * factorb = (double*)malloc(sizeof(double)*5);

    // factors for symplectic integrator.  they come from that sanz-serna paper.
    double fac,tntf = 1./3924.;
    factorB[0] = (642.+sqrt(471.))*tntf;
    factorB[1] = 121.*(12.-sqrt(471.))*tntf;
    factorB[2] = 1. - 2.*(factorB[0]+factorB[1]);
    factorB[3] = factorB[1];
    factorB[4] = factorB[0];

    factorb[0] = 6./11.;
    factorb[1] = .5 - factorb[0];
    factorb[2] = factorb[1];
    factorb[3] = factorb[0];
    factorb[4] = 0.;


    for ( int iter = 0; iter < total_iter; iter++ ) {

        // RK4
        // y(n+1) = y( n ) + 1/6 h ( k1 + 2k2 + 2k3 + k4 )
        // t(n+1) = t( n ) + h

        // k1     = f( t( n )         , y( n ) )
        BuildFock(&(Dre->pointer())[0][0],&(Dim->pointer())[0][0],iter * time_step);

        F_DGEMM('n','n',nmo_,nmo_,nmo_,1.0,&(Dim->pointer()[0][0]),nmo_,&(Fre->pointer()[0][0]),nmo_,0.0,k1re,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,1.0,&(Dre->pointer()[0][0]),nmo_,&(Fim->pointer()[0][0]),nmo_,1.0,k1re,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,-1.0,&(Fim->pointer()[0][0]),nmo_,&(Dre->pointer()[0][0]),nmo_,1.0,k1re,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,-1.0,&(Fre->pointer()[0][0]),nmo_,&(Dim->pointer()[0][0]),nmo_,1.0,k1re,nmo_);

        F_DGEMM('n','n',nmo_,nmo_,nmo_,-1.0,&(Dre->pointer()[0][0]),nmo_,&(Fre->pointer()[0][0]),nmo_,0.0,k1im,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,1.0,&(Dim->pointer()[0][0]),nmo_,&(Fim->pointer()[0][0]),nmo_,1.0,k1im,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,1.0,&(Fre->pointer()[0][0]),nmo_,&(Dre->pointer()[0][0]),nmo_,1.0,k1im,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,-1.0,&(Fim->pointer()[0][0]),nmo_,&(Dim->pointer()[0][0]),nmo_,1.0,k1im,nmo_);
        for (int i = 0; i < nmo_; i++) {
            for (int j = 0; j < nmo_; j++) {
                tempre[i*nmo_+j] = Dre->pointer()[i][j] + k1re[i*nmo_+j] * time_step / 2.0;
                tempim[i*nmo_+j] = Dim->pointer()[i][j] + k1im[i*nmo_+j] * time_step / 2.0;
            }
        }
        F_DGEMM('n','n',nmo_,nmo_,nmo_,-1.0,&(Dre->pointer()[0][0]),nmo_,&(Fre->pointer()[0][0]),nmo_,0.0,k3re,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,-1.0,&(Dim->pointer()[0][0]),nmo_,&(Fim->pointer()[0][0]),nmo_,0.0,k1re,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,1.0,&(Fre->pointer()[0][0]),nmo_,&(Dre->pointer()[0][0]),nmo_,0.0,k4re,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,1.0,&(Fim->pointer()[0][0]),nmo_,&(Dim->pointer()[0][0]),nmo_,0.0,k2re,nmo_);

        F_DGEMM('n','n',nmo_,nmo_,nmo_,1.0,&(Dim->pointer()[0][0]),nmo_,&(Fre->pointer()[0][0]),nmo_,0.0,k1re,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,1.0,&(Dre->pointer()[0][0]),nmo_,&(Fim->pointer()[0][0]),nmo_,1.0,k1re,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,-1.0,&(Fim->pointer()[0][0]),nmo_,&(Dre->pointer()[0][0]),nmo_,1.0,k1re,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,-1.0,&(Fre->pointer()[0][0]),nmo_,&(Dim->pointer()[0][0]),nmo_,1.0,k1re,nmo_);

        // k2     = f( t( n + 1/2 h ) , y( n ) + h/2 k1 )
        BuildFock(tempre,tempim,(iter+0.5)*time_step);

        F_DGEMM('n','n',nmo_,nmo_,nmo_,1.0,&(tempim[0]),nmo_,&(Fre->pointer()[0][0]),nmo_,0.0,k2re,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,1.0,&(tempre[0]),nmo_,&(Fim->pointer()[0][0]),nmo_,1.0,k2re,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,-1.0,&(Fim->pointer()[0][0]),nmo_,&(tempre[0]),nmo_,1.0,k2re,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,-1.0,&(Fre->pointer()[0][0]),nmo_,&(tempim[0]),nmo_,1.0,k2re,nmo_);

        F_DGEMM('n','n',nmo_,nmo_,nmo_,-1.0,&(tempre[0]),nmo_,&(Fre->pointer()[0][0]),nmo_,0.0,k2im,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,1.0,&(tempim[0]),nmo_,&(Fim->pointer()[0][0]),nmo_,1.0,k2im,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,1.0,&(Fre->pointer()[0][0]),nmo_,&(tempre[0]),nmo_,1.0,k2im,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,-1.0,&(Fim->pointer()[0][0]),nmo_,&(tempim[0]),nmo_,1.0,k2im,nmo_);

        for (int i = 0; i < nmo_; i++) {
            for (int j = 0; j < nmo_; j++) {
                tempre[i*nmo_+j] = Dre->pointer()[i][j] + k2re[i*nmo_+j] * time_step / 2.0;
                tempim[i*nmo_+j] = Dim->pointer()[i][j] + k2im[i*nmo_+j] * time_step / 2.0;
            }
        }
        // k3     = f( t( n + 1/2 h ) , y( n ) + h/2 k2 )
        BuildFock(tempre,tempim,(iter+0.5)*time_step);

        F_DGEMM('n','n',nmo_,nmo_,nmo_,1.0,&(tempim[0]),nmo_,&(Fre->pointer()[0][0]),nmo_,0.0,k3re,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,1.0,&(tempre[0]),nmo_,&(Fim->pointer()[0][0]),nmo_,1.0,k3re,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,-1.0,&(Fim->pointer()[0][0]),nmo_,&(tempre[0]),nmo_,1.0,k3re,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,-1.0,&(Fre->pointer()[0][0]),nmo_,&(tempim[0]),nmo_,1.0,k3re,nmo_);

        F_DGEMM('n','n',nmo_,nmo_,nmo_,-1.0,&(tempre[0]),nmo_,&(Fre->pointer()[0][0]),nmo_,0.0,k3im,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,1.0,&(tempim[0]),nmo_,&(Fim->pointer()[0][0]),nmo_,1.0,k3im,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,1.0,&(Fre->pointer()[0][0]),nmo_,&(tempre[0]),nmo_,1.0,k3im,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,-1.0,&(Fim->pointer()[0][0]),nmo_,&(tempim[0]),nmo_,1.0,k3im,nmo_);

        for (int i = 0; i < nmo_; i++) {
            for (int j = 0; j < nmo_; j++) {
                tempre[i*nmo_+j] = Dre->pointer()[i][j] + k3re[i*nmo_+j] * time_step;
                tempim[i*nmo_+j] = Dim->pointer()[i][j] + k3im[i*nmo_+j] * time_step;
            }
        }
        // k4     = f( t( n + h )     , y( n ) + h k3 )
        BuildFock(tempre,tempim,(iter+1)*time_step);

        F_DGEMM('n','n',nmo_,nmo_,nmo_,1.0,&(tempim[0]),nmo_,&(Fre->pointer()[0][0]),nmo_,0.0,k4re,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,1.0,&(tempre[0]),nmo_,&(Fim->pointer()[0][0]),nmo_,1.0,k4re,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,-1.0,&(Fim->pointer()[0][0]),nmo_,&(tempre[0]),nmo_,1.0,k4re,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,-1.0,&(Fre->pointer()[0][0]),nmo_,&(tempim[0]),nmo_,1.0,k4re,nmo_);

        F_DGEMM('n','n',nmo_,nmo_,nmo_,-1.0,&(tempre[0]),nmo_,&(Fre->pointer()[0][0]),nmo_,0.0,k4im,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,1.0,&(tempim[0]),nmo_,&(Fim->pointer()[0][0]),nmo_,1.0,k4im,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,1.0,&(Fre->pointer()[0][0]),nmo_,&(tempre[0]),nmo_,1.0,k4im,nmo_);
        F_DGEMM('n','n',nmo_,nmo_,nmo_,-1.0,&(Fim->pointer()[0][0]),nmo_,&(tempim[0]),nmo_,1.0,k4im,nmo_);

        // y(n+1) = y( n ) + 1/6 h ( k1 + 2k2 + 2k3 + k4 )

        for (int i = 0; i < nmo_; i++) {
            for (int j = 0; j < nmo_; j++) {
                Dre->pointer()[i][j] += 1.0 / 6.0 * time_step * ( k1re[i*nmo_+j] + 2.0 * k2re[i*nmo_+j] + 2.0 * k3re[i*nmo_+j] + k4re[i*nmo_+j] );
                Dim->pointer()[i][j] += 1.0 / 6.0 * time_step * ( k1im[i*nmo_+j] + 2.0 * k2im[i*nmo_+j] + 2.0 * k3im[i*nmo_+j] + k4im[i*nmo_+j] );
            }
        }



        // evaluate dipole moment:
        double dpre = C_DDOT(nmo_*nmo_,&(Dre->pointer())[0][0],1,&(dipole[0]->pointer())[0][0],1);
        dpre       += C_DDOT(nmo_*nmo_,&(Dre->pointer())[0][0],1,&(dipole[1]->pointer())[0][0],1);
        dpre       += C_DDOT(nmo_*nmo_,&(Dre->pointer())[0][0],1,&(dipole[2]->pointer())[0][0],1);
        double dpim = C_DDOT(nmo_*nmo_,&(Dim->pointer())[0][0],1,&(dipole[0]->pointer())[0][0],1);
        dpim       += C_DDOT(nmo_*nmo_,&(Dim->pointer())[0][0],1,&(dipole[1]->pointer())[0][0],1);
        dpim       += C_DDOT(nmo_*nmo_,&(Dim->pointer())[0][0],1,&(dipole[2]->pointer())[0][0],1);

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
                  + C_DDOT(nmo_*nmo_,&Dre->pointer()[0][0],1,&T->pointer()[0][0],1)
                  + C_DDOT(nmo_*nmo_,&Dre->pointer()[0][0],1,&V->pointer()[0][0],1);
        en += enuc_;

        outfile->Printf("@TDHF TIME %20.12lf %20.12lf %20.12lf %20.12lf %20.12lf\n",iter*time_step, dpre, dpim, field, en);

    }

    free(factorb);
    free(factorB);

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
