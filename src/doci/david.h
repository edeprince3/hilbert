#ifndef DAVID_H
#define DAVID_H

namespace psi{

// callback function
typedef void (*CallbackType)(size_t N, size_t maxdim,double **sigma, double **b, void * data);

int david_direct_redo(double *Adiag, int N, int M, double *eps, double **v, double cutoff, int print, CallbackType function, size_t & iter, void * data, int maxdim);

size_t david_direct(double *Adiag, size_t N, size_t M, double *eps, double **v, double cutoff,size_t print, CallbackType function, size_t & iter, void * data);
size_t david_in_core(double **A, size_t N, size_t M, double *eps, double **v,double cutoff, size_t print, size_t & iter);

}

#endif
