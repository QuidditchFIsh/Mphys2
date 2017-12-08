#include <fftw3.h>


#define REAL 0
#define IMAG 1

fftw_complex forwardTransform(vector<complex<double> > ,int );

fftw_complex backwardTransform(vector<complex<double> > ,int );