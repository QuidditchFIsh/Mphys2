#include <fftw3.h>
#include <vector>
#include <complex>
#include "constants.h"

using namespace std;


void forwardTransform(vector<complex<double> > &,int,double);

void backwardTransform(vector<complex<double> > &,int ,double);
