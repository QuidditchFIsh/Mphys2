#include <cmath>
#include <complex>
#include <vector>
//#include "constants.h"

using namespace std;

#define I complex<double>(0,1)
#define ONE complex<double>(1,0)
#define ZERO complex<double>(0,0)
#define PI2 6.28318530718

void Forward_Transform(vector<complex<double> > &,unsigned int ,double);
void Backward_Transform(vector<complex<double> > &,unsigned int ,double);