#include <stdlib.h>
#include <stdio.h>
#include <random>
#include <iostream>
#include <string>
#include <ctime>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include <fftw3.h>

#include "stastics.h"
#include "Fourier_Transform.h"
#include "constants.h"
#include "functions.h"

using namespace std;

void lattice_Evolution(unsigned int ,double ,unsigned int ,double,double,double,double,double);

double hmcAlgorithm(unsigned int length,double ,double,unsigned int,double &,double,double,vector<complex<double> > &,vector<complex<double> > &,vector<complex<double> > &,vector<complex<double> > &,double,int flag);




