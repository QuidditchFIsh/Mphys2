#include <vector>
#include <stdio.h>
#include <random>
#include <cmath>
#include <algorithm>
#include <complex>
#include <fftw3.h>
#include "constants.h"

using namespace std;

double avgX(vector<complex<double> > );

double avg_X_Sqd(vector<complex<double> > );

double avg_X_four(vector<complex<double> > );

double standard_Deviation(double , double ,double  );

double lattice_Hamiltonian(vector<complex<double> > ,vector<complex<double> > ,unsigned int ,double ,double ,double ,double ,double);

double lattice_Action(vector<complex<double> > ,unsigned int,double,double,double,double);

double lattice_KineticEnergy(vector<complex<double> >,unsigned int);

double Harmonic_hamiltonian(double ,double ,double,double,double,double);

double Harmonic_action(double,double,double ,double,double);

double Anarmonic_hamiltonian(double ,double ,double,double,double,double,double);

double Anarmonic_action(double,double,double,double,double,double);

double kinetic_Energy(double);