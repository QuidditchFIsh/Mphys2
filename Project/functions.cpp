/*
Author : Aneirin John Bake r
Date: 31/01/2018
Description: Script to hold all of the update method formule so that they can easily be editted.
For the Fourier transfromed algorithm.
*/

#include "functions.h"

complex<double> Harmonic_Potential(complex<double>  q,double m,double mu,double a,double length,double j)
{
	//usual harmonic potential
	return q * ((mu * a) + ((m/a) * 4 * sin((PI/length) * j) * sin((PI/length) * j)));
}
complex<double> Anharmonic_Potential(complex<double> q,double m,double a,double lamba,double f,double j,double length)
{
	//anharmonic potential with lamba and f
	return q * (((m/a) * 4 * sin((PI/length) * j) * sin((PI/length) * j)) + (4 * lamba * a * ((q*q)- f)));
}