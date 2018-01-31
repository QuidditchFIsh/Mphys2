/*
	Author: Aneirin John Baker 
	Date: 18/12/17
	Description: Script to hold all of the potentials and kinetic energies to be used in the alforithm. Here for ease of use and editing. 
*/

#include "functions.h"

double Harmonic_Potential(double q,double q_plus,double q_minus,double m,double mu,double a)
{
	//usual harmonic potential
	return (a * mu * q) + ((m/a)*(2*q - q_minus - q_plus));
}
double Anharmonic_Potential_f(double q,double q_plus,double q_minus,double m,double a,double lamba,double f)
{
	//anharmonic potential with lamba and f
	return (4 * a * lamba * q * (pow(q,2)-f)) + ((m/a)*(2*q - q_minus - q_plus));
}