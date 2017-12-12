/*
	Author: Aneirin John Baker 
	Date: 23/09/17
	Description: This script is where the functions to run the stastics are housed. It will include methods
		which will calculate the mean and other such stastical variables.

	Editting it in this new branch so it have the right functionality to take into account 
	the new complex<double> structuer.Also don't need periodic boundary conditions any more.


*/

#include "stastics.h"



double avgX(vector<complex<double> > results)
{	
	double sum =0;
	int length = results.size();

	for(int i=0;i<length;i++)
	{
		sum += results[i].real();
	}

	return sum/length;
}

double avg_X_Sqd(vector<complex<double> > results)
{	

	double sum =0;
	unsigned int length = results.size();

		for(unsigned int j=0;j<length;j++)
		{
			sum += results[j].real() * results[j].real();
		}

	return sum/(double)length;
}

double avg_X_four(vector<complex<double> > results)
{

	double sum =0;

	unsigned int length = results.size();

		for(unsigned int j=0;j<length;j++)
		{
			sum += results[j].real() * results[j].real() * results[j].real() * results[j].real();
		}

	return sum/(double)length;
}


double standard_Deviation(double avg_X_Sqd, double avgX,double length )
{
	double std_dev = (avg_X_Sqd - pow(avgX,2)) / (length -1.0);

	return sqrt(std_dev);
}

double lattice_Hamiltonian(vector<complex<double> > p,vector<comeplx<double> > q,unsigned int length,double mu,double lamba,double m,double a)
{
	double H=0;
	//loop for all sites which are not effected by periodic BC's
#if Oscillator_flip
	
	for(unsigned int i=0;i<length;i++)
	{
		H += Harmonic_hamiltonian(p[i].real(),q[i].real(),q[i+1].real(),mu,m,a);
	}


#endif

#if !Oscillator_flip
	for(unsigned int i=0;i<length;i++)
	{
		H += Anarmonic_hamiltonian(p[i].real(),q[i].real(),q[i+1].real(),mu,lamba,m,a);
	}
#endif

	return H;

}
double lattice_Action(vector<complex<double> > q,unsigned int length,double m,double a,double mu,double lamba)
{
	double S = 0;
#if Oscillator_flip
	for(unsigned int i =0; i<length;i++)
	{
		S += Harmonic_action(q[i].real(),q[i+1].real(),m,a,mu);
	}
#endif

#if !Oscillator_flip
	for(unsigned int i =0; i<length;i++)
	{
		S += Anarmonic_action(q[i].real(),q[i+1].real(),m,a,mu,lamba);
	}
#endif

	return S/(double)length;
}

double lattice_KineticEnergy(vector<complex<double> > p,unsigned int length)
{
	double KE =0;

	for(unsigned int i =0; i<length;i++)
	{
		KE += kinetic_Energy(p[i].real());
	}

	return KE/(double)length;
}


double Harmonic_hamiltonian(double p,double q,double q_plus,double mu,double m,double a)
{
	//regular version
	return (p*p*0.5*(1/m)) + ((pow((q_plus - q),2)*0.5*(m/a)) + (a*mu*0.5*q*q));
}

double Harmonic_action(double q, double q_plus,double m,double a,double mu)
{
	return ((pow((q_plus - q),2)*0.5*(m/a)) + (a*mu*0.5*q*q));
}
double Anarmonic_hamiltonian(double p,double q,double q_plus ,double mu,double lamba,double m,double a)
{
	return (p*p*0.5*(1/m)) + (pow((q_plus - q),2)*0.5*(m/a)) + (a*lamba * ((q*q)-f)* ((q*q-f)));
	//return (p*p*0.5*(1/m)) + ((pow((q_plus - q),2)*0.5*(m/a)) + (a*lamba*pow(q,4)) + (a*mu*0.5*pow(q,2)));
}

double Anarmonic_action(double q, double q_plus,double m,double a,double mu,double lamba)
{
	return (pow((q_plus - q),2)*0.5*(m/a)) + (a*lamba * ((q*q)-f)* ((q*q)-f));
}
double kinetic_Energy(double p)
{
	return (p * p * 0.5);
}
