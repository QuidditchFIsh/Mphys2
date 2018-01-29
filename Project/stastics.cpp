/*
	Author: Aneirin John Baker 
	Date: 23/09/17
	Description: This script is where the functions to run the stastics are housed. It will include methods
		which will calculate the mean and other such stastical variables.
*/

#include "stastics.h"



double avgX(vector<double> results)
{	
	double sum =0;
	int length = results.size();

	for(int i=0;i<length;i++)
	{
		sum += results[i];
	}

	return sum/length;
}

double avg_X_Sqd(vector<double> results)
{	

	double sum =0;
	unsigned int length = results.size();

		for(unsigned int j=0;j<length;j++)
		{
			sum += results[j]*results[j];
		}

	return sum/(double)length;
}

double avg_X_four(vector<double> results)
{

	double sum =0;

	unsigned int length = results.size();

		for(unsigned int j=0;j<length;j++)
		{
			sum += results[j]*results[j]*results[j]*results[j];
		}

	return sum/(double)length;
}


double standard_Deviation(double avg_X_Sqd, double avgX,double length )
{
	double std_dev = (avg_X_Sqd - pow(avgX,2)) / (length -1.0);

	return sqrt(std_dev);
}

double lattice_Hamiltonian(vector<vector<double> > state,unsigned int length,double mu,double lamba,double m,double a,double f)
{
	double H=0;
	//loop for all sites which are not effected by periodic BC's
#if !Oscillator_flip
	
	for(unsigned int i=0;i<length-1;i++)
	{
		H += Harmonic_hamiltonian(state[0][i],state[1][i],state[1][i+1],mu,m,a);
	}
	//Periodic BC sites
	H += Harmonic_hamiltonian(state[0][length-1],state[1][length-1],state[1][0],mu,m,a);

#endif

#if Oscillator_flip
	for(unsigned int i=0;i<length-1;i++)
	{
		H += Anarmonic_hamiltonian(state[0][i],state[1][i],state[1][i+1],1,m,a,f);
	}
	//Periodic BC sites
	H += Anarmonic_hamiltonian(state[0][length-1],state[1][length-1],state[1][0],1,m,a,f);
#endif

	return H;

}
double lattice_Action(vector<double> q,unsigned int length,double m,double a,double mu,double lamba)
{
	double S = 0;
#if !Oscillator_flip
	for(unsigned int i =0; i<length-1;i++)
	{
		S += Harmonic_action(q[i],q[i+1],m,a,mu);
	}

	S += Harmonic_action(q[length-1],q[0],m,a,mu);
#endif

#if Oscillator_flip
	for(unsigned int i =0; i<length-1;i++)
	{
		S += Anarmonic_action(q[i],q[i+1],m,a,mu,lamba);
	}

	S += Anarmonic_action(q[length-1],q[0],m,a,mu,lamba);
#endif

	return S/(double)length;
}

double lattice_KineticEnergy(vector<double> p,unsigned int length)
{
	double KE =0;

	for(unsigned int i =0; i<length;i++)
	{
		KE += kinetic_Energy(p[i]);
	}

	return KE/(double)length;
}


double Harmonic_hamiltonian(double p,double q,double q_plus,double mu,double m,double a)
{
	return (p*p*0.5) + ((pow((q_plus - q),2)*0.5*(m/a)) + (a*mu*0.5*q*q));
}

double Harmonic_action(double q, double q_plus,double m,double a,double mu)
{
	return ((pow((q_plus - q),2)*0.5*(m/a)) + (a*mu*0.5*q*q));
}
double Anarmonic_hamiltonian(double p,double q,double q_plus ,double lamba,double m,double a,double f)
{
	return (p*p*0.5) + (pow((q_plus - q),2)*0.5*(m/a)) + (a*lamba * ((q*q) - f) * ((q*q) - f));
}

double Anarmonic_action(double q, double q_plus,double m,double a,double mu,double lamba)
{
	return (pow((q_plus - q),2)*0.5*(m/a)) + (a*lamba * ((q*q)-4) * ((q*q)-4));
}
double kinetic_Energy(double p)
{
	return (p * p * 0.5);
}
