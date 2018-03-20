/*
Author:Aneirin John Baker
Deate: 20/03/2018
Description:Basic slow fourier transfrom class to make sure that I know what i am doing. It will be taken from 
https://www2.ph.ed.ac.uk/~s0948358/mysite/Intro-LQCD-basics.pdf
*/

#include "FourierTransform.h"

int main()
{
	vector<complex<double> > test(10,ONE);

	for (int i=0;i<10;i++)
	{
		printf("%f %f\n",test[i].real(),test[i].imag());
	}
	Forward_Transform(test,10,1);
	// for (int i=0;i<10;i++)
	// {
	// 	printf("%f %f\n",test[i].real(),test[i].imag());
	// }

	return 1;
}

void Forward_Transform(vector<complex<double> > &input,unsigned int length,double a)
{
	//create a tempoary complex vector to keep it in
	double real =0,complex=0,argument=0;
	double N = (double) length;

	//loop over the length of the latice for the fourier transform wiriting everything in terms of cos and sine 
	for(unsigned int i = 0;i < length;i++)
	{
		for(unsigned int j = 0;j<length;j++)
		{
			argument = PI2 * i * j * a *(1/N);
			real    += input[j].real() * cos(argument) + input[j].imag() * sin(argument);
			complex += input[j].imag() * cos(argument) - input[j].real() * sin(argument);
		}
		real *= a;
		complex *= a;
		printf("%f %f\n",real,complex);
		input[i] = (real * (ONE)) + (complex * I);
		real =0;
		complex=0;
	}
}
void Backward_Transform(vector<complex<double> > &input,unsigned int length,double a)
{

}