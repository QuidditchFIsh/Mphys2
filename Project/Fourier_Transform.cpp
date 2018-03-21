#include "Fourier_Transform.h"

void forwardTransform(vector<complex<double> > &in_vector,int length,double a)
{
/*
	//define array which is to be outputted
	fftw_complex out[length];
	fftw_complex in[length];

	for(int i=0;i<length;i++)
	{
		in[i][REAL] = in_vector[i].real();
		in[i][IMAG] = in_vector[i].imag();
	}
	//create plan
	fftw_plan forward = fftw_plan_dft_1d(length,in,out,FFTW_FORWARD,FFTW_ESTIMATE);

	//execute fourier transform 
	fftw_execute(forward);

	//clean up
	fftw_destroy_plan(forward);
	
	for(int i=0;i<length;i++)
	{
		out[i][REAL] /= length;
		out[i][IMAG] /= length;
	}

	for(int i=0;i<length;i++)
	{
		in_vector[i] = 	(out[i][REAL] * ONE) + (out[i][IMAG] * I);
	}
*/

	double real =0,complex=0,argument=0;
	double N = (double) length;
	//double a=1.0;
	double tmp1 [2][length] ;

	for(unsigned int i = 0;i < length;i++)
	{
		for(unsigned int j = 0;j<length;j++)
		{
			argument = i * j * (PI2/N) * a;
			real    += (in_vector[j].real() * cos(argument)) - (in_vector[j].imag() * sin(argument));
			complex += (in_vector[j].imag() * cos(argument)) + (in_vector[j].real() * sin(argument));
			//printf("one:%f %f %f\n", cos(argument),sin(argument),(in_vector[j].real() * cos(argument)) );
		}
		//printf("\n");
		real *= a;
		complex *= a;
		//printf("first:%g %g\n",real,complex);
		tmp1[0][i] = real;
		tmp1[1][i] = complex;
		real =0;
		complex=0;
	}
	for(int i=0;i<length;i++)
	{
		in_vector[i] = (tmp1[0][i] * ONE) + (tmp1[1][i] * I);
	}

}
void backwardTransform(vector<complex<double> > &in_vector,int length,double a)
{
	//define array which is to be outputted
	//fftw_complex out[length];
	// fftw_complex in[length];
	// double out[length];

	// for(int i=0;i<length;i++)
	// {
	// 	in[i][REAL] = in_vector[i].real();
	// 	in[i][IMAG] = in_vector[i].imag();
	// }

	// //create plan
	// //fftw_plan backward = fftw_plan_dft_1d(length,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
	// fftw_plan backward = fftw_plan_dft_c2r_1d(length,in,out,FFTW_ESTIMATE);
	
	// //execute fourier transform 
	// fftw_execute(backward);

	// //clean up
	// fftw_destroy_plan(backward);
	
	// for(int i=0;i<length;i++)
	// {
	// 	//out[i][REAL] /= (double)length;
	// 	//out[i][IMAG] /= (double)length;
	// 	out[i] /= (double)length;
	// }
	// for(int i=0;i<length;i++)
	// {
	// 	//in_vector[i] = 	(out[i][REAL] * ONE) + (out[i][IMAG] * I);
	// 	in_vector[i] = out[i] * ONE;

	// }

	double real =0,complex=0,argument=0;
	double N = (double) length;
	//double a=1.0;
	double tmp1 [2][length] ;

	for(unsigned int i = 0;i < length;i++)
	{
		for(unsigned int j = 0;j<length;j++)
		{
			argument = i * j * (PI2/N) * a;
			real    += (in_vector[j].real() * cos(argument)) + (in_vector[j].imag() * sin(argument));
			complex += (in_vector[j].imag() * cos(argument)) - (in_vector[j].real() * sin(argument));
			//printf("one:%f %f %f\n", cos(argument),sin(argument),(in_vector[j].real() * cos(argument)) );
		}
		//printf("\n");
		real *= a;
		complex *= a;
		//printf("first:%g %g\n",real,complex);
		tmp1[0][i] = real;
		tmp1[1][i] = complex;
		real =0;
		complex=0;
	}
	for(int i=0;i<length;i++)
	{
		in_vector[i] = ((1/N) * tmp1[0][i] * ONE) + ((1/N) * tmp1[1][i] * I);
	}
}
