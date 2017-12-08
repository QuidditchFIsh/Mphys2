#include "Fourier_Transform.h"


fftw_complex forwardTransform(vector<complex<double> > in_vector,int length)
{

	//define array which is to be outputted
	fftw_complex out[length];
	fftw_complex in[length];

	for(int i=0;i<length;i++)
	{
		in[i][REAL] = in_vector[i].real()
		in[i][IMAG] = in_vector[i].imag();
	}
	//create plan
	fftw_plan forward = fftw_plan_dft_1d(length,in,out,FFTW_FORWARD,FFTW_ESTIMATE);

	//execute fourier transform 
	fftw_execute(forward);

	//clean up
	fftw_destroy_plan(forward);

	return out;


}
fftw_complex backwardTransform(vector<complex<double> > in_vector,int length)
{
{
	//define array which is to be outputted
	fftw_complex out[length];
	fftw_complex in[length];

	for(int i=0;i<length;i++)
	{
		in[i][REAL] = in_vector[i].real()
		in[i][IMAG] = in_vector[i].imag();
	}

	//create plan
	fftw_plan backward = fftw_plan_dft_1d(length,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);

	//execute fourier transform 
	fftw_execute(backward);

	//clean up
	fftw_destroy_plan(backward);

	return out;
}