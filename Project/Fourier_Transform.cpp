#include "Fourier_Transform.h"


fftw_complex forwardTransform(fftw_complex in,int length)
{
	//define array which is to be outputted
	fftw_complex out[length];

	//create plan
	fftw_plan forward = fftw_plan_dft_1d(length,in,out,FFTW_FORWARD,FFTW_ESTIMATE);

	//execute fourier transform 
	fftw_execute(forward);

	//clean up
	fftw_destroy_plan(forward);

	return out;


}
fftw_complex backwardTransform(fftw_complex in,int length)
{
	//define array which is to be outputted
	fftw_complex out[length];

	//create plan
	fftw_plan backward = fftw_plan_dft_1d(length,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);

	//execute fourier transform 
	fftw_execute(backward);

	//clean up
	fftw_destroy_plan(backward);

	return out;
}