#include <fftw3.h>
#include <iostream>
#include <cmath>
#include <math.h>

#define REAL 0
#define IMAG 1
#define PI 3.141592654

int main()
{
	//define the length of the complex arrays
	int n = 20;
	//input artrays
	fftw_plan p = fftw_plan_dft_1d(n,fftw_complex x[n],fftw_complex y[n],FFTW_FORWARD,FFTW_MEASURE);

	//fftw_complex x[n]; //equivlant to x[n][2]
	//fftw_complex y[n];
	fftw_complex z[n];

	double in[n];
	fftw_complex out[n];

	double input=0;

	for(int i=0;i<n;i++)
	{
		x[i][REAL] = i;
	}

	printf("Origional Array:\n");
	for(int i=0;i<n;i++)
	{
		printf("%f +i%f\n",x[i][REAL],x[i][IMAG]);
	}

	//fftw_plan p = fftw_plan_dft_1d(n,x,y,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	//a backward transform will produce a scaled results so you will need to divide 
	//results of the transformation by the length of the array

	printf("Results are for Forward :\n");
	for(int i=0;i<n;i++)
	{
		printf("%f +i%f\n",y[i][REAL],y[i][IMAG]);
	}

	fftw_plan r = fftw_plan_dft_1d(n,y,z,FFTW_BACKWARD,FFTW_ESTIMATE);
	fftw_execute(r);
	for(int i=0;i<n;i++)
	{
		y[i][REAL] /= (double)n;
		y[i][IMAG] /= (double)n;

		z[i][REAL] /= (double)n;
		z[i][IMAG] /= (double)n;
	}


	//clean up
	//fftw_destroy_plan(p);
	//fftw_destroy_plan(r);
//	fftw_cleanup();

	printf("Origional Array:\n");
	for(int i=0;i<n;i++)
	{
		printf("%f +i%f\n",x[i][REAL],x[i][IMAG]);
	}

	printf("Results are for Forward :\n");
	for(int i=0;i<n;i++)
	{
		printf("%f +i%f\n",y[i][REAL],y[i][IMAG]);
	}

	printf("Results are for both directions:\n");
	for(int i=0;i<n;i++)
	{
		printf("%f +i%f\n",z[i][REAL],z[i][IMAG]);
	}

	return 1;
}