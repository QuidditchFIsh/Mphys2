#include "fftw_test.h"

#define REAL 0
#define IMAG 1
#define PI 3.141592654
/*
######NOTES######
basic recipie for DFT using FFTW:
-create arrays
-create the plan(s)
-execute the DFT
-re normalise the arrays by n as DFT is no normalised
-clean up

since we are not going for speed here FFTW_ESTIMATE will probably suffice here
however it seems that FFTW_MEASURE would be better but would need to look a little 
more into that before working with it.

remember some of the procedures will destroy the input arrays after they have used them
may need to edit some of the code

also i don't think that fftw has support for vectors so will need to spend some time 
re casting most of the program into double[][] or the fftw_complex 

also need to write methods for complex hamiltonians as some of the update equations 
are now complex and hence will need a complex square 

also need to work on how to pass fftw arrays into the stats class for the complex squaring 



*/  
int main()
{

	test1();

	fftw_complex x[20];

	return 1;
}

void test1()
{
#if 0
	//testing combinations of vectors and fftw_complex 
	fftw_complex x[5];

	vector<fftw_complex> State(2,x);

	#endif 
	#if 1

	//define the length of the complex arrays
	int n = 20;

	double T1[20][2] ;

	fftw_complex * in = reinterpret_cast<fftw_complex*>(T1);
	//input artrays

	fftw_complex x[n]; //equivlant to x[n][2]
	fftw_complex y[n];
	fftw_plan p = fftw_plan_dft_1d(n,in,y,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_plan q = fftw_plan_dft_1d(n,y,in,FFTW_BACKWARD,FFTW_ESTIMATE);


	fftw_complex z[n];

	//double in[n];
	fftw_complex out[n];

	double input=0;

	for(int i=0;i<n;i++)
	{
		T1[i][REAL] = i;
	}

	printf("Origional Array:\n");
	for(int i=0;i<n;i++)
	{
		printf("%f +i%f\n",T1[i][REAL],T1[i][IMAG]);
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

	fftw_execute(q);

	printf("Results are for Backward :\n");
	for(int i=0;i<n;i++)
	{
		printf("%f +i%f\n",T1[i][REAL]/20,T1[i][IMAG]/20);
	}


#if 0
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
	#endif
#endif
}

