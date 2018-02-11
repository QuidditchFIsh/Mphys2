#include "test.h"
/*
Testing weather reinterpret cast can be used on a vector.
*/

int main()
{

	default_random_engine generator;

	double varience =0,length=5,PI=3.141,a=1,m=1;

	vector<normal_distribution<double> > generators;

 	for(unsigned int j = 0;j<length;j++)
 	{
 		varience = 1/ ( 1 + ((4 * m / a) * sin((PI/length) * j) * sin((PI/length) * j)));
 		normal_distribution<double> distribution(0.0,varience);
 		generators.push_back( distribution);
 	}

 	printf("%f\n",generators[0](generator));
 	printf("%f\n",generators[3](generator));
 	printf("%f\n",generators[4](generator));


}

void test1(vector<complex<double> > p)
{
	for(int i=0;i<3;i++)
	{
		printf("%f\n",p[i].imag());
	}
}