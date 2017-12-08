#include "test.h"
/*
Testing weather reinterpret cast can be used on a vector.
*/

int main()
{


vector<complex<double> > test(5,ZERO);

for(int i=0;i<2;i++)
{
	test[i] += ONE;
}
for(int i=2;i<5;i++)
{
	test[i].real() = i;
	test[i].imag() = 2*i;
	test[i] += ONE;
}

complex<double> w (2,2);
complex<double> e (1,1);
w = 2.0 * w;
printf("%f + %fi\n",w.real(),w.imag());
//test1(test);
/*
for(int i=0;i<5;i++)
{
	printf("%f + %fi the absolute value of which is : %f\n",test[i].real(),test[i].imag(),abs(test[i]));
}
*/
}

void test1(vector<complex<double> > p)
{
	for(int i=0;i<3;i++)
	{
		printf("%f\n",p[i].imag());
	}
}