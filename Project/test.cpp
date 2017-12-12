#include "test.h"
/*
Testing weather reinterpret cast can be used on a vector.
*/

int main()
{

default_random_engine generator(random_device{}());
 normal_distribution<double> distribution(0.0,1.0);

uniform_real_distribution<double> Udistribution(0.0,1.0);
vector<complex<double> > test(5,ZERO);

for(int i=0;i<2;i++)
{
	test[i] += ONE;
}
for(int i=2;i<5;i++)
{
	test[i] = Udistribution(generator)*ONE;
	//test[i] += ONE;
} 
for(int i=0;i<5;i++)
{
	printf("%f + %fi\n",test[i].real(),test[i].imag());
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