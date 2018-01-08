#include "testing.h"

int main()
{

	vector<complex<double> > test(10,ONE);
	test[5] +=ONE;
	//test[6] -=ONE;

	printf("Average of this is %f \n",avg_X_Sqd(test));
	return 1;
}

double avg_X_Sqd(vector<complex<double> > results)
{	

	double sum =0;
	unsigned int length = results.size();

		for(unsigned int j=0;j<length;j++)
		{
			printf("square is: %f\n",results[j].real() * results[j].real());
			sum += results[j].real() * results[j].real();
		}
		printf("sum is: %f\n",sum);

	return sum/(double)length;
}
