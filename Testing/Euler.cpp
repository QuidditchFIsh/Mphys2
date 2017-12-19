#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>

using namespace std;

int main()
{
	//inital conditions q=0 , p=1
	double q=0,p=1,t_step=0.3;
	int t = 40;

	FILE * output;
	output = fopen("Test_Euler.dat","w");

	for (int i=0;i<t;i++)
	{
		fprintf(output,"%f %f %f\n",p,q,i*t_step);

		p = p - (t_step * q);
		q = q + (t_step * p);
	}

	return 1;
}