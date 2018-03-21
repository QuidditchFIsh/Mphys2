/*
	Author: Aneirin John Baker 
	Date: 23/09/17
	Description: The script is where the HMC algorithm will take place. The intergration methods will be 
	housed here and will be executed here and all of the stats functions will be executed here to create the Raw stats data.

	Fourier Accelerated Branch so only will be working with harmonic oscillator for now.

	will have to change from vectors to fftw_complex for the fourier acceleration to work. 
	Will now have to add another fftw_complex to the simulation as a 2D FFT works differently from
	2 1d ones. 

	Had to change from using vectors to doubles due to need of re interpret cast *

	*may actually have to change to complex class as the fourier transform will create complex numbers. 

	Will always need to have doubles in the multipication otherwise there will be errors thrown. 

	Notes:
	p-0,q-1


*/
#include "Monte_carlo.h"


void lattice_Evolution(unsigned int length,double t_step,unsigned int iterations,double mu,double lamba,double m,double a,double f)
{

printf("##########################\n");
printf("\n");

#if Oscillator_flip
printf("Running in Harmonic Mode\n");
#endif

#if !Oscillator_flip
printf("Running in Anharmonic Mode\n");
#endif

	FILE * output_stats;
	output_stats = fopen("HMC_Stats.dat","w");


	FILE * output_X;
	output_X = fopen("HMC_X.dat","w");



//create vectors of complex numbers one for each p and q for convience and clarity in the fourier transform 
	vector<complex<double> > p(length,ZERO);
	vector<complex<double> > q(length,ZERO);
	vector<complex<double> > p_temp(length,ZERO);
	vector<complex<double> > q_temp(length,ZERO);

	default_random_engine generator(random_device{}());
 	normal_distribution<double> distribution(0.0,1.0);

 	uniform_real_distribution<double> Udistribution(0.0,1.0);

 	vector<normal_distribution<double> > generators;
 	double varience =0;


//check there is now way to make this more efficent ie don't look over all of l only loop over half of it !!!
 	for(unsigned int j=0;j<length;j++)
 	{

 		varience = 1/ ( 1 + ((4/(a*a)) * sin((PI/length) * j) * sin((PI/length) * j)));
 		normal_distribution<double> distribution(0.0,varience);
 		generators.push_back(distribution);
 	}


 	double acceptance =0,delta_H_Average=0,avgx=0,avgx2=0,temp_avgx=0,temp_avgx2=0,temp_avgx4=0,avgx4=0,dH_avg=0;

 	unsigned int steps =20,burn=0;

//initalise the first state of the siulation 
 	for(unsigned int j=0;j<length;j++)
 	{ 
 	q[j] = Udistribution(generator) * ONE;
 		if( j % 2 == 0)
 		{
 		 	//state[1][j]= state[1][j] * -1;
 			//q[j] *= -1.0 * ONE;
 		}
	}
//run main algorithm
 	for(unsigned int i = 0; i<iterations;i++)
 	{
 		default_random_engine generator(random_device{}());

 		//NEED TO ONLY DO THIS FOR L/2	
 		for(unsigned int j = 0; j<length/2;j++)
 		{
 			//p[j] = distribution(generator) * ONE;
 			p[j] = (generators[j](generator) * ONE) + (generators[j](generator) * I);
 		}


 		for(unsigned int j=0; j<length/2 ;j++)

 		{
 			p[length/2 + j] = (p[j].real() * ONE) - (p[j].imag() * I);
 		}
 		// for(int j=0;j<length;j++)
 		// {
 		// 	printf("%f %f ",p[j].real(),p[j].imag());
 		// }
 		// printf("\n");

//Start the main algorithm 
 

 		acceptance += hmcAlgorithm(length,t_step,1.0,steps,delta_H_Average,1.0,1.0,p,q,p_temp,q_temp,f);
 		if(i%10 ==0)
 		{
 			printf("%d\n",i);
 		}

//perform the stats calculations for the raw data

		temp_avgx = avgX(q);
		temp_avgx2 = avg_X_Sqd(q);
		temp_avgx4 = avg_X_four(q);
		dH_avg += delta_H_Average;

		if(i>burn)
		{
 		avgx +=temp_avgx;
 		avgx2 +=temp_avgx2;
 		avgx4 += temp_avgx4;
 	 	fprintf(output_stats,"%d %f %f %f %f %f \n",i,temp_avgx,delta_H_Average,temp_avgx2,lattice_Action(q,length,m,a,mu,lamba),lattice_KineticEnergy(p,length));
 	 	}
 	//  	for(unsigned int l=0;l<length;l++)
		// {
 	//  		fprintf(output_X,"%f ",q[l].real());
 	//  	}
 	//  	fprintf(output_X,"\n");
 		

 	}
 	
 	// 	for(unsigned int l=0;l<length;l++)
		// {
 	// 		fprintf(output_X1,"%f\n",q[l].real());
 	// 	}
 	// 	fprintf(output_X1,"\n");


 	double stdx=0,stdx2=0;

 	stdx  = sqrt(((avgx2 / (iterations - burn)) - pow(avgx  / (iterations - burn),2)) / (iterations - burn - 1));
 	stdx2 = sqrt(((avgx4 / (iterations - burn)) - pow(avgx2 / (iterations - burn),2)) / (iterations - burn - 1));

//Output the Data to the Terminal To save Calcuation time in Python
 	printf("########## Data ##########\n");
 	printf("\n");
 	printf("Acceptance: %f%%\n",(acceptance*100)/(double) iterations);
 	printf("The Average Delta H was %f\n",dH_avg/iterations);
 	printf("Average x:	%f +/-%f\n",avgx/(iterations-burn),stdx);
 	printf("Average x^2: %f +/-%f\n",avgx2/(iterations-burn),stdx2);
 	printf("Average x^4: %f\n",avgx4/(iterations-burn));
 	double GroundState=0;

 	GroundState = (mu*avgx2/(iterations-burn));

 	printf("Ground State Energy: %f\n",GroundState);


}
/*
double hmcAlgorithm(unsigned int length,double t_step,double mu,unsigned int steps,double &delta_H_Average,double m ,double a,vector<complex<double> > &p,vector<complex<double> > &q,vector<complex<double> > &p_temp,vector<complex<double> > &q_temp,double f)
{

	double min=0,H_old=0,H_new=0,H_inter=0;


	H_old=lattice_Hamiltonian(p,q,length,mu,1.0,m,a,f);

//	Fourier transform the arrays
	//forwardTransform(p,length);

	//backward transform the p's
	backwardTransform(p,length);

	for(unsigned int j = 0;j<length;j++)
	{
		p_temp[j] = p[j] - (0.5 * t_step * Harmonic_Potential(q[j],m,mu,a,length,j));

		q_temp[j] = q[j];
	}

//full step in p and q for n steps
	for(unsigned int i = 0;i<steps;i++)
	{
//update all q's
		//transform to fourie soace 
		forwardTransform(p_temp,length);

			//CALCULATE THE P's IN FOURIER SPACE
		for(int j=0;j<length;j++)
		{
			p_temp[j] = p_temp[j] / (1 + ((4/a) * sin((PI/length) * j) * sin((PI/length) * j)));
		}

		backwardTransform(p_temp,length);

	//transform them back to positon space to update the 
		for(unsigned int j = 0;j<length;j++)
		{
			q_temp[j] = q_temp[j] + (t_step * p_temp[j]);
		}

		if(i != steps-1)
		{
			for(unsigned int j = 0;j<length;j++)
			{
				p_temp[j] = p_temp[j] - (t_step  * Harmonic_Potential(q_temp[j],m,mu,a,length,j));

			}
		} 	

	}
//half step in the p
	for(unsigned int j = 0;j<length;j++)
	{
		p_temp[j] = p_temp[j] - (0.5 * t_step  * Harmonic_Potential(q_temp[j],m,mu,a,length,j));
	}

//#########backward fourier transform goes here#############

	//backwardTransform(p_temp,length);
	//backwardTransform(p,length);

	forwardTransform(p_temp,length);

	H_new = lattice_Hamiltonian(p_temp,q_temp,length,mu,1.0,m,a,f);
	//printf("%f %f\n",H_old,H_new);
//metroplis update
	double r = ((double) rand() / (RAND_MAX));

	min = (1 < exp(H_old - H_new)) ? 1 : exp(H_old - H_new);
	if(r < min)
	{
//accept
		for(unsigned int i = 0;i<length;i++)
		{
			q[i] = q_temp[i];
			
		}

		delta_H_Average= H_old - H_new;

		return 1;
		
	}
	delta_H_Average = H_old - H_new;

	return 0;
	
}
*/

double hmcAlgorithm(unsigned int length,double t_step,double mu,unsigned int steps,double &delta_H_Average,double m ,double a,vector<complex<double> > &p,vector<complex<double> > &q,vector<complex<double> > &p_temp,vector<complex<double> > &q_temp,double f)
{

	double min=0,H_old=0,H_new=0;


// 	printf("Before\n");
// 	printf("p: ");
// for(int j=0;j<length;j++)
// {
// 	printf("%f %f ",p[j].real(),p[j].imag());
// }
// printf("\n");
// printf("q: ");
// for(int j=0;j<length;j++)
// {
// 	printf("%f %f ",q[j].real(),q[j].imag());
// }
// printf("\n");
// printf("############################\n");


	backwardTransform(p,length,a);

	H_old=lattice_Hamiltonian(p,q,length,mu,1.0,m,1.0,f);
	//printf("%f\n",H_old);

	forwardTransform(p,length,a);

//calculate the new p's
for (unsigned int j=0;j<length;j++)
{
	p_temp[j] = abs(p[j]) / (1 + ((4/(a*a)) * sin((PI/length) * j) * sin((PI/length) * j)));
	q_temp[j] = q[j];
}

	backwardTransform(p_temp,length,a);

// half update the q's
// 		printf("p: ");
// for(int j=0;j<length;j++)
// {
// 	printf("%f %f ",p[j].real(),p[j].imag());
// }
// printf("\n");
for(unsigned int j=0;j<length;j++)
{
	q_temp[j] = q_temp[j] + (0.5 * t_step * p_temp[j]);
}



for(unsigned int k = 0;k < steps;k++)
{
//update the p's
	for(unsigned int j=0;j<length;j++)
	{
		//p_temp[j] = p_temp[j] - 4 * t_step * a * 1 * q_temp[j] * ((q_temp[j] * q_temp[j]) - f);
		p_temp[j] = p_temp[j] - (t_step* mu * a * q_temp[j]);
	}

	if( k != steps-1)
	{
		forwardTransform(p_temp,length,a);
	

//calculate the new p's
		for(unsigned int j=0;j<length;j++)
		{
			p_temp[j] = abs(p_temp[j]) / (1 + ((4/(a*a)) * sin((PI/length) * j) * sin((PI/length) * j)));
		
		}

			backwardTransform(p_temp,length,a);

//update the q's

		for(unsigned int j=0;j<length;j++)
		{
			q_temp[j] = q_temp[j] + (t_step * p_temp[j]);
		}
	}
// 	printf("p:");
// 	for(int j=0;j<length;j++)
// {
// 	printf("%f %f ",p_temp[j].real(),p_temp[j].imag());
// }
// printf("\n");
// printf("q: ");
// for(int j=0;j<length;j++)
// {
// 	printf("%f %f ",q_temp[j].real(),q_temp[j].imag());
// }
// printf("\n");
// printf("############################\n");
}

	forwardTransform(p_temp,length,a);

	//calculate the new p's
	for (unsigned int j=0;j<length;j++)
	{
		p_temp[j] = abs(p_temp[j]) / (1 + ((4 /(a*a)) * sin((PI/length) * j) * sin((PI/length) * j)));

	}

	backwardTransform(p_temp,length,a);

	// half update the q's
	for(unsigned int j=0;j<length;j++)
	{
		q_temp[j] = q_temp[j] + (0.5 * t_step * p_temp[j]);
	}
	forwardTransform(p_temp,length,a);

	H_new = lattice_Hamiltonian(p_temp,q_temp,length,mu,1.0,m,1.0,f);
	//printf("%f\n",H_new);
	//printf("##################\n");

//metroplis update
	double r = ((double) rand() / (RAND_MAX));


// 	printf("After\n");
// 	printf("p: ");
// for(int j=0;j<length;j++)
// {
// 	printf("%f %f ",p_temp[j].real(),p_temp[j].imag());
// }
// printf("\n");
// printf("q: ");
// for(int j=0;j<length;j++)
// {
// 	printf("%f %f ",q_temp[j].real(),q_temp[j].imag());
// }
// printf("\n");
// printf("############################\n");
// printf("-----------------------\n");



	min = (1 < exp(H_old - H_new)) ? 1 : exp(H_old - H_new);
	if(r < min)
	{
//accept
		for(unsigned int i = 0;i<length;i++)
		{
			q[i] = q_temp[i];
			
		}

		delta_H_Average= H_old - H_new;

		return 1;
		
	}
	delta_H_Average = H_old - H_new;

	return 0;
	
}