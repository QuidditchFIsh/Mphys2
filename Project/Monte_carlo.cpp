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

	11/02/18-Note	
	Changing the algorithm a little because you can half step either p or q first 
	But for this fourier transformed version we will need to calculate the new p's before we 
	update the q's using the p's. We will also need to fourier transform into and out of
	k space here.

	The new hamiltonian is 

	H =0.5 * p^2/(1-d^2) + lamba*(q^2 - f^2)^2


*/
#include "Monte_carlo.h"


void lattice_Evolution(unsigned int length,double t_step,unsigned int iterations,double mu,double lamba,double m,double a,double f)
{

printf("##########################\n");
printf("\n");



printf("Running in Anharmonic Mode\n");


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

//perhaps find a different way to do this.

 	vector<normal_distribution<double> > generators;
 	double varience =0;

 	for(unsigned int j = 0;j<length;j++)
 	{
 		varience = 1/ ( 1 + ((4 * m / a) * sin((PI/length) * j) * sin((PI/length) * j)));
 		normal_distribution<double> distribution(0.0,varience);
 		generators.push_back( distribution);
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
	ProgressBar progressBar(iterations, 70);
//run main algorithm
 	for(unsigned int i = 0; i<iterations;i++)
 	{
 			++progressBar;
 			progressBar.display();
 		
 		//default_random_engine generator(random_device{}());
 		for(unsigned int j = 0; j<length;j++)
 		{
 			p[j] = generators[j](generator) * ONE;
 		}

//Start the main algorithm 
 
 		acceptance += hmcAlgorithm(length,t_step,mu,steps,delta_H_Average,m,a,p,q,p_temp,q_temp,f);

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

#if Output_Stats

 	 	fprintf(output_stats,"%d %f %f %f %f %f \n",i,temp_avgx,delta_H_Average,temp_avgx2,lattice_Action(q,length,m,a,mu,lamba),lattice_KineticEnergy(p,length));

#endif

 	 	}
#if Output_X
 		for(unsigned int l=0;l<length;l++)
		{
 			fprintf(output_X,"%f ",q[l].real());
 		}
 		fprintf(output_X,"\n");
#endif
 		

 	}
 	progressBar.done();
 	


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

double hmcAlgorithm(unsigned int length,double t_step,double mu,unsigned int steps,double &delta_H_Average,double m ,double a,vector<complex<double> > &p,vector<complex<double> > &q,vector<complex<double> > &p_temp,vector<complex<double> > &q_temp,double f)
{

	double min=0,H_old=0,H_new=0;


	backwardTransform(p,length);

	H_old=lattice_Hamiltonian(p,q,length,mu,1.0,m,a,f);

	forwardTransform(p,length);

//calculate the new p's
for (unsigned int j=0;j<length;j++)
{
	p_temp[j] = p[j] / ((4 * m / a) * sin((PI/length) * j) * sin((PI/length) * j));
	q_temp[j] = q[j];
}

	backwardTransform(p,length);

// half update the q's
for(unsigned int j=0;j<length;j++)
{
	q_temp[j] = q_temp[j] + (0.5 * t_step * p_temp[j]);
}

for(unsigned int k = 0;k < steps;k++)
{
//update the p's
	for(unsigned int j=0;j<length;j++)
	{
		p_temp[j] = p_temp[j] - 4 * t_step * a * 1 * q_temp[j] * ((q_temp[j] * q_temp[j]) - f);
	}

	if( k != steps-1)
	{
		forwardTransform(p,length);
	

//calculate the new p's
		for(unsigned int j=0;j<length;j++)
		{
			p_temp[j] = p_temp[j] / ((4 * m / a) * sin((PI/length) * j) * sin((PI/length) * j));
		
		}

			backwardTransform(p,length);

// update the q's

		for(unsigned int j=0;j<length;j++)
		{
			q_temp[j] = q_temp[j] + (t_step * p_temp[j]);
		}
	}
}

	forwardTransform(p,length);

	//calculate the new p's
	for (unsigned int j=0;j<length;j++)
	{
		p_temp[j] = p_temp[j] / ((4 * m / a) * sin((PI/length) * j) * sin((PI/length) * j));

	}

	backwardTransform(p,length);

	// half update the q's
	for(unsigned int j=0;j<length;j++)
	{
		q_temp[j] = q_temp[j] + (0.5 * t_step * p_temp[j]);
	}

	H_new = lattice_Hamiltonian(p_temp,q_temp,length,mu,1.0,m,a,f);


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



