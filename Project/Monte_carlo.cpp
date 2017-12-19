/*
	Author: Aneirin John Baker 
	Date: 23/09/17
	Description: The script is where the HMC algorithm will take place. The intergration methods will be 
	housed here and will be executed here and all of the stats functions will be executed here to create the Raw stats data.
*/
#include "Monte_carlo.h"
#define Oscillator_flip 1
<<<<<<< HEAD
//0 = Harmonic
//1 = Anharmonic
=======
//1 = Harmonic
//0 = Anharmonic
#define Anharmonic_flip 0
//1 = lamba and mu
//0 = lamba and f
>>>>>>> 63622eac36ac696e8972a64a076a6b58df1d5301

void lattice_Evolution(vector<vector<double> > &lattice,unsigned int length,double t_step,unsigned int iterations,double mu,double lamba,double m,double a)
{
printf("##########################\n");
printf("\n");
#if Oscillator_flip
	printf("Running in Harmonic Mode\n\n");
#endif

#if !Oscillator_flip
	printf("Running in Anharmonic Mode\n");
#endif


	//FILE * out;
	//out = fopen("HMC_LeapFrog_H","w");

	FILE * output_stats;
	output_stats = fopen("HMC_Stats.dat","w");

	FILE * output_X;
	output_X = fopen("HMC_X.dat","w");

	FILE * output_X1;
	output_X1 = fopen("HMC_X1.dat","w");

	// p-0,q-1
	vector<double> v(length,0);
	vector<vector<double> > State(2,v);
	vector<vector<double> >temp_State(2,v);
	vector<vector<double> >Energy_save(3,v);
	vector<vector<double> >first_state(2,v);

	vector<double> square_state(length,0);

	vector<double> H_store(201,0);
	H_store[0]=0;

	default_random_engine generator(random_device{}());
 	normal_distribution<double> distribution(0.0,1.0);

 	uniform_real_distribution<double> Udistribution(0.0,1.0);

 	double acceptance =0,delta_H_Average=0,avgx=0,avgx2=0,error_x2=0,error_x=0,temp_avgx=0,temp_avgx2=0,temp_avgx4=0,avgx4=0,dH_avg=0;
 	unsigned int steps =10,burn=2000;


//initalise the first state of the siulation 
 	for(unsigned int j=0;j<length;j++)
 	{

 	State[1][j] = Udistribution(generator);
 		if(j % 2 == 0)
 		{
 			State[1][j]= State[1][j] * -1;
 		}
	}

 	//run main algorithm
 	for(unsigned int i = 0; i<iterations;i++)
 	{
 		default_random_engine generator(random_device{}());

 		for(unsigned int j = 0; j<length;j++)
 		{
 			State[0][j] = distribution(generator);
 		}

 		//Maind HMC algorithm.
 		acceptance += hmcAlgorithm_Harmonic(length,t_step,State,temp_State,H_store,mu,steps,delta_H_Average,m,a);

		temp_avgx = avgX(State[1]);
		temp_avgx2 = avg_X_Sqd(State[1]);
		temp_avgx4 = avg_X_four(State[1]);
		dH_avg += delta_H_Average;

		if(i>burn)
		{

 		avgx +=temp_avgx;
 		avgx2 +=temp_avgx2;
 		avgx4 += temp_avgx4;

 		fprintf(output_stats,"%d %f %f %f %f %f %f\n",i,temp_avgx,delta_H_Average,temp_avgx2,error_x2,lattice_Action(State[1],length,m,a,mu,lamba),lattice_KineticEnergy(State[0],length));
 		}
 		for(unsigned int l=0;l<length;l++)
		{
 			fprintf(output_X,"%f ",State[1][l]);
 		}
 		fprintf(output_X,"\n");

 	}
 		for(unsigned int l=0;l<length;l++)
		{
 			fprintf(output_X1,"%f ",State[1][l]);
 		}
 		fprintf(output_X1,"\n");


 	double stdx=0,stdx2=0;

 	stdx= sqrt(((avgx2/(iterations-burn)) - pow(avgx/(iterations-burn),2))/(iterations-burn-1));
 	stdx2= sqrt(((avgx4/(iterations-burn)) - pow(avgx2/(iterations-burn),2))/(iterations-burn-1));

//Output the Data to the Terminal To save Calcuation time in Python
 	printf("########## Data ##########\n");
 	printf("\n");
 	printf("Acceptance: %f%%\n",(acceptance*100)/(double) iterations);
 	printf("The Average Delta H was %f\n",dH_avg/iterations);
 	printf("Average x:	%f +/-%f\n",avgx/(iterations-burn),stdx);
 	printf("Average x^2: %f +/-%f\n",avgx2/(iterations-burn),stdx2);
 	printf("Average x^4: %f\n",avgx4/(iterations-burn));
 	double GroundState=0;
#if Oscillator_flip
 	GroundState = (mu*avgx2/(iterations-burn));
#endif

#if !Oscillator_flip
 	GroundState = (mu*avgx2/(iterations-burn)) + (3 * lamba * (avgx4/(iterations-burn)));
#endif
 	printf("Ground State Energy: %f\n",GroundState);

}

double hmcAlgorithm(unsigned int length,double t_step,vector<vector<double> > &old_state,vector<vector<double> > &temp_State,vector<double> &H_store,double mu,unsigned int steps,double &delta_H_Average,double m ,double a)
{
	double min=0,H_old=0,H_new=0,H_inter=0,f=1;

	H_old=lattice_Hamiltonian(old_state,length,mu,1,m,a);

	//half step in the p
	#if Oscillator_flip

	temp_State[0][0] = old_state[0][0] -  (0.5*t_step * Anharmonic_Potential_f(old_state[1][0],old_state[1][1],old_state[1][length-1],m,a,1,f));

	#endif

	#if !Oscillator_flip

	temp_State[0][0] = old_state[0][0] -  (0.5*t_step * Harmonic_Potential(old_state[1][0],old_state[1][1],old_state[1][length-1],m,mu,a));

	#endif

	temp_State[1][0] = old_state[1][0];

	for(unsigned int j = 1;j<length-1;j++)
	{
		#if Oscillator_flip

		temp_State[0][j] = old_state[0][j] - (0.5*t_step * Anharmonic_Potential_f(old_state[1][j],old_state[1][j+1],old_state[1][j-1],m,a,1,f));

		#endif

		#if !Oscillator_flip

		temp_State[0][j] = old_state[0][j] - (0.5*t_step * Harmonic_Potential(old_state[1][j],old_state[1][j+1],old_state[1][j-1],m,mu,a));

		#endif

		temp_State[1][j] = old_state[1][j];
	}

	#if Oscillator_flip

	temp_State[0][length-1] = old_state[0][length-1] - (0.5*t_step * Anharmonic_Potential_f(old_state[1][length-1],old_state[1][0],old_state[1][length-2],m,a,1,f));

	#endif

	#if !Oscillator_flip

	temp_State[0][length-1] = old_state[0][length-1] - (0.5*t_step * Harmonic_Potential(old_state[1][length-1],old_state[1][0],old_state[1][length-2],m,mu,a));

	#endif

	temp_State[1][length-1] = old_state[1][length-1];

	//full step in p and q for n steps
	for(unsigned int i = 0;i<steps;i++)
	{
		//update all q's
		for(unsigned int j = 0;j<length;j++)
		{
			temp_State[1][j] = temp_State[1][j] + ((t_step/m) * temp_State[0][j]);
		}

//a full step for when running the algorithm normally
		if(i != steps-1)
		{
			#if Oscillator_flip
			temp_State[0][0] = temp_State[0][0] -  (t_step * Anharmonic_Potential_f(temp_State[1][0],temp_State[1][1],temp_State[1][length-1],m,a,1,f));
			#endif

			#if !Oscillator_flip
			temp_State[0][0] = temp_State[0][0] -  (t_step * Harmonic_Potential(temp_State[1][0],temp_State[1][1],temp_State[1][length-1],m,mu,a));
			#endif

			for(unsigned int j = 1;j<length-1;j++)
			{
				#if Oscillator_flip
				temp_State[0][j] = temp_State[0][j] - (t_step * Anharmonic_Potential_f(temp_State[1][j],temp_State[1][j+1],temp_State[1][j-1],m,a,1,f));
				#endif

				#if !Oscillator_flip
				temp_State[0][j] = temp_State[0][j] - (t_step * Harmonic_Potential(temp_State[1][j],temp_State[1][j+1],temp_State[1][j-1],m,mu,a));
				#endif

			}
			#if Oscillator_flip
			temp_State[0][length-1] = temp_State[0][length-1] - (t_step * Anharmonic_Potential_f(temp_State[1][length-1],temp_State[1][0],temp_State[1][length-2],m,a,1,f));
			#endif

			#if !Oscillator_flip
			temp_State[0][length-1] = temp_State[0][length-1] - (t_step * Harmonic_Potential(temp_State[1][length-1],temp_State[1][0],temp_State[1][length-2],m,mu,a));
			#endif
		}

	}
	//half step in the p
	#if Oscillator_flip

	temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * Anharmonic_Potential_f(temp_State[1][0],temp_State[1][1],temp_State[1][length-1],m,a,1,f));

	#endif

	#if !Oscillator_flip

	temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * Harmonic_Potential(temp_State[1][0],temp_State[1][1],temp_State[1][length-1],m,mu,a));

	#endif

			for(unsigned int j = 1;j<length-1;j++)
			{
				#if Oscillator_flip

				temp_State[0][j] = temp_State[0][j] - (0.5*t_step * Anharmonic_Potential_f(temp_State[1][j],temp_State[1][j+1],temp_State[1][j-1],m,a,1,f));

				#endif

				#if !Oscillator_flip

				temp_State[0][j] = temp_State[0][j] - (0.5*t_step * Harmonic_Potential(temp_State[1][j],temp_State[1][j+1],temp_State[1][j-1],m,mu,a));

				#endif

			}
			#if Oscillator_flip

			temp_State[0][length-1] = temp_State[0][length-1] - (0.5*t_step * Anharmonic_Potential_f(temp_State[1][length-1],temp_State[1][0],temp_State[1][length-2],m,a,1,f));

			#endif 

			#if !Oscillator_flip

			temp_State[0][length-1] = temp_State[0][length-1] - (0.5*t_step * Harmonic_Potential(temp_State[1][length-1],temp_State[1][0],temp_State[1][length-2],m,mu,a));

			#endif
	H_new = lattice_Hamiltonian(temp_State,length,mu,1,m,a);

	//metroplis update
	double r = ((double) rand() / (RAND_MAX));

	min = (1 < exp(H_old - H_new)) ? 1 : exp(H_old - H_new);
	if(r < min)
	{
		//accept
		for(unsigned int i = 0;i<length;i++)
		{
			old_state[1][i] = temp_State[1][i];
			old_state[0][i] = temp_State[0][i];
		}

		delta_H_Average= H_old - H_new;

		return 1;
		
	}
	delta_H_Average = H_old - H_new;

	return 0;
	
}
