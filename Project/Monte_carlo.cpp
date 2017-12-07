/*
	Author: Aneirin John Baker 
	Date: 23/09/17
	Description: The script is where the HMC algorithm will take place. The intergration methods will be 
	housed here and will be executed here and all of the stats functions will be executed here to create the Raw stats data.
*/
#include "Monte_carlo.h"
#define Oscillator_flip 1
//1 = Harmonic
//0 = Anharmonic
#define Anharmonic_flip 0
//1 = lamba and mu
//0 = lamba and f

void lattice_Evolution(vector<vector<double> > &lattice,unsigned int length,double t_step,unsigned int iterations,double mu,double lamba,double m,double a)
{
printf("##########################\n");
printf("\n");
#if Oscillator_flip
	printf("Running in Harmonic Mode\n\n");
#endif

#if !Oscillator_flip
	printf("Running in Anharmonic Mode\n");

	#if Anharmonic_flip
	printf("Running with Mu and Lamba\n\n");
	#endif

	#if !Anharmonic_flip
	printf("Running with Lamba and f\n\n");
	#endif

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
 	unsigned int steps =15,burn=0;


//initalise the first state of the siulation 
 	for(unsigned int j=0;j<length;j++)
 	{
 	//first_state[1][j]=State[1][j];
 	State[1][j] = Udistribution(generator);
 	if(j % 2== 0)
 	{
 		State[1][j]= State[1][j] * -1;
 	}
 	//State[1][j] = 1;
	}

 	//run main algorithm
 	for(unsigned int i = 0; i<iterations;i++)
 	{
 		default_random_engine generator(random_device{}());
 		for(unsigned int j = 0; j<length;j++)
 		{
 			State[0][j] = distribution(generator);
 			//State[0][j] = 1;

 		}
#if Oscillator_flip 
 		acceptance += hmcAlgorithm_Harmonic(length,t_step,State,temp_State,H_store,mu,steps,delta_H_Average,m,a);
#endif

#if !Oscillator_flip
 		acceptance += hmcAlgorithm_Anharmonic(length,t_step,State,temp_State,H_store,mu,lamba,steps,delta_H_Average,m,a);

#endif
//perform the stats calculations for the raw data
 	// 	for(unsigned int k = 0;k<length;k++)
 	// 	{
		// 	square_state[k] = State[1][k] * State[1][k];
		// }
 		//avgx = avgX(square_state);
 		//avgx2 = avg_X_Sqd(square_state);
		//error_x = standard_Deviation(avgx2,avgx,length);
		temp_avgx = avgX(State[1]);
		temp_avgx2 = avg_X_Sqd(State[1]);
		temp_avgx4 = avg_X_four(State[1]);
		dH_avg += delta_H_Average;
		if(i>burn)
		{
 		avgx +=temp_avgx;
 		avgx2 +=temp_avgx2;
 		avgx4 += temp_avgx4;
 		//error_x2 = standard_Deviation(avgx2,avgx,length);
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

double hmcAlgorithm_Harmonic(unsigned int length,double t_step,vector<vector<double> > &old_state,vector<vector<double> > &temp_State,vector<double> &H_store,double mu,unsigned int steps,double &delta_H_Average,double m ,double a)
{

	double min=0,H_old=0,H_new=0,H_inter=0;

	H_old=lattice_Hamiltonian(old_state,length,mu,0,m,a);

	//half step in the p
	temp_State[0][0] = old_state[0][0] -  (0.5*t_step * ((a*mu*old_state[1][0]) - ((m/a)*(old_state[1][1]+old_state[1][length-1]-(2*old_state[1][0])))));
	temp_State[1][0] = old_state[1][0];
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = old_state[0][j] - (0.5*t_step * ((a*mu*old_state[1][j]) - ((m/a)*(old_state[1][j+1]+old_state[1][j-1]-(2*old_state[1][j])))));
		temp_State[1][j] = old_state[1][j];
	}
	temp_State[0][length-1] = old_state[0][length-1] - (0.5*t_step * ((a*mu*old_state[1][length-1]) - ((m/a)*(old_state[1][0]+old_state[1][length-2]-(2*old_state[1][length-1])))));
	temp_State[1][length-1] = old_state[1][length-1];

	//full step in p and q for n steps
	for(unsigned int i = 0;i<steps;i++)
	{
		//update all q's
		for(unsigned int j = 0;j<length;j++)
		{
			temp_State[1][j] = temp_State[1][j] + ((t_step/m) * temp_State[0][j]);
		}

#if 1
//a full step for when running the algorithm normally
		if(i != steps-1)
		{
			temp_State[0][0] = temp_State[0][0] -  (t_step * ((a*mu*temp_State[1][0]) - ((m/a)*(temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));

			for(unsigned int j = 1;j<length-1;j++)
			{
				temp_State[0][j] = temp_State[0][j] -  (t_step * ((a*mu*temp_State[1][j]) - ((m/a)*(temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
			}

			temp_State[0][length-1] = temp_State[0][length-1] - (t_step * ((a*mu*temp_State[1][length-1]) - ((m/a)*(temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1])))));
		}
#endif 

#if 0
//two half steps for running when checking the algorithm still is equivalent to one full step only with a calcuation of the hamiltonian in the middle
		temp_State[0][0] = temp_State[0][0] -  (0.5 * t_step * (temp_State[0][1] - ((temp_State[1][1]+temp_State[length-1][1]-(2*temp_State[0][1])))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[j][0] = temp_State[j][0] -  (0.5 * t_step * (temp_State[j][1] - ((temp_State[j+1][1]+temp_State[j-1][1]-(2*temp_State[j][1])))));
		}

		temp_State[length-1][0] = temp_State[length-1][0] - (0.5 * t_step * (temp_State[length-1][1] - ((temp_State[0][1]+temp_State[length-2][1]-(2*temp_State[length-1][1])))));
		//calcuate the hamiltonian here.
		H_store[i+1] +=lattice_Hamiltonian(temp_State,length)-H_old;


		temp_State[0][0] = temp_State[0][0] -  (0.5 * t_step * (temp_State[0][1] - ((temp_State[1][1]+temp_State[length-1][1]-(2*temp_State[0][1])))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[j][0] = temp_State[j][0] -  (0.5 * t_step * (temp_State[j][1] - ((temp_State[j+1][1]+temp_State[j-1][1]-(2*temp_State[j][1])))));
		}

		temp_State[length-1][0] = temp_State[length-1][0] - (0.5 * t_step * (temp_State[length-1][1] - ((temp_State[0][1]+temp_State[length-2][1]-(2*temp_State[length-1][1])))));

#endif
	}
	//half step in the p
	temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * ((a*mu*temp_State[1][0]) - ((m/a)*(temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = temp_State[0][j] - (0.5*t_step * ((a*mu*temp_State[1][j]) - ((m/a)*(temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
	}
	temp_State[0][length-1] = temp_State[0][length-1] - (0.5*t_step * ((a*mu*temp_State[1][length-1]) - ((m/a)*(temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1])))));

	H_new = lattice_Hamiltonian(temp_State,length,mu,0,m,a);

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
		//acceptance=acceptance+1;
		delta_H_Average= H_old - H_new;
		//printf("%f %f %f\n",r,exp(H_old - H_new),H_old - H_new);
		return 1;
		
	}
	delta_H_Average = H_old - H_new;
//	printf("%f %f %f\n",r,exp(H_old - H_new),H_old - H_new);
	return 0;
	
}
/*
double hmcAlgorithm_Anharmonic(unsigned int length,double t_step,vector<vector<double> > &old_state,vector<vector<double> > &temp_State,vector<double> &H_store,double mu,double lamba,unsigned int steps,double &delta_H_Average,double m,double a)
{
	//printf("welp\n");
	lamba = lamba * 4;
	double f=0,H_old=0,H_new=0,H_inter=0,min=0;
	H_old=lattice_Hamiltonian(old_state,length,mu,lamba/4,m,a);
	
	//half step in the p

//old potential ie q^4
#if Anharmonic_flip
	temp_State[0][0] = old_state[0][0] -  (0.5*t_step * (a*lamba*pow(old_state[1][0],3) + (a*mu* old_state[1][0]) - ((m/a)*(old_state[1][1]+old_state[1][length-1]-(2*old_state[1][0])))));
	temp_State[1][0] = old_state[1][0];
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = old_state[0][j] - (0.5*t_step * (a*lamba*pow(old_state[1][j],3) + (a*mu*old_state[1][j]) - ((m/a)*(old_state[1][j+1]+old_state[1][j-1]-(2*old_state[1][j])))));
		temp_State[1][j] = old_state[1][j];
	}
	temp_State[0][length-1] = old_state[0][length-1] - (0.5*t_step * (a*lamba*pow(old_state[1][length-1],3) + (a*mu *old_state[1][length-1]) - ((m/a)*(old_state[1][0]+old_state[1][length-2]-(2*old_state[1][length-1])))));
	temp_State[1][length-1] = old_state[1][length-1];
#endif
	//modified potential V = (q^2 - f^2)^2
#if !Anharmonic_flip
	temp_State[0][0] = old_state[0][0] -  (0.5*t_step * ((a*lamba*old_state[1][0]*(pow(old_state[1][0],2)-f)) - ((m/a)*(old_state[1][1]+old_state[1][length-1]-(2*old_state[1][0])))));
	temp_State[1][0] = old_state[1][0];
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = old_state[0][j] - (0.5*t_step * ((a*lamba*old_state[1][j]*(pow(old_state[1][j],2)-f))- ((m/a)*(old_state[1][j+1]+old_state[1][j-1]-(2*old_state[1][j])))));
		temp_State[1][j] = old_state[1][j];
		//printf("%f %f\n",temp_State[j][0],temp_State[j][1]);
	}
	temp_State[0][length-1] = old_state[0][length-1] - (0.5*t_step * ((a*lamba*old_state[1][length-1]*(pow(old_state[1][length-1],2)-f))- ((m/a)*(old_state[1][0]+old_state[1][length-2]-(2*old_state[1][length-1])))));
	temp_State[1][length-1] = old_state[1][length-1];
#endif

	//full step in p and q for n steps
	for(unsigned int i = 0;i<steps;i++)
	{
		//update all q's

		for(unsigned int j = 0;j<length;j++)
		{
			temp_State[1][j] = temp_State[1][j] + ((t_step/m) *temp_State[0][j]);
		}

#if 1
//a full step for when running the algorithm normally

		if(i != steps-1)
		{
		
		//potential as V=(x^2 - f^2)^2

#if Anharmonic_flip		
		temp_State[0][0] = temp_State[0][0] -  (t_step * (a*lamba*pow(temp_State[1][0],3) + (a*mu*temp_State[1][0]) - ((m/a)*(temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[0][j] = temp_State[0][j] -  (t_step * (a*lamba*pow(temp_State[1][j],3) + (a*mu*temp_State[1][j]) - ((m/a)*(temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
		}

		temp_State[0][length-1] = temp_State[0][length-1] - (t_step *(a*lamba*pow(temp_State[1][length-1],3) + (a*mu*temp_State[1][length-1]) - ((m/a)*(temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1])))));
#endif

#if !Anharmonic_flip
		temp_State[0][0] = temp_State[0][0] -  (t_step * ((a*lamba*temp_State[1][0]*(pow(temp_State[1][0],2)-f)) - ((m/a)*(temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[0][j] = temp_State[0][j] -  (t_step * ((a*lamba*temp_State[1][j]*(pow(temp_State[1][j],2)-f))  - ((m/a)*(temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
		}

		temp_State[0][length-1] = temp_State[0][length-1] - (t_step * ((a*lamba*temp_State[1][length-1]*(pow(temp_State[1][length-1],2)-f))  - ((m/a)*(temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1])))));
#endif
		}
#endif 

#if 0
//two half steps for running when checking the algorithm still is equivalent to one full step only with a calcuation of the hamiltonian in the middle
		temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * (pow(temp_State[1][0],3) + temp_State[1][0] - ((temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[0][j] = temp_State[0][j] -  (0.5*t_step * (pow(temp_State[1][j],3) + temp_State[1][j] - ((temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
		}

		temp_State[0][length-1] = temp_State[0][length-1] - (0.5*t_step *(pow(temp_State[1][length-1],3) + (temp_State[1][length-1] - ((temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1]))))));

		//calcuate the hamiltonian here.
		H_store[i+1] +=lattice_Hamiltonian(temp_State,length)-H_old;


		temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * (pow(temp_State[1][0],3) + temp_State[1][0] - ((temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[0][j] = temp_State[0][j] -  (0.5*t_step * (pow(temp_State[1][j],3) + temp_State[1][j] - ((temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
		}

		temp_State[0][length-1] = temp_State[0][length-1] - (0.5 * t_step *(pow(temp_State[1][length-1],3) + t_step * (temp_State[1][length-1] - ((temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1]))))));

#endif
	}
	//half step in the p
#if Anharmonic_flip
	temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * (a*lamba*pow(temp_State[1][0],3) + (a*mu *temp_State[1][0]) - ((m/a)*(temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = temp_State[0][j] - (0.5*t_step * (a*lamba*pow(temp_State[1][j],3) + (a*mu *temp_State[1][j]) - ((m/a)*(temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
	}

	temp_State[0][length-1] = temp_State[0][length-1] - (0.5*t_step * (a*lamba*pow(temp_State[1][length-1],3) + (a*mu*temp_State[1][length-1]) - ((m/a)*(temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1])))));
#endif
//Potential at V=(x^2-f^2)^2
#if !Anharmonic_flip
	temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * ((a*lamba*temp_State[1][0]*(pow(temp_State[1][0],2)-f)) - ((m/a)*(temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = temp_State[0][j] - (0.5*t_step * ((a*lamba*temp_State[1][j]*(pow(temp_State[1][j],2)-f)) - ((m/a)*(temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
	}
	temp_State[0][length-1] = temp_State[0][length-1] - (0.5*t_step * ((a*lamba*temp_State[1][length-1]*(pow(temp_State[1][length-1],2)-f)) - ((m/a)*(temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1])))));
#endif
	H_new = lattice_Hamiltonian(temp_State,length,mu,lamba/4,m,a);

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
		//acceptance = acceptance +1;
	//	printf("welp1\n");
		delta_H_Average= H_old - H_new;
		return 1;

	}
	//printf("welp\n");

	delta_H_Average= H_old - H_new;
	return 0;
	
}
*/
double hmcAlgorithm_Anharmonic(unsigned int length,double t_step,vector<vector<double> > &old_state,vector<vector<double> > &temp_State,vector<double> &H_store,double mu,double lamba,unsigned int steps,double &delta_H_Average,double m,double a)
{
	//printf("welp\n");
	double f=2,H_old=0,H_new=0,H_inter=0,min=0;
	H_old=lattice_Hamiltonian(old_state,length,mu,lamba,m,a);
	
	//half step in the p
	temp_State[0][0] = old_state[0][0] -  (0.5*t_step * ((4*a*lamba*old_state[1][0]*(pow(old_state[1][0],2)-f)) - ((m/a)*(old_state[1][1]+old_state[1][length-1]-(2*old_state[1][0])))));
	temp_State[1][0] = old_state[1][0];
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = old_state[0][j] - (0.5*t_step * ((4*a*lamba*old_state[1][j]*(pow(old_state[1][j],2)-f)) - ((m/a)*(old_state[1][j+1]+old_state[1][j-1]-(2*old_state[1][j])))));
		temp_State[1][j] = old_state[1][j];
	}
	temp_State[0][length-1] = old_state[0][length-1] - (0.5*t_step * ((4*a*lamba*old_state[1][length-1]*(pow(old_state[1][length-1],2)-f)) - ((m/a)*(old_state[1][0]+old_state[1][length-2]-(2*old_state[1][length-1])))));
	temp_State[1][length-1] = old_state[1][length-1];

	//full step in p and q for n steps
	for(unsigned int i = 0;i<steps;i++)
	{
		//update all q's
		for(unsigned int j = 0;j<length;j++)
		{
			temp_State[1][j] = temp_State[1][j] + ((t_step/m) * temp_State[0][j]);
		}

#if 1
//a full step for when running the algorithm normally
		if(i != steps-1)
		{
			temp_State[0][0] = temp_State[0][0] -  (t_step * ((4*a*lamba*temp_State[1][0]*(pow(temp_State[1][0],2)-f)) - ((m/a)*(temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));

			for(unsigned int j = 1;j<length-1;j++)
			{
				temp_State[0][j] = temp_State[0][j] -  (t_step * ((4*a*lamba*temp_State[1][j]*(pow(temp_State[1][j],2)-f)) - ((m/a)*(temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
			}

			temp_State[0][length-1] = temp_State[0][length-1] - (t_step * ((4*a*lamba*temp_State[1][length-1]*(pow(temp_State[1][length-1],2)-f)) - ((m/a)*(temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1])))));
		}
#endif 

	}
	//half step in the p
	temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * ((4*a*lamba*temp_State[1][0]*(pow(temp_State[1][0],2)-f)) - ((m/a)*(temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = temp_State[0][j] - (0.5*t_step * ((4*a*lamba*temp_State[1][j]*(pow(temp_State[1][j],2)-f)) - ((m/a)*(temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
	}
	temp_State[0][length-1] = temp_State[0][length-1] - (0.5*t_step * ((4*a*lamba*temp_State[1][length-1]*(pow(temp_State[1][length-1],2)-f)) - ((m/a)*(temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1])))));

	H_new = lattice_Hamiltonian(temp_State,length,mu,lamba,m,a);

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
	delta_H_Average= H_old - H_new;
	return 0;
	
}

double hmcAlgorithm_Anharmonic_uncoupled(unsigned int length,double t_step,vector<vector<double> > &old_state,vector<vector<double> > &temp_State,vector<double> &H_store,double mu,double lamba,unsigned int steps,double m,double a)
{
	mu = mu * 2;
	lamba = lamba * 4;
	double f=2;
	double H_old=0,H_new=0,H_inter=0,min=0;

	H_old=lattice_Hamiltonian(old_state,length,mu,lamba,m,a);
	
	//half step in the p

//old potential ie q^4
	/*
	temp_State[0][0] = old_state[0][0] -  (0.5*t_step * (lamba*pow(old_state[1][0],3) + (mu* old_state[1][0]) - (old_state[1][1]+old_state[1][length-1]-(2*old_state[1][0]))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = old_state[0][j] - (0.5*t_step * (lamba*pow(old_state[1][j],3) + (mu*old_state[1][j]) - (old_state[1][j+1]+old_state[1][j-1]-(2*old_state[1][j]))));
		temp_State[1][j] = old_state[1][j];
	}
	temp_State[0][length-1] = old_state[0][length-1] - (0.5*t_step * (lamba*pow(old_state[1][length-1],3) + (mu *old_state[1][length-1]) - (old_state[1][0]+old_state[1][length-2]-(2*old_state[1][length-1]))));
*/
	//modified potential V = (q^2 - f^2)^2
	
	temp_State[0][0] = old_state[0][0] -  (0.5*t_step * ((4*old_state[1][0]*(pow(old_state[1][0],2)-f))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = old_state[0][j] - (0.5*t_step * ((4*old_state[1][j]*(pow(old_state[1][j],2)-f))));
		temp_State[1][j] = old_state[1][j];
		//printf("%f %f\n",temp_State[j][0],temp_State[j][1]);
	}
	temp_State[0][length-1] = old_state[0][length-1] - (0.5*t_step * ((4*old_state[1][length-1]*(pow(old_state[1][length-1],2)-f))));

	//full step in p and q for n steps
	for(unsigned int i = 0;i<steps;i++)
	{
		//update all q's
		for(unsigned int j = 0;j<length;j++)
		{
			temp_State[1][j] = temp_State[1][j] + (t_step * temp_State[0][j]);
		}

#if 1
//a full step for when running the algorithm normally

		if(i != steps-1)
		{
		
		//potential as V=(x^2 - f^2)^2
		temp_State[0][0] = temp_State[0][0] -  (t_step * ((4*temp_State[1][0]*(pow(temp_State[1][0],2)-f))));


		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[0][j] = temp_State[0][j] -  (t_step * ((4*temp_State[1][j]*(pow(temp_State[1][j],2)-f))));
		}


		temp_State[0][length-1] = temp_State[0][length-1] - (t_step * ((4*temp_State[1][length-1]*(pow(temp_State[1][length-1],2)-f))));
		/*
		temp_State[0][0] = temp_State[0][0] -  (t_step * (lamba*pow(temp_State[1][0],3) + (mu*temp_State[1][0]) - ((temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[0][j] = temp_State[0][j] -  (t_step * (lamba*pow(temp_State[1][j],3) + (mu*temp_State[1][j]) - ((temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
		}

		temp_State[0][length-1] = temp_State[0][length-1] - (t_step *(lamba*pow(temp_State[1][length-1],3) + (mu*temp_State[1][length-1]) - ((temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1])))));
		*/
		}
#endif 

#if 0
//two half steps for running when checking the algorithm still is equivalent to one full step only with a calcuation of the hamiltonian in the middle
		temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * (pow(temp_State[1][0],3) + temp_State[1][0] - ((temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[0][j] = temp_State[0][j] -  (0.5*t_step * (pow(temp_State[1][j],3) + temp_State[1][j] - ((temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
		}

		temp_State[0][length-1] = temp_State[0][length-1] - (0.5*t_step *(pow(temp_State[1][length-1],3) + (temp_State[1][length-1] - ((temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1]))))));

		//calcuate the hamiltonian here.
		H_store[i+1] +=lattice_Hamiltonian(temp_State,length)-H_old;


		temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * (pow(temp_State[1][0],3) + temp_State[1][0] - ((temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[0][j] = temp_State[0][j] -  (0.5*t_step * (pow(temp_State[1][j],3) + temp_State[1][j] - ((temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
		}

		temp_State[0][length-1] = temp_State[0][length-1] - (0.5 * t_step *(pow(temp_State[1][length-1],3) + t_step * (temp_State[1][length-1] - ((temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1]))))));

#endif
	}
	//half step in the p
	/*
	temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * (lamba*pow(temp_State[1][0],3) + (mu *temp_State[1][0]) - (temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0]))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = temp_State[0][j] - (0.5*t_step * (lamba*pow(temp_State[1][j],3) + (mu *temp_State[1][j]) - (temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j]))));
		temp_State[1][j] = temp_State[1][j];
		//printf("%f %f\n",temp_State[j][0],temp_State[j][1]);
	}
	temp_State[0][length-1] = temp_State[0][length-1] - (0.5*t_step * (lamba*pow(temp_State[1][length-1],3) + (mu*temp_State[1][length-1]) - (temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1]))));
	*/
//Potential at V=(x^2-f^2)^2
	temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * ((4*temp_State[1][0]*(pow(temp_State[1][0],2)-f))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = temp_State[0][j] - (0.5*t_step * ((4*temp_State[1][j]*(pow(temp_State[1][j],2)-f))));
	}
	temp_State[0][length-1] = temp_State[0][length-1] - (0.5*t_step * ((4*temp_State[1][length-1]*(pow(temp_State[1][length-1],2)-f))));
	
	H_new = lattice_Hamiltonian(temp_State,length,mu,lamba,m,a);

	//metroplis update
	double r = ((double) rand() / (RAND_MAX));

	min = (1 < exp(H_old - H_new)) ? 1 : exp(H_old - H_new);
	if(r < min)
	{
		//accept
		for(unsigned int i = 0;i<length;i++)
		{
			old_state[1][i] = temp_State[1][i];

		}
		//return H_old - H_new;
		return 1;
	}
	//return H_old - H_new;
	return 0;
}
