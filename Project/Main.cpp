/*
	Author: Aneirin John Baker 
	Date: 23/09/17
	Description: Simulation of a Quantum Harmonic/Anharmonic Oscillator using Hamiltonian Monte Carlo 
		techniques to advance the simulation for my Masters Project.This is the main class where all of 
		the main functions will be called from and the final results collected and outputted from. Here
		I shall be using natural units such that h_bar = c = 1. Beginning with a 1d system. 
*/
#include "Main.h"

int main(){
	printf("\n");
	printf("Beginning Simulation Initalising System\n\n");
	clock_t t1,t2;

	//Import the data from the file which holds all the data about the simulation 
	//in form 
	//iterations length t_step mu lamba 

	vector<double> input;
	ifstream inputFile;
	double number=0;
	inputFile.open("input.txt");
	if(!inputFile)
	{
		printf("Error:The file Can't be Openend Aborting Simulation");
	}
	else
	{
		while(!inputFile.eof())
		{
			inputFile >> number;
			input.push_back(number);
		}
	}

	//Number of iterations of the HMC algorithm to be performed, and number of times the algoirthm is going to loop
	unsigned int iterations = (unsigned int)input[0],length = (unsigned int)input[1];

	double t_step=input[2],mu=input[3],lamba=input[4];

	printf("##########Simulation Parameters##########\n");
	printf("\n");
	printf("Oscillators:	%d\nIterations:	%d\nT_Step:         %f\n",length,iterations,t_step);
	printf("\n");
	printf("##########Equation Parameters############\n");
	printf("\n");
	printf("Lattice Spacing:	%d\nMass:   		%d\nmu^2:			%f\nLamba:  		%f\n",1,1,mu,lamba);

	vector<double> v2(length,0);
	vector<vector<double> >lattice(iterations,v2);

	t1=clock();
	lattice_Evolution(lattice,length,t_step,iterations,mu,lamba);
	t2=clock();
	
	float seconds =((float)t2-(float)t1)/(CLOCKS_PER_SEC);
	printf("Simulation Completed in %f seconds\n\n",seconds);
	printf("########## End Simulation ##########\n\n");

}
