double hmcAlgorithm_Anharmonic(unsigned int length,double t_step,vector<vector<double> > &old_state,vector<vector<double> > &temp_State,vector<double> &H_store,double mu,double lamba,unsigned int steps)
{
	mu = mu * 2;
	lamba = lamba * 4;
	double f=-1,H_old=0,H_new=0,H_inter=0,min=0;

	H_old=lattice_Hamiltonian(old_state,length);
	
	//half step in the p

//old potential ie q^4
#if anharmonic_flip
	temp_State[0][0] = old_state[0][0] -  (0.5*t_step * (lamba*pow(old_state[1][0],3) + (mu* old_state[1][0]) - (old_state[1][1]+old_state[1][length-1]-(2*old_state[1][0]))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = old_state[0][j] - (0.5*t_step * (lamba*pow(old_state[1][j],3) + (mu*old_state[1][j]) - (old_state[1][j+1]+old_state[1][j-1]-(2*old_state[1][j]))));
		temp_State[1][j] = old_state[1][j];
	}
	temp_State[0][length-1] = old_state[0][length-1] - (0.5*t_step * (lamba*pow(old_state[1][length-1],3) + (mu *old_state[1][length-1]) - (old_state[1][0]+old_state[1][length-2]-(2*old_state[1][length-1]))));
#endif
	//modified potential V = (q^2 - f^2)^2
#if !anharmonic_flip
	temp_State[0][0] = old_state[0][0] -  (0.5*t_step * ((lamba*old_state[1][0]*(pow(old_state[1][0],2)-f)) - (old_state[1][1]+old_state[1][length-1]-(2*old_state[1][0]))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = old_state[0][j] - (0.5*t_step * ((lamba*old_state[1][j]*(pow(old_state[1][j],2)-f))- (old_state[1][j+1]+old_state[1][j-1]-(2*old_state[1][j]))));
		temp_State[1][j] = old_state[1][j];
		//printf("%f %f\n",temp_State[j][0],temp_State[j][1]);
	}
	temp_State[0][length-1] = old_state[0][length-1] - (0.5*t_step * ((lamba*old_state[1][length-1]*(pow(old_state[1][length-1],2)-f))- (old_state[1][0]+old_state[1][length-2]-(2*old_state[1][length-1]))));
#endif
	//full step in p and q for n steps
	for(unsigned int i = 0;i<steps;i++)
	{
		//update all q's



		for(unsigned int j = 0;j<length;j++)
		{
			temp_State[1][j] = temp_State[1][j] + (t_step *temp_State[0][j]);
		}

#if 1
//a full step for when running the algorithm normally

		if(i != steps-1)
		{
		
		//potential as V=(x^2 - f^2)^2
#if !anharmonic_flip
		temp_State[0][0] = temp_State[0][0] -  (t_step * ((lamba*temp_State[1][0]*(pow(temp_State[1][0],2)-f)) - ((temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[0][j] = temp_State[0][j] -  (t_step * ((lamba*temp_State[1][j]*(pow(temp_State[1][j],2)-f))  - ((temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
		}

		temp_State[0][length-1] = temp_State[0][length-1] - (t_step * ((lamba*temp_State[1][length-1]*(pow(temp_State[1][length-1],2)-f))  - ((temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1])))));
#endif
#if anharmonic_flip		
		temp_State[0][0] = temp_State[0][0] -  (t_step * (lamba*pow(temp_State[1][0],3) + (mu*temp_State[1][0]) - ((temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[0][j] = temp_State[0][j] -  (t_step * (lamba*pow(temp_State[1][j],3) + (mu*temp_State[1][j]) - ((temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
		}

		temp_State[0][length-1] = temp_State[0][length-1] - (t_step *(lamba*pow(temp_State[1][length-1],3) + (mu*temp_State[1][length-1]) - ((temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1])))));
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
#if anharmonic_flip
	temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * (lamba*pow(temp_State[1][0],3) + (mu *temp_State[1][0]) - (temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0]))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = temp_State[0][j] - (0.5*t_step * (lamba*pow(temp_State[1][j],3) + (mu *temp_State[1][j]) - (temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j]))));
		temp_State[1][j] = temp_State[1][j];
		//printf("%f %f\n",temp_State[j][0],temp_State[j][1]);
	}
	temp_State[0][length-1] = temp_State[0][length-1] - (0.5*t_step * (lamba*pow(temp_State[1][length-1],3) + (mu*temp_State[1][length-1]) - (temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1]))));
#endif
//Potential at V=(x^2-f^2)^2
#if !anharmonic_flip
	temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * ((lamba*temp_State[1][0]*(pow(temp_State[1][0],2)-f)) - (temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0]))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = temp_State[0][j] - (0.5*t_step * ((lamba*temp_State[1][j]*(pow(temp_State[1][j],2)-f)) - (temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j]))));
	}
	temp_State[0][length-1] = temp_State[0][length-1] - (0.5*t_step * ((lamba*temp_State[1][length-1]*(pow(temp_State[1][length-1],2)-f)) - (temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1]))));
#endif
	H_new = lattice_Hamiltonian(temp_State,length);

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






double hmcAlgorithm_Anharmonic(unsigned int length,double t_step,vector<vector<double> > &old_state,vector<vector<double> > &temp_State,vector<double> &H_store,double mu,double lamba,unsigned int steps)
{
	lamba = lamba * 4;
	double f=-1,H_old=0,H_new=0,H_inter=0,min=0;

	H_old=lattice_Hamiltonian(old_state,length);
	
	//half step in the p
	temp_State[0][0] = old_state[0][0] -  (0.5*t_step * (lamba*pow(old_state[1][0],3) + (mu* old_state[1][0]) - (old_state[1][1]+old_state[1][length-1]-(2*old_state[1][0]))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = old_state[0][j] - (0.5*t_step * (lamba*pow(old_state[1][j],3) + (mu*old_state[1][j]) - (old_state[1][j+1]+old_state[1][j-1]-(2*old_state[1][j]))));
		temp_State[1][j] = old_state[1][j];
	}
	temp_State[0][length-1] = old_state[0][length-1] - (0.5*t_step * (lamba*pow(old_state[1][length-1],3) + (mu *old_state[1][length-1]) - (old_state[1][0]+old_state[1][length-2]-(2*old_state[1][length-1]))));


	//full step in p and q for n steps
	for(unsigned int i = 0;i<steps;i++)
	{
		//update all q's
		for(unsigned int j = 0;j<length;j++)
		{
			temp_State[1][j] = temp_State[1][j] + (t_step *temp_State[0][j]);
		}


//a full step for when running the algorithm normally

		if(i != steps-1)
		{
		
		temp_State[0][0] = temp_State[0][0] -  (t_step * (lamba*pow(temp_State[1][0],3) + (mu*temp_State[1][0]) - ((temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[0][j] = temp_State[0][j] -  (t_step * (lamba*pow(temp_State[1][j],3) + (mu*temp_State[1][j]) - ((temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
		}

		temp_State[0][length-1] = temp_State[0][length-1] - (t_step *(lamba*pow(temp_State[1][length-1],3) + (mu*temp_State[1][length-1]) - ((temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1])))));

		}

	}
	//half step in the p

	temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * (lamba*pow(temp_State[1][0],3) + (mu *temp_State[1][0]) - (temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0]))));

	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = temp_State[0][j] - (0.5*t_step * (lamba*pow(temp_State[1][j],3) + (mu *temp_State[1][j]) - (temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j]))));
		temp_State[1][j] = temp_State[1][j];
	}

	temp_State[0][length-1] = temp_State[0][length-1] - (0.5*t_step * (lamba*pow(temp_State[1][length-1],3) + (mu*temp_State[1][length-1]) - (temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1]))));

	H_new = lattice_Hamiltonian(temp_State,length);

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



//harmonic with the lattice size in not sure if works 

double hmcAlgorithm_Harmonic(unsigned int length,double t_step,vector<vector<double> > &old_state,vector<vector<double> > &temp_State,vector<double> &H_store,double mu,unsigned int steps)
{

	double min=0,H_old=0,H_new=0,H_inter=0,a=1;

	H_old=lattice_Hamiltonian(old_state,length);

	//half step in the p
	temp_State[0][0] = old_state[0][0] -  (0.5*t_step * a*((mu*old_state[1][0]) - ((old_state[1][1]+old_state[1][length-1]-(2*old_state[1][0]))*(1/a*a))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = old_state[0][j] - (0.5*t_step * a*((mu*old_state[1][j]) - ((old_state[1][j+1]+old_state[1][j-1]-(2*old_state[1][j]))*(1/a*a))));
		temp_State[1][j] = old_state[1][j];
	}
	temp_State[0][length-1] = old_state[0][length-1] - (0.5*t_step * a*((mu*old_state[1][length-1]) - ((old_state[1][0]+old_state[1][length-2]-(2*old_state[1][length-1]))*(1/a*a))));
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
		temp_State[0][0] = temp_State[0][0] -  (t_step * a*((mu*temp_State[1][0]) - ((temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0]))*(1/a*a))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[0][j] = temp_State[0][j] -  (t_step * a*((mu*temp_State[1][j]) - ((temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j]))*(1/a*a))));
		}

		temp_State[0][length-1] = temp_State[0][length-1] - (t_step * a*((mu*temp_State[1][length-1]) - ((temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1]))*(1/a*a))));
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
	temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * a*((mu*temp_State[1][0]) - ((1/a*a)*(temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = temp_State[0][j] - (0.5*t_step * a*((mu*temp_State[1][j]) - ((temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j]))*(1/a*a))));
	}
	temp_State[0][length-1] = temp_State[0][length-1] - (0.5*t_step * a*((mu*temp_State[1][length-1]) - ((temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1]))*(1/a*a))));

	H_new = lattice_Hamiltonian(temp_State,length);

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
		return H_old - H_new;
		//return 1;
	}
	return H_old - H_new;
	//return 0;
}



//harmonic algorithm without lattics size DEFINATELY works

double hmcAlgorithm_Harmonic(unsigned int length,double t_step,vector<vector<double> > &old_state,vector<vector<double> > &temp_State,vector<double> &H_store,double mu,unsigned int steps)
{

	double min=0,H_old=0,H_new=0,H_inter=0,a=1;

	H_old=lattice_Hamiltonian(old_state,length);

	//half step in the p
	temp_State[0][0] = old_state[0][0] -  (0.5*t_step * ((mu*old_state[1][0]) - ((old_state[1][1]+old_state[1][length-1]-(2*old_state[1][0])))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = old_state[0][j] - (0.5*t_step * ((mu*old_state[1][j]) - ((old_state[1][j+1]+old_state[1][j-1]-(2*old_state[1][j])))));
		temp_State[1][j] = old_state[1][j];
	}
	temp_State[0][length-1] = old_state[0][length-1] - (0.5*t_step * ((mu*old_state[1][length-1]) - ((old_state[1][0]+old_state[1][length-2]-(2*old_state[1][length-1])))));
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
		temp_State[0][0] = temp_State[0][0] -  (t_step * ((mu*temp_State[1][0]) - ((temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[0][j] = temp_State[0][j] -  (t_step * ((mu*temp_State[1][j]) - ((temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
		}

		temp_State[0][length-1] = temp_State[0][length-1] - (t_step * ((mu*temp_State[1][length-1]) - ((temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1])))));
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
	temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * ((mu*temp_State[1][0]) - ((temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = temp_State[0][j] - (0.5*t_step * ((mu*temp_State[1][j]) - ((temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
	}
	temp_State[0][length-1] = temp_State[0][length-1] - (0.5*t_step * ((mu*temp_State[1][length-1]) - ((temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1])))));

	H_new = lattice_Hamiltonian(temp_State,length);

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
		return H_old - H_new;
		//return 1;
	}
	return H_old - H_new;
	//return 0;
}

double hmcAlgorithm_Anharmonic(unsigned int length,double t_step,vector<vector<double> > &old_state,vector<vector<double> > &temp_State,vector<double> &H_store,double mu,double lamba,unsigned int steps,double &delta_H_Average,double m,double a)
{
	lamba = lamba * 4;
	double f=0,H_old=0,H_new=0,H_inter=0,min=0;

	H_old=lattice_Hamiltonian(old_state,length,mu,lamba,m,a);
	
	//half step in the p

//old potential ie q^4
#if anharmonic_flip
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
#if !anharmonic_flip
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
#if !anharmonic_flip
		temp_State[0][0] = temp_State[0][0] -  (t_step * ((a*lamba*temp_State[1][0]*(pow(temp_State[1][0],2)-f)) - ((m/a)*(temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[0][j] = temp_State[0][j] -  (t_step * ((a*lamba*temp_State[1][j]*(pow(temp_State[1][j],2)-f))  - ((m/a)*(temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
		}

		temp_State[0][length-1] = temp_State[0][length-1] - (t_step * ((a*lamba*temp_State[1][length-1]*(pow(temp_State[1][length-1],2)-f))  - ((m/a)*(temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1])))));
#endif
#if anharmonic_flip		
		temp_State[0][0] = temp_State[0][0] -  (t_step * (a*lamba*pow(temp_State[1][0],3) + (a*mu*temp_State[1][0]) - ((m/a)*(temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));

		for(unsigned int j = 1;j<length-1;j++)
		{
			temp_State[0][j] = temp_State[0][j] -  (t_step * (a*lamba*pow(temp_State[1][j],3) + (a*mu*temp_State[1][j]) - ((m/a)*(temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
		}

		temp_State[0][length-1] = temp_State[0][length-1] - (t_step *(a*lamba*pow(temp_State[1][length-1],3) + (a*mu*temp_State[1][length-1]) - ((m/a)*(temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1])))));
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
#if anharmonic_flip
	temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * (a*lamba*pow(temp_State[1][0],3) + (a*mu *temp_State[1][0]) - ((m/a)*(temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = temp_State[0][j] - (0.5*t_step * (a*lamba*pow(temp_State[1][j],3) + (a*mu *temp_State[1][j]) - ((m/a)*(temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
	}

	temp_State[0][length-1] = temp_State[0][length-1] - (0.5*t_step * (a*lamba*pow(temp_State[1][length-1],3) + (a*mu*temp_State[1][length-1]) - ((m/a)*(temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1])))));
#endif
//Potential at V=(x^2-f^2)^2
#if !anharmonic_flip
	temp_State[0][0] = temp_State[0][0] -  (0.5*t_step * ((a*lamba*temp_State[1][0]*(pow(temp_State[1][0],2)-f)) - ((m/a)*(temp_State[1][1]+temp_State[1][length-1]-(2*temp_State[1][0])))));
	for(unsigned int j = 1;j<length-1;j++)
	{
		temp_State[0][j] = temp_State[0][j] - (0.5*t_step * ((a*lamba*temp_State[1][j]*(pow(temp_State[1][j],2)-f)) - ((m/a)*(temp_State[1][j+1]+temp_State[1][j-1]-(2*temp_State[1][j])))));
	}
	temp_State[0][length-1] = temp_State[0][length-1] - (0.5*t_step * ((a*lamba*temp_State[1][length-1]*(pow(temp_State[1][length-1],2)-f)) - ((m/a)*(temp_State[1][0]+temp_State[1][length-2]-(2*temp_State[1][length-1])))));
#endif
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
		//acceptance = acceptance +1;
	//	printf("welp1\n");
		delta_H_Average= H_old - H_new;
		return 1;

	}
	//printf("welp\n");

	delta_H_Average= H_old - H_new;
	return 0;
	
}

 		printf("Iteration:%d\n\n",i);
 		for(int j=0;j<length;j++)
 		{
 			printf("%f ",State[0][j]);
 		}
 		printf("\n");
 		for(int j=0;j<length;j++)
 		{
 			printf("%f ",State[1][j]);
 		}
 		printf("\n\n");

 		acceptance += hmcAlgorithm_Anharmonic(length,t_step,State,temp_State,H_store,mu,lamba,steps,delta_H_Average,m,a);

 		for(int j=0;j<length;j++)
 		{
 			printf("%f ",State[0][j]);
 		}
 		printf("\n");
 		for(int j=0;j<length;j++)
 		{
 			printf("%f ",State[1][j]);
 		}
 		printf("\n\n");

 		acceptance += hmcAlgorithm_Anharmonic(length,-1*t_step,State,temp_State,H_store,mu,lamba,steps,delta_H_Average,m,a);

 		for(int j=0;j<length;j++)
 		{
 			printf("%f ",State[0][j]);
 		}
 		printf("\n");
 		for(int j=0;j<length;j++)
 		{
 			printf("%f ",State[1][j]);
 		}
 		printf("\n\n");







		p_temp[j] = p[j] - (0.5 * t_step * q[j] * 0.5 * a * (mu + (4 * sin((PI2/length) * j * a * 0.5) * sin((PI2/length) * j * a * 0.5))));
		p_temp[j] = p[j] - (0.5 * (t_step / (a * a)) * q[j] * ((a * a * mu) + (4 * sin((PI2/length) * j * a * 0.5) * sin((PI2/length) * j * a * 0.5))));




		p_temp[j] = p_temp[j] - (t_step  * q_temp[j] * 0.5 * a * (mu + (4 * sin((PI2/length) * j * a * 0.5) * sin((PI2/length) * j * a * 0.5))));
		p_temp[j] = p_temp[j] - (t_step / (a * a)) * q_temp[j] * ((a * a * mu) + (4 * sin((PI2/length) * j * a * 0.5) * sin((PI2/length) * j * a * 0.5)));

		p_temp[j] = p_temp[j] - (0.5 * t_step * q_temp[j] * 0.5 * a *(mu + (4 * sin((PI2/length) * j * a * 0.5) * sin((PI2/length) * j * a * 0.5))));
		p_temp[j] = p_temp[j] - (0.5 * (t_step / (a * a)) * q_temp[j] * ((a * a * mu) + (4 * sin((PI2/length) * j * a * 0.5) * sin((PI2/length) * j * a * 0.5))));