#include <stdlib.h>
#include <stdio.h>
#include <random>
#include <iostream>
#include <string>
#include <ctime>
#include <vector>
#include <cmath>
#include <algorithm>
#include "stastics.h"
#include "functions.h"
#include "constants.h"

using namespace std;

void lattice_Evolution(vector<vector<double> > &,unsigned int ,double ,unsigned int ,double,double,double,double,double);

double hmcAlgorithm(unsigned int ,double ,vector<vector<double> > &,vector<vector<double> > & ,vector<double> &,double,unsigned int,double &,double,double,double);
