import math
from numpy import *
import numpy as np

mu2 = 1
a  = 1

omega = sqrt(mu2*(1+(a*a*mu2)/(4)))


R = 1 + (a*a*mu2)/(2) - a*omega

j=3
N=1000
N_j=N-j
norm =0
save =0

for j in range(0,10):
	avgij = 1/(2*omega) * (R**j + R**N_j)/(1-R**N)
	if (j==0):
		save = avgij
	norm=avgij/save
	print(str(j) + " " + str(avgij) + " " + str(norm))