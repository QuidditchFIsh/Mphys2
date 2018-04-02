from numpy import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
import statistics 
from array import *
import math
from matplotlib import gridspec



#initalise arrays and variables
Mavgx2=[];Mavgx21=[]
avgx2=[];avgx21=[]
sum2=0;sum21=0




file  = open("/home/s1403955/Mphys2/Project/HMC_Stats.dat",'r')
file1  = open("HMC_Stats.dat",'r')

for line in file:
	a,b,c,d,e,f,gg = line.split(' ', 6)
	avgx2.append(float(d))

for line in file1:
	a1,b1,c1,d1,e1,f1,gg1 = line.split(' ', 6)
	avgx21.append(float(d1))


length = len(avgx2)
variation = [];variation1 = []

for i in range(1,length):
	sum2 += avgx2[i-1]
	Mavgx2.append(sum2/i)
	variation.append((Mavgx2[i-1] - 0.4472135955)**2)

	sum21 += avgx21[i-1]
	Mavgx21.append(sum21/i)
	variation1.append((Mavgx21[i-1] - 0.4472135955)**2)

file = open("variation.dat","w")

for i in range(0,length-1):
	file.write(str(i) + " " +str(variation[i]) + " " + str(variation1[i]) + "\n")
file.close()

print("Fourier Acceleration")
print(min(variation))
print(max(variation))

print("Regular HMC")
print(min(variation1))
print(max(variation1))

print("Sum of variation under curve")
FAsum  = np.trapz(variation)
HMCsum = np.trapz(variation1)
print("Fourier Accelerated:" + str(FAsum))
print("Regular HMC:" + str(HMCsum))

maxX = 1000

FA = plt.figure(figsize=(8,6))
plt.ylim(-0.0001,max(variation))
plt.xlim(0,maxX)
plt.plot(variation)
FA.savefig("graphs/Fourier HMC Variation <x^2>.pdf")

HMC = plt.figure(figsize=(8,6))
plt.ylim(-0.0001,max(variation1))
plt.xlim(0,maxX)
plt.plot(variation1)
HMC.savefig("graphs/HMC Variation <x^2>.pdf")

COMP = plt.figure(figsize=(8,6))
plt.ylim(-0.0001,max(variation1))
plt.xlim(0,maxX)
plt.plot(variation1,label = "No Fourier Acceleration")
plt.plot(variation, label = "Fourier Accelerated")
plt.legend(loc="upper right")
plt.title("Squared Difference Comparison of Fourier Acclerated HMC to HMC")
plt.xlabel("Monte Carlo Iterations")
plt.ylabel("(<x^2> - 0.44721)^2")
COMP.savefig("graphs/Comparison.pdf")








