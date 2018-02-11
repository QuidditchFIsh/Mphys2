from numpy import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
import statistics 
from array import *
import math
from matplotlib import gridspec



#initalise arrays and variables
avgx=[];avgx2=[];action=[];KE=[];delta_h=[];i=[];avgx4=[];jjj=[]
Mavgx=[];Mavgx2=[];Maction=[];MKE=[];Mdelta_h=[];Mavgx4=[]
avgxerr=[];avgx2err=[];actionerr=[];KEerr=[];delta_herr=[];avgx4err=[]
sum1=0;sum2=0;sum3=0;sum4=0;sum5=0;sum6=0
sum12=0;sum22=0;sum32=0;sum42=0;sum52=0;sum62=0

#import all data from the file
file  = open("HMC_Stats.dat",'r')
#file1 = open("HMC_X.dat","r")

for line in file:
	a,b,c,d,e,f,gg = line.split(' ', 6)
	i.append(float(a))
	avgx.append(float(b))
	#avgxerr.append(float(c))
	avgx2.append(float(d))
	avgx2err.append(float(e))
	action.append(float(f))
	#KE.append(float(gg))
	delta_h.append(float(c))
	#avgx4.append(float(c))

#dataX=np.genfromtxt("HMC_Final_X.dat", unpack=True)
#print(data)

#stats calculations
rm=0
for j in range(rm,len(i)):
	sum1 += avgx[j]
	sum2 += avgx2[j]
	sum3 += action[j]
	#sum4 += KE[j]
	sum5 += delta_h[j]
	#sum5 += avgx4[j]

	
	sum12 += avgx[j] * avgx[j]
	sum22 += avgx2[j] * avgx2[j]
	sum32 += action[j] * action[j]
	#sum42 += KE[j] * KE[j]
	sum52 += delta_h[j] * delta_h[j]

	k=j+1-rm
	Mavgx.append(sum1/k)
	Mavgx2.append(sum2/k)
	Maction.append(sum3/k)
	#MKE.append(sum4/k)
	Mdelta_h.append((sum5/k))
	Mavgx4.append(sum5/k)


del avgx2err[:rm]


iter =len(i)
length=2000
mu=1
f=1
oscillator_flip=1

mean = Mavgx[-1]
sum=0
sum1=0
temp1=0
error=[]
data=[]
std1=[]
std2=[]
l=7


po=[]
#Plotting
g=plt.figure(figsize=(8,6))
gs1 = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
gx2=plt.subplot(gs1[0])
plt.title('Average X')
plt.ylabel('<X>')
#gx2.axhline(y=0.0,lw=1,color='red',label='Theoretical AvgX')
gx2.plot(Mavgx,label='Simulated AvgX')
#gx2.set_ylim([-0.01,0.01])
plt.legend(loc='lower right')

gx=plt.subplot(gs1[1])
plt.ylabel('<X> error')
plt.xlabel('Iterations')
gx.plot(avgxerr,label='Simulated AvgX Error',color='green')
if(oscillator_flip==1):
	g.savefig("graphs/Average_X_Harmonic_"+str(iter) +"_"+str(length)+"_"+str(mu)+".pdf")
else:
	g.savefig("graphs/Average_X_Anharmonic_"+str(iter) +"_"+str(length)+"_"+str(f)+".pdf")




h=plt.figure(figsize=(8,6))
gs2 = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
hx2=plt.subplot(gs2[0])
plt.title('Average X^2')
plt.ylabel('<X^2>')
#hx2.axhline(y=0.4472135955,lw=1,color='red',label='Theoretical AvgX^2')
hx2.plot(Mavgx2,label='Simulated AvgX^2')
#hx2.set_ylim([0.4,0.48])
plt.legend(loc='center right')

hx=plt.subplot(gs2[1])
plt.ylabel('<X^2> error')
plt.xlabel('Iterations')
hx.plot(avgx2err,label='Simulated AvgX^2 Error',color='green')
if(oscillator_flip==1):
	h.savefig("graphs/Average_X^2_Harmonic_"+str(iter) +"_"+str(length)+"_"+str(mu)+".pdf")
else:
	h.savefig("graphs/Average_X^2_Anharmonic_"+str(iter) +"_"+str(length)+"_"+str(f)+".pdf")


o=plt.figure()
plt.ylabel('<X^4>')
plt.xlabel('Monte Carlo Iterations')
plt.title('Average X^4')
#hx2=h.add_subplot(211)
plt.plot(po,label='Simulated AvgX^4')
#plt.axhline(y=0.4472135955,lw=1,color='red',label='Theoretical AvgX^2')
plt.legend(loc='center right')
#hx=h.add_subplot(212)
#hx.plot(avgx2err,label='Simulated AvgX^2 Error')
if(oscillator_flip==1):
	o.savefig("graphs/Average_X^4_Harmonic_"+str(iter) +"_"+str(length)+"_"+str(mu)+".pdf")
else:
	o.savefig("graphs/Average_X^4_Anharmonic_"+str(iter) +"_"+str(length)+"_"+str(mu)+".pdf")

k=plt.figure()
axes=plt.gca()
plt.xlabel('Monte Carlo Iterations')
plt.ylabel('Action Per Lattice Site(S)')
plt.title('Average Action')
plt.plot(Maction,label='Simulated Action')
#plt.axhline(y=0.5,lw=1,color='red',label='Theoretical Action')
plt.legend(loc='lower right')
#axes.set_ylim([0.49,0.51])
if(oscillator_flip==1):
	k.savefig("graphs/Average_Action_Harmonic_"+str(iter) +"_"+str(length)+"_"+str(mu)+".pdf")
else:
	k.savefig("graphs/Average_Action_Anharmonic_"+str(iter) +"_"+str(length)+"_"+str(f)+".pdf")


l=plt.figure()
axes=plt.gca()
plt.xlabel('Monte Carlo Iterations')
plt.ylabel('Kinetic Energy Per Lattice Site')
plt.title('Average Kinetic Energy')
plt.plot(MKE,label='Simulated Average Kinetic Energy')
#plt.axhline(y=0.5,lw=1,color='red',label='Theoretical Average Kinetic Energy')
plt.legend(loc='center right')
#axes.set_ylim([0.49,0.51])
if(oscillator_flip==1):
	l.savefig("graphs/Average_KE_Harmonic_"+str(iter) +"_"+str(length)+"_"+str(mu)+".pdf")
else:
	l.savefig("graphs/Average_KE_Anharmonic_"+str(iter) +"_"+str(length)+"_"+str(f)+".pdf")

m=plt.figure()
plt.xlabel('Monte Carlo Iterations')
plt.ylabel('<(Delta_H)>')
plt.title('<(Delta_H)> vs Iterations')
plt.plot(Mdelta_h,label='Simulated <(Delta_H)>')
#plt.axhline(y=1,lw=1,color='red',label='Theoretical <exp(Delta_H)>')
plt.legend(loc='center right')
if(oscillator_flip==1):
	m.savefig("graphs/Average_DeltaH_Harmonic_"+str(iter) +"_"+str(length)+"_"+str(mu)+".pdf")
else:
	m.savefig("graphs/Average_DeltaH_Anharmonic_"+str(iter) +"_"+str(length)+"_"+str(f)+".pdf")

n=plt.figure()
#w = 1.00124922
w=1
x = np.linspace(-3,3,200) # 100 linearly spaced numbers
y = ((w/(3.141592654))**0.5)*np.exp(-w*(x**2))
plt.plot(x,y)
plt.hist(dataX,bins=50,normed=1)
plt.xlabel("x")
plt.ylabel("|psi|^2")
plt.title("Anharmonic Wavefunction for mu=-4 lamba =0.1")
if(oscillator_flip==1):
	n.savefig("graphs/Average_Wavefunction_Harmonic_"+str(iter) +"_"+str(length)+"_"+str(mu)+".pdf")
else:
	n.savefig("graphs/Average_Wavefunction_Anharmonic_"+str(iter) +"_"+str(length)+"_"+str(f)+".pdf")
