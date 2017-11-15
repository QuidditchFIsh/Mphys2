from numpy import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from array import *
import linecache
import statistics 


data1=[];data2=[];data3=[];data4=[];data5=[]
data=[]
var0=0
avgx=0
avgx1=0
avgx2=0
var1=0
for j in range(0,9000):
	x  =linecache.getline("HMC_X.dat",1+j)
	x   = x.split(" ")
	del x   [-1]
	x    = [float(a) for a in x   ]

	y  =linecache.getline("HMC_X.dat",2+j)
	y   = y.split(" ")
	del y   [-1]
	y    = [float(b) for b in y   ]

	z  =linecache.getline("HMC_X.dat",3+j)
	z   = z.split(" ")
	del z   [-1]
	z    = [float(b) for b in z   ]

	aa  =linecache.getline("HMC_X.dat",4+j)
	aa   = aa.split(" ")
	del aa   [-1]
	aa    = [float(aaa) for aaa in aa   ]

	b1  =linecache.getline("HMC_X.dat",5+j)
	b1  = b1.split(" ")
	del b1   [-1]
	b1    = [float(bb) for bb in b1   ]

	for i in range(0,len(x)):
		avgx2 += x[i]*x[i]
		avgx  += x[i]
	avgx2 = avgx2/len(x)
	avgx = avgx/len(x)
	var0= avgx2 - avgx**2
	data1.append(1.0)
	avgx=0
	agvx2=0

	for i in range (0,len(x)):
		avgx2 += x[i]*y[i]
		avgx  += x[i]
		avgx1 +=y[i]
	avgx2 = avgx2/len(x)
	avgx = avgx/len(x)
	avgx1 = avgx1/len(x)
	var1= avgx2 - (avgx*avgx1)
	data2.append(var1/var0)
	avgx=0
	agvx2=0
	avgx1=0

	for i in range (0,len(x)):
		avgx2 += x[i]*z[i]
		avgx  += x[i]
		avgx1 +=z[i]
	avgx2 = avgx2/len(x)
	avgx = avgx/len(x)
	avgx1 = avgx1/len(x)
	var1= avgx2 - (avgx*avgx1)
	data3.append(var1/var0)
	avgx=0
	agvx2=0

	for i in range (0,len(x)):
		avgx2 += x[i]*aa[i]
		avgx  += x[i]
		avgx1 +=aa[i]
	avgx2 = avgx2/len(x)
	avgx = avgx/len(x)
	avgx1 = avgx1/len(x)
	var1= avgx2 - (avgx*avgx1)
	data4.append(var1/var0)
	avgx=0
	agvx2=0

	for i in range (0,len(x)):
		avgx2 += x[i]*b1[i]
		avgx  += x[i]
		avgx1 +=b1[i]
	avgx2 = avgx2/len(x)
	avgx = avgx/len(x)
	avgx1 = avgx1/len(x)
	var1= avgx2 - (avgx*avgx1)
	data5.append(var1/var0)

	avgx=0
	agvx2=0
	var0=0
	avgx=0
	avgx1=0
	avgx2=0
	var1=0


data.append(np.mean(data1))
data.append(np.mean(data2))
data.append(np.mean(data3))
data.append(np.mean(data4))
data.append(np.mean(data5))

print(data)

plt.plot(data)
plt.show()






