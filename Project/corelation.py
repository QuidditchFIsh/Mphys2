from numpy import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from array import *
import linecache
import statistics 


data=[]
length=1


one=0
zero=0
two=0
three=0
four=0
five=0
six=0
seven=0
eight=0
nine=0
ten=0

sum1 = 0
sum2 = 0
sum3 = 0
sum4 = 0
sum5 = 0
sum6 = 0
sum7 = 0
sum8 = 0
sum9 = 0
sum10 = 0
sum11 =0

ones=0
zeros=0
twos=0
threes=0
fours=0
fives=0
sixs=0
sevens=0
eights=0
nines=0
tens=0

length = 2000
start = 4000
for i in range(start,start+length):
	x  =linecache.getline("HMC_X.dat",1+i)
	x   = x.split(" ")
	del x   [-1]
	xn    = [float(a) for a in x   ]

	y  =linecache.getline("HMC_X.dat",2+i)
	y   = y.split(" ")
	del y   [-1]
	yn    = [float(b) for b in y   ]

	z  =linecache.getline("HMC_X.dat",3+i)
	z   = z.split(" ")
	del z   [-1]
	zn    = [float(b) for b in z   ]

	aa  =linecache.getline("HMC_X.dat",4+i)
	aa   = aa.split(" ")
	del aa   [-1]
	aan    = [float(aaa) for aaa in aa   ]

	b1  =linecache.getline("HMC_X.dat",5+i)
	b1  = b1.split(" ")
	del b1   [-1]
	b1n    = [float(bb) for bb in b1   ]

	c1  =linecache.getline("HMC_X.dat",6+i)
	c1   = c1.split(" ")
	del c1   [-1]
	c1n    = [float(b1) for b1 in c1   ]

	z1  =linecache.getline("HMC_X.dat",7+i)
	z1   = z1.split(" ")
	del z1   [-1]
	z1n    = [float(b) for b in z1   ]

	aa1  =linecache.getline("HMC_X.dat",8+i)
	aa1   = aa1.split(" ")
	del aa1   [-1]
	aa1n    = [float(aaa) for aaa in aa1   ]

	b11  =linecache.getline("HMC_X.dat",9+i)
	b11  = b11.split(" ")
	del b11   [-1]
	b11n    = [float(bb) for bb in b11   ]

	c11  =linecache.getline("HMC_X.dat",10+i)
	c11   = c11.split(" ")
	del c11   [-1]
	c11n    = [float(b1) for b1 in c11   ]

	c12  =linecache.getline("HMC_X.dat",11+i)
	c12   = c12.split(" ")
	del c12   [-1]
	c12n    = [float(b1) for b1 in c12   ]




	for i in range(0,len(x)-10):
		zero +=xn[i]*xn[i]
		one+=xn[i]*xn[i+1]
		two +=xn[i]*xn[i+2]
		three +=xn[i]*xn[i+3]
		four +=xn[i]*xn[i+4]
		five +=xn[i]*xn[i+5]
		six +=xn[i] *xn[i+6]
		seven +=xn[i] *xn[i+7]
		eight +=xn[i] *xn[i+8]
		nine +=xn[i] *xn[i+9]
		ten +=xn[i] *xn[i+10]

	sum1+=(zero/(len(x)-10))
	sum2+=(one/(len(x)-10))
	sum3+=(two/(len(x)-10))
	sum4+=(three/(len(x)-10))
	sum5+=(four/(len(x)-10))
	sum6+=(five/(len(x)-10))
	sum7+=(six/(len(x)-10))
	sum8+=(seven/(len(x)-10))
	sum9+=(eight/(len(x)-10))
	sum10+=(nine/(len(x)-10))
	sum11+=(ten/(len(x)-10))


	one=0
	zero=0
	two=0
	three=0
	four=0
	five=0
	six=0
	seven=0
	eight=0
	nine=0
	ten=0

data.append((sum1/sum1))
data.append((sum2/sum1))
data.append((sum3/sum1))
data.append((sum4/sum1))
data.append((sum5/sum1))
data.append((sum6/sum1))
data.append((sum7/sum1))
data.append((sum8/sum1))
data.append((sum9/sum1))
data.append((sum10/sum1))
data.append((sum11/sum1))

print(data)
plt.plot(data)
#plt.yscale('log')
plt.show()


'''



	zeros+=(zero/zero)
	ones+=(one/zero)
	twos+=(two/zero)
	threes+=(three/zero)
	fours+=(four/zero)
	fives+=(five/zero)
	six+=(six/zero)
	seven+=(seven/zero)
	eights+=(eight/zero)
	nine+=(nine/zero)
	ten+=(ten/zero)

	one=0
	zero=0
	two=0
	three=0
	four=0
	five=0
	six=0
	seven=0
	eight=0
	nine=0
	ten=0

ones=ones/length
zeros=zeros/length
twos=twos/length
threes=threes/length
fours=fours/length
fives=fives/length
sixs=sixs/length
sevens=sevens/length
eights=eights/length
nines=nines/length
tens=tens/length

print(zeros)
print(ones)
print(twos)
print(threes)
print(fours)
print(fives)
print(sixs)
print(sevens)
print(eights)
print(nines)
print(tens)
'''
