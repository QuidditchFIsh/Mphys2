set term png
set output 'wavefunction.png'
set key off

binwidth = 0.05
set boxwidth binwidth
sum = 0

set title 'Wavefunction for Harmonic Oscillator a=0.05 mu=1.0'
set xlabel 'position [x]'
set ylabel '|phi(x)|^2'

s(x)          = ((sum=sum+1), 0)
bin(x, width) = width*floor(x/width) + binwidth/2.0

plot "HMC_X.dat" u ($1):(s($1))
set output 'WF1.png'
w=1.0
f(x) = ((w/(3.141592654))**0.5)*exp(-w*(x**2))
plot "HMC_X.dat" u (bin($1, binwidth)):(1.0/(binwidth*sum)) smooth freq w boxes,f(x) w l
