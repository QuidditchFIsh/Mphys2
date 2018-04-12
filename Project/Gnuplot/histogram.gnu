set term png 
set key off

set xrange[-3:3]

binwidth = 0.05
set boxwidth binwidth
sum = 0

s(x)          = ((sum=sum+1), 0)
bin(x, width) = width*floor(x/width) + binwidth/2.0

set output './graphs/WaveFunction_1.png'
plot "HMC_X.dat" u ($1):(s($1))

set title 'Wavefunction for Harmonic Oscillator'
set xlabel 'x'
set ylabel '|\psi(x)|^2'
w = 1.00
pi = 3.141592654
f(x) = (w/pi)**0.5 * exp(-w*(x**2))
set output './graphs/WaveFunction.png'
plot "HMC_X.dat" u (bin($1, binwidth)):(1.0/(binwidth*sum)) smooth freq w boxes,f(x) with lines
