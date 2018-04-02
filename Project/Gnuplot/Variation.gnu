set term png
set output 'Variation.png'

set title 'HMC vs Fourier Accelerated HMC'
set xlabel 'Monte Carlo Iterations'
set ylabel 'Difference from Theoritical Value'

set xrange[0:20]
set yrange[-0.0001:0.01]

plot '../variation.dat' using 1:2 title 'Fourier Accelerated' with lines, '../variation.dat' using 1:3 title 'Non Fourier Accelerated' with lines
