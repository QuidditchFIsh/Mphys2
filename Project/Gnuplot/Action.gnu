set term png
set output 'Action_fourier.png'
set key off

set title 'Action vs Leapfrog Steps'
set xlabel 'Leapfrog Iterations'
set ylabel 'Action'
set xrange[0:5000]

plot '../HMC_H' with lines
