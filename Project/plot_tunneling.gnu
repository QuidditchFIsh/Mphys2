set term png 
set output 'Tunneling.png'
set key off

set title 'Position of Oscillators along Lattice'
set xlabel 'time [t]'
set ylabel 'position [x]'

plot 'HMC_X1.dat' with lines
