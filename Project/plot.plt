set term png 

set output "/home/s1403955/MPhys/graphs/Leap Frog vs Hamiltonian Fourier-Harmonic"
set title "Leap Frog Iterations vs Hamiltonian"
set xlabel "Leap Frog Iterations"
set ylabel "Hamiltonian"

plot "HMC_H" with lines

set term x11
