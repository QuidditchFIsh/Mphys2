set term png 

set title "Leaf Frog Hamiltonian Evolution-Harmonic"
set xlabel "Leap Frog Iterations"
set ylabel "Hamiltonian"

set output "/home/s1403955/MPhys/graphs/Leap Frog vs Hamiltonian-Harmonic"
plot "HMC_H_Harmonic" with lines

set title "Leaf Frog Hamiltonian Evolution-Anharmonic"
set output "/home/s1403955/MPhys/graphs/Leap Frog vs Hamiltonian-Anharmonic"
plot "HMC_H_Anharmonic" with lines

set term x11


