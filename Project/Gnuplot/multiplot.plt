set term png
set nokey


set title "Harmonic Oscillator"
set xlabel "Leap Frog Iterations"
set ylabel "Hamiltonian"
set output "Multi_out_Harmonic_2000.png"
plot 'HMC_H_No_Fourier_2000' with lines


set title "Anharmonic Oscillator"
set output "Multi_out_Anharmonic_2000.png"
plot 'HMC_H_No_Fourier_2000_Anharmonic' with lines

unset multiplot
set term x11
set key
