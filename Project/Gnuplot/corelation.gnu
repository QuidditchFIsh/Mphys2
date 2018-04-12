set term png
set output 'Corelation.png'

set xrange[0:6]

set title 'Corelation function for Harmonic Oscillator'
set xlabel 't'
set ylabel '<x(0)x(t)>'


plot 'corelation_out.txt' u 1:2 with lines title 'Simulated Corelation Function', 'corelation_theory_out.txt'u 1:2 title 'Theoritical Corelation Function' 
