set term png
set output 'Anharmonic_Potential.png'
set key off

set title 'Anharmonic Potential'
set xlabel 'x'
set ylabel 'V(x) = 3(x^2 -4)^2'

set xrange[-5:5]
set yrange[0:100]

f(x) = 3*(x**2 - 4)**2

plot f(x)
