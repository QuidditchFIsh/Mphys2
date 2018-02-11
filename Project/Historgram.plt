clear
reset
set key off
set border 3

# Add a vertical dotted line at x=0 to show centre (mean) of distribution.
set yzeroaxis

# Each bar is half the (visual) width of its x-range.
set boxwidth 0.05 absolute
set style fill solid 1.0 noborder

bin_width = 0.1;
sum = 0

bin_number(x) = floor(x/bin_width)

rounded(x) = bin_width * ( bin_number(x) + 0.5 )

s(x) = ((sum = sum+1),0)

set term png
set output 'WaveFunction1.png'
plot 'HMC_X.dat' using ($1):(s($1))

set output 'WaveFunction.png'
plot 'HMC_X.dat' using (rounded($1)):(1.0/(sum*bin_width)) smooth frequency with boxes


