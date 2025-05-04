# gnuplot
# load "plot_spectralIrradiance.plot"

set style line 1 lt rgb "#7F7F7F" lw 3 pt 6

# plot data
set xlabel "Lambda"
set ylabel "Spectral Irradiance"
set arrow from 680,0 to 680,2.5 nohead linecolor "red"
set arrow from 550,0 to 550,2.5 nohead linecolor "green"
set arrow from 440,0 to 440,2.5 nohead linecolor "blue"
plot "plot_spectralIrradiance.dat" using 1:2 with linespoints ls 1
