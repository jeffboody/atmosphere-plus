# gnuplot
# load "plot_spectralToRGB.plot"

set style line 1 lt rgb "#FF0000" lw 3 pt 6
set style line 2 lt rgb "#00FF00" lw 3 pt 6
set style line 3 lt rgb "#0000FF" lw 3 pt 6

# plot data
set xlabel "Lambda"
set ylabel "Intensity"
set arrow from 680,0 to 680,2.5 nohead linecolor "red"
set arrow from 550,0 to 550,2.5 nohead linecolor "green"
set arrow from 440,0 to 440,2.5 nohead linecolor "blue"
set multiplot layout 3,1
plot "plot_spectralToRGB.dat" using 1:2  with linespoints ls 1 title 'R', \
     "plot_spectralToRGB.dat" using 1:3  with linespoints ls 2 title 'G', \
     "plot_spectralToRGB.dat" using 1:4  with linespoints ls 3 title 'B'
plot "plot_spectralToRGB.dat" using 1:5  with linespoints ls 1 title 'R (Sum)', \
     "plot_spectralToRGB.dat" using 1:6  with linespoints ls 2 title 'G (Sum)', \
     "plot_spectralToRGB.dat" using 1:7  with linespoints ls 3 title 'B (Sum)'
plot "plot_spectralToRGB.dat" using 1:8  with linespoints ls 1 title 'R (Peak)', \
     "plot_spectralToRGB.dat" using 1:9  with linespoints ls 2 title 'G (Peak)', \
     "plot_spectralToRGB.dat" using 1:10 with linespoints ls 3 title 'B (Peak)'
unset multiplot
