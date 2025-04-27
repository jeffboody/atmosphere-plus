# gnuplot
# load "plot_densityRM.plot"

set style line 1 lt rgb "#0000FF" lw 3 pt 6
set style line 2 lt rgb "#7F7F7F" lw 3 pt 6

set xlabel "h (0-100000 meters)"
set ylabel "density"

plot "plot_densityRM.dat" using 0:1 with linespoints ls 1 title 'densityR', \
     "plot_densityRM.dat" using 0:2 with linespoints ls 2 title 'densithM'
