# gnuplot
# load "plot_density.plot"

set style line 1 lt rgb "#0000FF" lw 3 pt 6
set style line 2 lt rgb "#7F7F7F" lw 3 pt 6
set style line 3 lt rgb "#FF0000" lw 3 pt 6

set xlabel "h (0-100000 meters)"
set ylabel "density"

plot "plot_density.dat" using 1:2 with linespoints ls 1 title 'pR',  \
     "plot_density.dat" using 1:3 with linespoints ls 2 title 'pM',  \
     "plot_density.dat" using 1:4 with linespoints ls 3 title 'pO'
