# gnuplot
# load "plot_avgT.plot"

# plot data
set xlabel "phi (0-180 degrees)"
set ylabel "h (0-100000 meters)"
set zlabel "avgT"
set pm3d
set hidden3d
splot "plot_avgT.dat" matrix
