#
# Gnuplot commands
#
# Usage: gnuplot ./plot_bands.gp
#
#set xr [1:71]
#set yr [-4:3]
set grid
plot 'bands_plot.dat' u 1:2 w l lc rgb 'red'
pause -1
