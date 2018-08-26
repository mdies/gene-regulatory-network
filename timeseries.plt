set terminal postscript enhanced lw 1.5 landscape color
set output 'timeseries.ps'
set xlabel "Time (h)" font "Helvetica, 18"
set tics in
set xtics mirror
set ytics mirror
set title " Stochastic simulation for p53-Mdm2 system " font "Helvetica, 18"
set ylabel " Number of molecules " font "Helvetica, 18"
set xrange [ 0.0 : 43200.0 ] noreverse nowriteback
set xtics norotate ("0" 0, "3" 10800, "6" 21600, "9" 32400, "12" 43200)
plot "out_p53.dat" using 1:2 with line lt 1 lc rgb "red" lw 2 title "p53",  "out_p53.dat" using 1:3 with line lt 3 lc rgb "blue" lw 2 title "Mdm2"
