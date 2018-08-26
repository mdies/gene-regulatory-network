#!/bin/ksh

./anderson\_Model2\_Hurdle < p53\_data.dat > exec.log 2>&1
gnuplot timeseries.plt
ps2pdf timeseries.ps
