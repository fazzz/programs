#!/bin/sh

if [ -z "$1" ]
then
    echo Usage: ./$( basename $0 ) filename
    exit
fi
           
filename=$1

cat << eof1 > plot_contour.gp
set size square
set contour base
set cntrparam cubicspline
set cntrparam levels 20
#set palette rgbformulae 33,13,10
set nosurface
set view 0,0
#set nokey
#set tics font "Helvetica,15"
set xrange [-180.0:180.0]
set yrange [-180.0:180.0]
sp "${filename}" w l
set term table 
set output "${filename}.table"
replot
set term x11
p "${filename}.table" i 0:0 u 1:2 w l, \
  "${filename}.table" i 1:1 u 1:2 w l, \
  "${filename}.table" i 2:2 u 1:2 w l, \
  "${filename}.table" i 3:3 u 1:2 w l, \
  "${filename}.table" i 4:4 u 1:2 w l, \
  "${filename}.table" i 5:5 u 1:2 w l, \
  "${filename}.table" i 6:6 u 1:2 w l, \
  "${filename}.table" i 7:7 u 1:2 w l, \
  "${filename}.table" i 8:8 u 1:2 w l, \
  "${filename}.table" i 9:9 u 1:2 w l, \
  "${filename}.table" i 10:10  u 1:2 w l
pause -1
eof1

gnuplot plot_contour.gp
