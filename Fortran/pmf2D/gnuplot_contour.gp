#!/bin/sh

if [ -z "$1" ]
then
    echo Usage: ./$( basename $0 ) filename minx maxx miny maxy output
    exit
fi
           
filename=$1

cat << eof1 > plot_contour.gp
set terminal postscript enhanced background rgb 'white' color size 12cm, 12cm "Times" 10
set output "$6"
set pm3d
set size square
set contour
set cntrparam cubicspline
set cntrparam levels 20
#set palette rgbformulae 33,13,10
set palette rgbformulae 33,13,10
set nosurface
set view 0,0
set nokey
set xrange [$2:$3]
set yrange [$4:$5]
sp "${filename}" w l notitle
quit
eof1

gnuplot plot_contour.gp

