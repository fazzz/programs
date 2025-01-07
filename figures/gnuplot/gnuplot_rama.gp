#!/bin/sh

if [ -z "$1" ]
then
    echo Usage: ./$( basename $0 ) rama_id 
    exit
fi

                                                                                           
rama_id=$1

cat << eof1 > plot_rama.gp
set pm3d
set size square
set contour
set cntrparam cubicspline
set cntrparam levels 20
#set palette defined (0 "black", 1.0 "red", 2.0 "orange", 4.0 "yellow", 4.5 "white" )
#set palette rgbformulae 22,13,-31 #-nijiiro-
set palette rgbformulae 33,13,10
set nosurface
set view 0,0
#set nokey
#set tics font "Helvetica,15"
set xrange [-180.0:180.0]
set yrange [-180.0:180.0]
sp "${rama_id}" w l
pause -1
eof1

gnuplot plot_rama.gp
