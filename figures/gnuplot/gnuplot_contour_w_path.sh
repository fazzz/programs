#!/bin/sh

if [ -z $2 ]; then
    echo "Usage: ./$( basename $0 ) filename(map) filename(path) filename(path) ......"
    exit
fi
           
filename=$1
nfile=`expr $# - 1 `
for num in `seq 0 ${nfile}`; do
    filename2[$num]=$2
    shift 1
done

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
#set output
set term x11
set nokey
p "${filename}.table" i 0:0 u 1:2 w l
rep  "${filename}.table" i 1:1 u 1:2 w l
rep  "${filename}.table" i 2:2 u 1:2 w l
rep  "${filename}.table" i 3:3 u 1:2 w l 
rep "${filename}.table" i 4:4 u 1:2 w l
rep  "${filename}.table" i 5:5 u 1:2 w l
rep  "${filename}.table" i 6:6 u 1:2 w l
rep  "${filename}.table" i 7:7 u 1:2 w l
rep  "${filename}.table" i 8:8 u 1:2 w l
rep  "${filename}.table" i 9:9 u 1:2 w l
rep  "${filename}.table" i 10:10  u 1:2 w l
eof1

num=`expr ${num} - 1`
for num2 in `seq 0  ${num}`; do
   cat << eof2 >> plot_contour.gp
rep  "${filename2[$num2]}" u 1:2 w l lw 2
eof2
done

cat << eof3 >> plot_contour.gp
pause -1
eof3

gnuplot plot_contour.gp
