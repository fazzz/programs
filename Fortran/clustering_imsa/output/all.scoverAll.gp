#!/usr/local/bin/gnuplot

set terminal postscript eps enhanced background rgb 'white' color size 24cm , 8cm  "Times" 16
set output "all.scoverAll.eps"

set encoding iso_8859_1
set size square

set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0

set tics out

set multiplot layout 1,3

set ylabel "score(sim)" offset +1,0

set xlabel "score(col)"

set origin 0.1, 0.1
set size   0.25, 0.8

set xrange [0:100]
set yrange [0:100]

p "all.scoverAll" u 2:5 w p pt 7 lc rgb "black" notitle, \
  "all.scoverAll" u 2:8 w p pt 7 lc rgb "red" notitle, \
  "all.scoverAll" u 2:11 w p pt 7 lc rgb "green" notitle, \
  x w l lw 2 lc rgb "black" notitle

set origin 0.40, 0.1
set size   0.25, 0.8

set xlabel "score(dep)"
p "all.scoverAll" u 3:6 w p pt 7 lc rgb "black" notitle, \
  "all.scoverAll" u 3:9 w p pt 7 lc rgb "red" notitle, \
  "all.scoverAll" u 3:12 w p pt 7 lc rgb "green" notitle, \
   x w l lw 2 lc rgb "black" notitle

set origin 0.70, 0.1
set size   0.25, 0.8

set xrange [0:10]
set yrange [0:10]

set xlabel "score(stamp)"
p "all.scoverAll" u 4:7 w p pt 7 lc rgb "black" notitle, \
  "all.scoverAll" u 4:10 w p pt 7 lc rgb "red" notitle, \
  "all.scoverAll" u 4:13 w p pt 7 lc rgb "green" notitle, \
   x w l lw 2 lc rgb "black" notitle

quit

