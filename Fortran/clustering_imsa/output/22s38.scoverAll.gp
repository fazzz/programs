#!/usr/local/bin/gnuplot

set terminal postscript eps enhanced background rgb 'white' color size 21cm , 7cm  "Times" 14
set output "output/22s38.scoverAll.eps"

set encoding iso_8859_1

set tics out

set size square

set multiplot layout 1,3

set ylabel "score(sim)" offset +1,0

set xlabel "score(col)"

p "output/22s38.scoverAll" u 2:9 w p pt 7 lc rgb "black" notitle, \
  "output/22s38.scoverAll" u 2:10 w p pt 6 lc rgb "red" notitle, \
  "output/22s38.scoverAll" u 2:11 w p pt 6 lc rgb "green" notitle

set xlabel "score(dep)"

p "output/22s38.scoverAll" u 3:9 w p pt 7 lc rgb "black" notitle, \
  "output/22s38.scoverAll" u 3:10 w p pt 6 lc rgb "red" notitle, \
  "output/22s38.scoverAll" u 3:11 w p pt 6 lc rgb "green" notitle

set xlabel "score(stamp)"

p "output/22s38.scoverAll" u 4:9 w p pt 7 lc rgb "black" notitle, \
  "output/22s38.scoverAll" u 4:10 w p pt 6 lc rgb "red" notitle, \
  "output/22s38.scoverAll" u 4:11 w p pt 6 lc rgb "green" notitle

quit
