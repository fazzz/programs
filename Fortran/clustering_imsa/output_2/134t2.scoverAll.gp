#!/usr/local/bin/gnuplot

set terminal postscript eps enhanced background rgb 'white' color size 20cm , 16cm  "Times" 20
set output "output_2/134t2.scoverAll.eps"

set encoding iso_8859_1

set tics out

set size square

set multiplot layout 1,3

set ylabel "score(sim)" offset +1,0

set xlabel "score(col)"
p "output_2/134t2.scoverAll" u 2:9 w p pt 7 lc rgb "black" notitle,   "output_2/134t2.scoverAll" u 2:10 w p pt 7 lc rgb "red" notitle,   "output_2/134t2.scoverAll" u 2:11 w p pt 7 lc rgb "green" notitle
set xlabel "score(dep)"
p "output_2/134t2.scoverAll" u 3:9 w p pt 7 lc rgb "black" notitle,   "output_2/134t2.scoverAll" u 3:10 w p pt 7 lc rgb "red" notitle,   "output_2/134t2.scoverAll" u 3:11 w p pt 7 lc rgb "green" notitle
set xlabel "score(stamp)"
p "output_2/134t2.scoverAll" u 4:9 w p pt 7 lc rgb "black" notitle,   "output_2/134t2.scoverAll" u 4:10 w p pt 7 lc rgb "red" notitle,   "output_2/134t2.scoverAll" u 4:11 w p pt 7 lc rgb "green" notitle

quit
