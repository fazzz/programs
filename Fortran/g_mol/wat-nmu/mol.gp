
set terminal postscript enhanced color size 10cm , 9cm "Arial" 12
set output "mol.eps"

set grid
set xrange [0:1]
set xlabel "time(ns)"

set ylabel "M"
set ytics 0.1
set yrange [2.5:3.5]
plot "mol.txt" u ($1/100):4 w l lc rgb "black" notitle

quit
