set terminal postscript enhanced background rgb 'white' color size 12cm , 12cm  "Times" 12
set output "contribution_apo.eps"

set encoding iso_8859_1

set tics out

set size square

set style fill solid border lc rgb "black"

set tics font "Times-Roman,12"
set xlabel font "Times-Roman,20"
set ylabel font "Times-Roman,20"
set label font "Times-Roman,12"
set key font "Times-Roman,12"

cont_pc="cont_apo.txt"

set xlabel "Index of the Eigenvalue"
set xrange [-1:16]
set ylabel "Contribution(%)"
set yrange [-1:101]

plot cont_pc     u 1:3 w lp pt 7 lc rgb "red"   notitle, \

quit

