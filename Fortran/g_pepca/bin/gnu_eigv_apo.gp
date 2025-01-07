set terminal postscript enhanced background rgb 'white' color size 16cm , 16cm  "Times" 12
set output "eig_apo.eps"

set encoding iso_8859_1

set tics out

set size square

set multiplot layout 1,2

set style fill solid border lc rgb "black"

set tics font "Times-Roman,12"
set xlabel font "Times-Roman,12"
set ylabel font "Times-Roman,12"
set label font "Times-Roman,12"
set key font "Times-Roman,12"

eigv="eig_apo.txt"

set xlabel "Index of the Factors"
set xrange [0:16]
set ylabel "Components of 1^{st} PC"
set yrange [-1:1]

plot eigv u 1:2 w impulses lw 4 lc rgb "red" notitle

set ylabel "Components of 2^{nd} PC"

plot eigv u 1:3 w impulses lw 4 lc rgb "red" notitle

quit


