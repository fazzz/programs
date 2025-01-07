
set terminal postscript enhanced color size 10cm , 9cm "Arial" 12
set output "rdf.eps"

set grid
set yrange [0:5]
set xlabel "radius(nm)"

set ylabel "g(r)"
set yrange [0:3]
plot "rdf.xvg" u 1:2 w l lc rgb "black" notitle

quit
