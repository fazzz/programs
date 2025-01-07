set terminal postscript enhanced background rgb 'white' color size 12cm, 12cm "Times" 10
set output "pmf.eps"
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
set xrange [-1:3]
set yrange [0:3]
sp "pmf_pca.txt" w l notitle
quit
