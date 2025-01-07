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
sp "../FASYS2/trunk/FASYS_w_-180-180_-180_180.map" w l
set term table 
set output "../FASYS2/trunk/FASYS_w_-180-180_-180_180.map.table"
replot
#set output
set term x11
set nokey
p "../FASYS2/trunk/FASYS_w_-180-180_-180_180.map.table" i 0:0 u 1:2 w l
rep  "../FASYS2/trunk/FASYS_w_-180-180_-180_180.map.table" i 1:1 u 1:2 w l
rep  "../FASYS2/trunk/FASYS_w_-180-180_-180_180.map.table" i 2:2 u 1:2 w l
rep  "../FASYS2/trunk/FASYS_w_-180-180_-180_180.map.table" i 3:3 u 1:2 w l 
rep "../FASYS2/trunk/FASYS_w_-180-180_-180_180.map.table" i 4:4 u 1:2 w l
rep  "../FASYS2/trunk/FASYS_w_-180-180_-180_180.map.table" i 5:5 u 1:2 w l
rep  "../FASYS2/trunk/FASYS_w_-180-180_-180_180.map.table" i 6:6 u 1:2 w l
rep  "../FASYS2/trunk/FASYS_w_-180-180_-180_180.map.table" i 7:7 u 1:2 w l
rep  "../FASYS2/trunk/FASYS_w_-180-180_-180_180.map.table" i 8:8 u 1:2 w l
rep  "../FASYS2/trunk/FASYS_w_-180-180_-180_180.map.table" i 9:9 u 1:2 w l
rep  "../FASYS2/trunk/FASYS_w_-180-180_-180_180.map.table" i 10:10  u 1:2 w l
rep  "FASYS_w_-60_60_120_60.dtrj_fin" u 1:2 w l lw 2
rep  "FASYS_w_-60_60_120_60.dtrj_ini" u 1:2 w l lw 2
pause -1
