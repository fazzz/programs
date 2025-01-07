#!/bin/sh
                                                                                                                                                                                  if [ -z "$3" ]
then
echo Usage: ./$( basename $0 ) ph_map p h
exit
fi
                                                                                                                                                                                              
ph_map=$1
p=$2
h=$3

cat << eof1 > plot_ph.gp
set xrange [-180.0:180.0]
set yrange [-180.0:180.0]
pl "${ph_map}" u $p $h
pause -1
eof1

gnuplot plot_ph.gp
