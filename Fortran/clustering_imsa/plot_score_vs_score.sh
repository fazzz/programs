#!/bin/sh

opt=( dummy scorefil )
nopt=${#opt[*]}
if [ $# -le `expr ${nopt} - 2` ]; then
    echo "USAGE: $0" ${opt[*]:1:${nopt}}
    echo $*
    exit
fi

num=1
while [ $num -le `expr ${nopt} - 1` ]; do
    eval ${opt[$num]}=$1
    shift 1
    num=`expr $num + 1`
done

cat <<EOF > ${scorefil}.gp
#!/usr/local/bin/gnuplot

set terminal postscript eps enhanced background rgb 'white' color size 21cm , 7cm  "Times" 14
set output "${scorefil}.eps"

set encoding iso_8859_1

set tics out

set size square

set multiplot layout 1,3

set ylabel "score(sim)" offset +1,0

set xlabel "score(col)"

p "${scorefil}" u 2:9 w p pt 7 lc rgb "black" notitle, \\
  "${scorefil}" u 2:10 w p pt 6 lc rgb "red" notitle, \\
  "${scorefil}" u 2:11 w p pt 6 lc rgb "green" notitle

set xlabel "score(dep)"

p "${scorefil}" u 3:9 w p pt 7 lc rgb "black" notitle, \\
  "${scorefil}" u 3:10 w p pt 6 lc rgb "red" notitle, \\
  "${scorefil}" u 3:11 w p pt 6 lc rgb "green" notitle

set xlabel "score(stamp)"

p "${scorefil}" u 4:9 w p pt 7 lc rgb "black" notitle, \\
  "${scorefil}" u 4:10 w p pt 6 lc rgb "red" notitle, \\
  "${scorefil}" u 4:11 w p pt 6 lc rgb "green" notitle

quit
EOF

chmod +x ${scorefil}.gp
./${scorefil}.gp




