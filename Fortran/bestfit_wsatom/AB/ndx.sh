#!/bin/sh

echo "[$1]"

for i in `seq $2 $3`; do
    printf "%4d " $i 

    j=`expr $i % 15`
    if [ $j == 0 ]; then
	printf "\n"
    fi
done






printf "\n"
