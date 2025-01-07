gawk '{if ($3!="****"){print $1, " ", $2, " ", $3/0.6, " " }else{print $0}}' < >
