#!/bin/gawk


$1 ~ /ATOM/ { if ( $3 ~ /1(H*)/  ) { print $0 } }

$1 ~ /ENDMDL/ { exit }
