#!/bin/gawk


$1 ~ /CONECT/ { print $0 }
$1 ~ /ATOM/ { print $0 }

$1 ~ /ENDMDL/ { exit }
