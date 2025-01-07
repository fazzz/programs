#!/bin/gawk

$1 ~ /ATOM/ { if ($5 == A)  {print $0} }


