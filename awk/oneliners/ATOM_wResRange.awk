#!/bin/gawk

$1 ~ /ATOM/ { if ($6 >= Ini && $6 <= Fin )  {print $0} }


