#!/bin/gawk -f

BEGIN{
    i=0
  printf("step ene kine pote \n")
}

$1 ~ /(Etot)/{
  ++i
  if ((i % intv)==0) {
      printf("%d %lf  %lf  %lf\n",i,$3,$6,$9)
  }
}



