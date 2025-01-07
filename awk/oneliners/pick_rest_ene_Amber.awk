#!/bin/gawk -f

BEGIN{
    i=0
  printf("step restene \n")
}

$1 ~ /(EELEC)/{
  ++i
  if ((i % intv)==0) {
      printf("%d %lf\n",i,$9)
  }
}
