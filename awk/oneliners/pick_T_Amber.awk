#!/bin/gawk -f

BEGIN{
    i=0
    printf("step Temp \n") 
}

$1 ~ /(NSTEP)/{
  ++i
  if ((i % intv)==0) {
      printf("%d %lf \n",i,$9)
  }
}

