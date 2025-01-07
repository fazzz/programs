#!/bin/gawk -f

BEGIN{
    i=0
  printf("step T K N \n")
}

$1 ~ /(Etot)/{
  ++i
  if ((i % intv)==0) {
      K=$6
  }
}
$1 ~ /(NSTEP)/{
  if ((i % intv)==0) {
      T=$9
      N=2.0/1.98723e-3*K/T
      printf("%d %lf %lf %lf \n",i,T,K,N)
  }
}
