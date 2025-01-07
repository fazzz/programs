#!/bin/gawk -f

BEGIN{
    i=0
  printf("step T K N \n")
}


$1 ~ /(kinetic_energy8)/{
  ++i
  if ((i % intv)==0) {
      K=$3
  }
}
$1 ~ /(T_kelvin)/{
  if ((i % intv)==0) {
      T=$3
      N=2.0/1.98723e-3*K/T
      printf("%d %lf %lf %lf \n",i,T,K,N)
  }
}
