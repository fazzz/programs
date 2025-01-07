#!/bin/gawk -f

BEGIN{
    i=0
  printf("step pote(dihed) pote(e1-4) pote(LJ1-4) pote(e1-5) pote(LJ1-5) pote(d+1-4+1-5) \n")
}

$1 ~ /(BOND)/{
  ++i
  if ((i % intv)==0) {
      p1=$9
      printf("%d %lf  ",i,$9)
  }
}
$1 ~ /(1-4)/{
  if ((i % intv)==0) {
      p2=$4
      p3=$8
      p4=$11
      printf("%lf  %lf  %lf ",$4,$8,$11)
  }
}
$1 ~ /(EELEC)/{
  if ((i % intv)==0) {
      p5=$3
      p=p1+p2+p3+p4+p5
      printf("%lf  %lf  \n",$3,p)
  }
}

# NSTEP =   100000   TIME(PS) =     100.000  TEMP(K) =     2.16  PRESS =     0.0
# Etot   =         0.0712  EKtot   =         0.0708  EPtot      =         0.0795
# BOND   =         0.0000  ANGLE   =         0.1911  DIHED      =         0.0345
# 1-4 NB =         0.2114  1-4 EEL =         0.4839  VDWAALS    =         0.0632
# EELEC  =         0.5327  EHBOND  =         0.0000  RESTRAINT  =         0.0000
