
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]) {
  int i,j,num;
  double pi,interval,x,y;
  
  pi=acos(-1.0);
  num=20;
  interval=(2*pi)/num;
  for (i=0;i<=num;++i) {
    for (j=0;j<num;++j) {
      printf("%e %e %e\n",x=i*interval-pi,y=j*interval-pi,sin(x)-cos(y));
    }
    printf("\n");
  }
  
  return 0;
}

