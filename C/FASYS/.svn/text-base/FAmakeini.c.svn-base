
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"

int main(int argc, char *argv[]) {
  int i,j,k;

  double pi;
  double p1,h1,p2,h2;
  double r1[5][3];
  double l_eq=1.53;
  double a_eq;

  pi=acos(-1.0);
  a_eq=111.0/180.0*pi;

  char *inifilename;
  FILE *inifile;

  if (argc < 2) {
    printf("USAGE: %s inifilename p1 h1 p2 h2 \n",argv[0]);
    exit(1);
  }
  inifilename  = *++argv;
  inifile = efopen(inifilename,"w");
  for (i=0;i<5;++i)
    for (j=0;j<3;++j)
      r1[i][j]=0.0;
  r1[1][0]=l_eq;

  r1[2][0]=r1[1][0]+l_eq*cos(pi-a_eq);
  r1[2][1]=r1[1][1]+l_eq*sin(pi-a_eq);

  r1[3][0]=r1[2][0]+l_eq;
  r1[3][1]=r1[2][1];

  r1[4][0]=r1[3][0]+l_eq*cos(pi-a_eq);
  r1[4][1]=r1[3][1]+l_eq*sin(pi-a_eq);




  for (i=0;i<2;++i)
    for (j=0;j<5;++j)
      fprintf(inifile,"%12.8lf %12.8lf %12.8lf\n",r1[j][0],r1[j][1],r1[j][2]);

  fclose(inifile);
  
  return 0;
}



