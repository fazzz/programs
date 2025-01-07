
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <math.h>

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j;

  int N_bin=20;

  double prob;
  double x,y;
  double dx,dy;
  double minx=-3.0,maxx=3.0,miny=-3.0,maxy=3.0;
  double x0,y0;
  double hx,hy;

  char *filename;
  FILE *file;

  if (argc < 5) {
    printf("Error\n");
    exit(1);
  }  
  x0 = atof(*++argv);
  y0 = atof(*++argv);
  hx = atof(*++argv);
  hy = atof(*++argv);
  filename = *++argv;

  file=fopen(filename,"w");
  for (i=0;i<N_bin;++i) {
    x=(maxx-minx)/(double)N_bin*i+minx;
    for (j=0;j<N_bin;++j) {
      y=(maxy-miny)/(double)N_bin*j+miny;

      dx=x-x0;
      dy=y-y0;

      prob=exp(-0.5*dx*dx/hx)*exp(-0.5*dy*dy/hy);

      fprintf(file,"%10.8lf %10.8lf %10.8lf\n",x,y,prob);
    }
    fprintf(file,"\n");
  }

  fclose(file);
}
