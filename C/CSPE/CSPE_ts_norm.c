
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>

#include "fftw3.h"

#include "IO.h"
#include "EF.h"
#include "SPE.h"
#include "mymath.h"
#include "netcdf_mine.h"

#define ON 1
#define OFF 0

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,n;
  double f;

  int vflag=OFF,normflag=ON;

  int numstep,num;
  int num_of_split=1,numini=1,numfin,numinirow,numfinrow,numrow;

  double c=2.999792e-2;
  double kb=1.98723e-3*4.18407*100.0;

  double dt,pi,sum=0.0;
  double **spe,**dat,**dat_norm,**speave;

  char *inputfilename,*outputfilename,*outputfilename2,*progname;
  FILE *inputfile,*outputfile,*outputfile2;

  int ncid,indexofsim_dimid,fene_varid;
  size_t start[1],count[1];

  char *line,*dummy;
  size_t len=0;

  int c2;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c2=getopt(argc,argv,"hovi:f:j:g:d:s:r:"))!=-1) {
    switch(c2) {
    case 'v':
      vflag=ON;
      break;
    case 'i':
      numini=atof(optarg)-1;
      break;
    case 'f':
      numfin=atoi(optarg)-1;
      break;
    case 'j':
      numinirow=atoi(optarg)-1;
      break;
    case 'g':
      numfinrow=atoi(optarg)-1;
      break;
    case 'd':
      dt=atof(optarg);
      break;
    case 's':
      num_of_split=atoi(optarg);
      break;
    case 'r':
      numrow=atoi(optarg);
      break;
    case 'o':
      normflag=OFF;
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  argc-=optind;
  argv+=optind;

  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  inputfilename  = *argv;
  outputfilename = *++argv;
  outputfilename2 = *++argv;

  inputfile=efopen(inputfilename,"r");
  for (i=0;i<numini;++i)  getline(&line,&len,inputfile);

  numstep=(int)((numfin-numini+1)/num_of_split);
  num=numfinrow-numinirow+1;
  if (num<1) {
    printf("error: numfinrow must be larger or equal to numinirow\n");
    exit(1);
  }

  speave=(double **)gcemalloc(sizeof(double *)*num);
  for (i=0;i<num;++i)
    speave[i]=(double *)gcemalloc(sizeof(double)*numstep);
  for (i=0;i<num;++i)
    for (j=0;j<numstep;++j)
      speave[i][j]=0.0;

  for (i=0;i<num_of_split;++i) {
    spe=(double **)gcemalloc(sizeof(double)*num);
    dat=(double **)gcemalloc(sizeof(double)*num);
    dat_norm=(double **)gcemalloc(sizeof(double)*num);
    for (j=0;j<num;++j) {
      spe[j]=(double *)gcemalloc(sizeof(double)*numstep);
      dat[j]=(double *)gcemalloc(sizeof(double)*numstep);
      dat_norm[j]=(double *)gcemalloc(sizeof(double)*numstep);
    }

    for (j=0;j<numstep;++j) {
      l=0;
      for (k=0;k<numrow;++k) {
	if (k>=numinirow && k<=numfinrow ) {
	  fscanf(inputfile,"%lf",&dat[l][j]);
	  ++l;
	}
	else {
	  fscanf(inputfile,"%lf",&f);
	}
      }
    }

    for (j=0;j<num;++j) {
      if (normflag==ON) {
	dat_norm[j]=ts_normalize_od(dat[j],numstep);
	FT_ts(numstep,dat_norm[j],spe[j]);
      }
      else {
	FT_ts(numstep,dat[j],spe[j]);
      }
      for (k=0;k<numstep;++k) {
	speave[j][k]=(i*speave[j][k]+spe[j][k])/(i+1);
      }
    }
  }
  fclose(inputfile);
  
  pi=acos(-1.0);
  for (i=0;i<num;++i) {
    if (vflag==OFF) 
      for (j=0;j<numstep;++j) 
	speave[i][j]=(2.0*pi*j/numstep/dt)*(2.0*pi*j/numstep/dt)*speave[i][j];
    sum=0.0;
    for (j=0;j<(int)(numstep/2);++j) sum+=speave[i][j];
    for (j=0;j<numstep;++j) speave[i][j]=speave[i][j]/sum;  
  }

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    fprintf(outputfile,"%e %e ",(double)i/numstep/dt/c,(double)numstep*dt/i);
    for (j=0;j<num;++j) {
      fprintf(outputfile,"%e ",speave[j][i]);
    }
    fprintf(outputfile,"\n");
  }
  fclose(outputfile);

  outputfile2=efopen(outputfilename2,"w");
  for (i=0;i<numstep;++i) {
    for (j=0;j<num;++j) {
      fprintf(outputfile2,"%e ",dat_norm[j][i]);
    }
    fprintf(outputfile2,"\n");
  }
  fclose(outputfile2);

  return 0;
}

void USAGE(char *progname) {
  printf("-v -- derivative mode  \n");
  printf("-i -- numini\n");
  printf("-f -- numfinp\n");
  printf("-j -- numinirow\n");
  printf("-g -- numfinrow\n");
  printf("-d -- deltat\n");
  printf("-s -- num of split s\n");
  printf("-r -- numrow\n");
  printf("-h -- help\n");
  printf("USAGE:%s [-n] [-v] [-i] [-j] [-g] [-d] [-s] [-r] [-h] inputfilename outputfile\n", progname);
}
