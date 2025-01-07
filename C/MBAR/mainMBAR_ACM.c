
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "netcdf_mine.h"
#include "EF.h"
#include "IO.h"

#include "MBAR.h"

#define ON 1
#define OFF 0

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,t;
  int numarg;
  int flag=OFF,aflag=OFF,bflag=OFF,cflag=OFF,nflag=OFF;
  int interval=1;
  int d=0,num;
  double f;
  double din;
  int n_sim;
  int *n;

  int n_total=0;
  double *W,*WT;

  double criteria_BAR=1.0e-7;
  int MAXITE=1000;

  double *fene;
  double ***enek;
  double *covma,*covmb,*covmc;

  double *c1,*c2;

  double numt=0.0;
  double *fene_n;

  char *progname;
  char *enekfilename,*eneklistfilename,*fenefilename;
  char *outputfilenamea,*outputfilenameb,*outputfilenamec;

  FILE *enekfile,*eneklistfile,*fenefile;
  FILE *outputfile;

  int ncid,fene_varid;
  size_t start[1],count[1];

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"hatxbcmnp:"))!=-1) {
    switch(c) {
    case 'm':
      flag=ON;
      break;
    case 'n':
      nflag=ON;
      break;
    case 'b':
      bflag=ON;
      break;
    case 'c':
      cflag=ON;
      break;
    case 't':
      criteria_BAR=atof(optarg);
      break;
    case 'x':
      MAXITE=atoi(optarg);
      break;
    case 'p':
      interval=atoi(optarg);
      break;
    case 'a':
      aflag=atoi(optarg);
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  if ( nflag==ON && flag==ON ) {
    printf("option error: -m -n");
    exit(1);
  }

  if ( aflag==OFF && bflag==OFF && cflag==OFF ) {
    aflag=ON;
    exit(1);
  }

  argc-=optind;
  argv+=optind;

  numarg=4;
  if (aflag==ON && bflag==ON && cflag==OFF)
    ++numarg;
  if (aflag==ON && bflag==OFF && cflag==ON)
    ++numarg;
  if (aflag==OFF && bflag==ON && cflag==ON)
    ++numarg;
  if (aflag==ON && bflag==ON && cflag==ON)
    numarg+=2;
  if (argc < numarg) {
    USAGE(progname);
    exit(1);
  }
  n_sim = atoi(*argv);
  fenefilename = *++argv;
  eneklistfilename = *++argv;
  if (cflag==ON && aflag==OFF && bflag==OFF)
    outputfilenamea = *++argv;
  if (cflag==OFF && aflag==ON && bflag==OFF)
    outputfilenameb = *++argv;
  if (cflag==OFF && aflag==OFF && bflag==ON)
    outputfilenameb = *++argv;
  if (cflag==ON && aflag==ON && bflag==OFF) {
    outputfilenamea = *++argv; 
    outputfilenamec = *++argv;
 }
  if (cflag==ON && aflag==OFF && bflag==ON) {
    outputfilenamea = *++argv; 
    outputfilenameb = *++argv;
 }
  if (cflag==OFF && aflag==ON && bflag==ON) {
    outputfilenameb = *++argv; 
    outputfilenamec = *++argv;
 }
  if (cflag==ON && aflag==ON && bflag==ON) {
    outputfilenamea = *++argv;
    outputfilenameb = *++argv;
    outputfilenamec = *++argv;
 }

  fene=(double *)gcemalloc(sizeof(double)*n_sim);
  fene_n=(double *)gcemalloc(sizeof(double)*n_sim);
  if (flag==OFF && nflag ==OFF) {
    fenefile=efopen(fenefilename,"r");
    getline(&line,&len,fenefile);
    for (i=0;i<n_sim;++i)
      fscanf(fenefile,"%d %lf",&d,&fene[i]);
    fclose(fenefile);
  }
  if (nflag==ON) {
    start[0]=0;
    count[0]=n_sim;
    enc_open(fenefilename,NC_NOWRITE,&ncid);
    nc_inq_varid(ncid,"fene",&fene_varid);
    nc_get_vara_double(ncid,fene_varid,start,count,fene);
    encclose(ncid);
  }

  n=(int *)gcemalloc(sizeof(int)*n_sim);
  enek=(double ***)gcemalloc(sizeof(double **)*n_sim);
  for (i=0;i<n_sim;++i) enek[i]=(double **)gcemalloc(sizeof(double *)*n_sim);
  for (i=0;i<n_sim;++i) for (j=0;j<n_sim;++j) enek[i][j]=(double *)gcemalloc(sizeof(double));

  eneklistfile=efopen(eneklistfilename,"r");
  for (i=0;i<n_sim;++i) {
    for (j=0;j<n_sim;++j) {
      getline(&line,&len,eneklistfile);
      line[strlen(line)-1]='\0';
      enekfile=efopen(line,"r");
      getline(&line,&len,enekfile);
      k=0;
      num=0;
      d = 1;
      while ( d != -1  )  {
	d=fscanf(enekfile,"%lf",&f);
	d=fscanf(enekfile,"%lf",&f);
	if (k%interval == 0) {
	  enek[j][i]=(double *)gcerealloc(enek[j][i],sizeof(double)*(num+1));
	  enek[j][i][num]=f;
	  ++num;
	}
	++k;
      } 
      fclose(enekfile);
      n[i]=num-1;
    }
  }
  fclose(eneklistfile);

  for (i = 0; i < n_sim; ++i) {
    numt=0.0;
    for (j = 0; j < n_sim; ++j) {
      for (t = 0; t < n[j]; ++t) {
	din=0.0;
	for (k = 0; k < n_sim; ++k) {	  
	  din+=n[k]*exp(fene[k]-enek[k][j][t]);
	}
	numt+=exp(-enek[i][j][t])/din;
      }
    }
    fene_n[i]=-log(numt);
  }

  if (flag==ON) MBAR_ite(fene,enek,n_sim,n,criteria_BAR,MAXITE);

  for (i=0;i<n_sim;++i) n_total+=n[i];

  W=gcemalloc(sizeof(double)*n_total*n_sim); 
  WT=(double *)gcemalloc(sizeof(double)*n_total*n_sim);
  //  MBAR_setw(fene,enek,n_sim,n,W,WT);
  l=-1;
  for (i=0;i<n_sim;++i) {
    for (t=0;t<n[i];++t) {
      ++l;
      for (j=0;j<n_sim;++j) {
  	din=0.0;
  	for (k=0;k<n_sim;++k) {
  	  din+=n[k]*exp(fene[k]-enek[k][i][t]);
  	}
  	W[l*n_sim+j]=exp(fene[j]-enek[j][i][t])/din;
	WT[j*n_total+l] =W[l*n_sim+j];
      }
    }
  }

  if ( cflag==ON ) covmc=MBAR_ACM0(n_sim,n,W,WT);
  if ( aflag==ON ) covma=MBAR_ACM(n_sim,n,W);
  if ( bflag==ON ) covmb=MBAR_ACM2(n_sim,n,W,WT);

  if (flag==ON) {
    fenefile=efopen(fenefilename,"w");
    fprintf(fenefile,"# fene\n");
    for (i=0;i<n_sim;++i)
      fprintf(fenefile,"%d %12.10lf\n",i,fene[i]);
    fclose(fenefile);
  }

  if ( aflag==ON ) {
    outputfile=efopen(outputfilenamea,"w");
    for (i=0;i<n_sim;++i)
      fprintf(outputfile,"%d %e\n",i,covma[i*n_sim+i]-2.0*covma[i*n_sim]+covma[0]);
    fprintf(outputfile,"\n");
    fclose(outputfile);
  }
  if ( bflag==ON ) {
    outputfile=efopen(outputfilenameb,"w");
    for (i=0;i<n_sim;++i)
      fprintf(outputfile,"%d %e\n",i,covmb[i*n_sim+i]-2.0*covmb[i*n_sim]+covmb[0]);
    fprintf(outputfile,"\n");
    fclose(outputfile);
  }
  if ( cflag==ON ) {
    outputfile=efopen(outputfilenamec,"w");
    for (i=0;i<n_sim;++i)
      fprintf(outputfile,"%d %e\n",i,covmc[i*n_sim+i]-2.0*covmc[i*n_sim]+covmc[0]);
    fprintf(outputfile,"\n");
    fclose(outputfile);
  }

  return 0;
}

void USAGE(char *progname) {
  printf("[-m] -- with MBAR calculation\n");
  printf("[-t criteria_BAR ] -- criteria of MBAR iteration (default 10-7)\n");
  printf("[-x max_ite_step_of_BAR ] -- max num of MBAR iteration\n");
  printf("[-n] -- netcdf input\n");
  printf("[-a] -- calculation of error will be done by the method avoiding GIM calculation\n");
  printf("[-b] -- calculation of error will be done by the both methods\n");
  printf("[-c] -- calculation of error will be done by the three methods\n");
  printf("[-p num] -- interval of reading data\n");
  printf("[-h] -- help\n");
  printf("USAGE:%s [-m] [-n] [-a] [-b] [-c] [-p num ] [-h] n_sim fenefilename enelistfilename outputfilename1 ( outputfilename2 outputfilename3 )  \n", progname);
}
