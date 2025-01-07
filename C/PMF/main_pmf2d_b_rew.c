
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PMF.h"
#include "PT.h"
#include "EF.h"

#define ON 1
#define OFF 0

double *io_scandcoldata2(FILE *inputfile,int numi,int numcol,int xcol,int ycol,int *numstep,double *data);

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j;

  int minmaxflag=OFF;
  int udflag=OFF;

  int numstep;
  double *data,*w;
  int numi=0,numcol=2,xcol=1,ycol=2;
  double maxx=1.0,maxy=1.0,minx=0.0,miny=0.0;
  int framex,framey;
  double width;
  double *pmf,pi,pmf_min=0.0,pmf_max=0.0;

  double T=0.0,beta=1.0;
  double k_B=1.98723e-3;
  double UNITT=418.4070;
  
  char *inputfilename,*inputwfilename,*outputfilename,*c;
  FILE *inputfile1,*inputwfile,*outputfile;

  char *line;
  size_t len=0;

  int d;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *progname,*logfilename;

  FILE *logfile;

  int opt_idx=1;

  struct option long_opt[] = {
    {"width",1,NULL,'w'},
    {"u",1,NULL,'u'},
    {"temp",1,NULL,'t'},
    {"beta",1,NULL,'b'},
    {"numi",1,NULL,'i'},
    {"numcol",1,NULL,'c'},
    {"xcol",1,NULL,'x'},
    {"ycol",1,NULL,'y'},
    {"maxx",1,NULL,'A'},
    {"minx",1,NULL,'I'},
    {"maxy",1,NULL,'@'},
    {"miny",1,NULL,'j'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((d=getopt_long(argc,argv,"huw:i:t:b:c:x:y:A:I:@:j:",long_opt,&opt_idx))!=-1) {
    switch(d) {
    case 'w':
      width=atof(optarg);
      break;
    case 'u':
      udflag=ON;
      break;
    case 't':
      T=atof(optarg);
      break;
    case 'b':
      beta=atof(optarg);
      break;
    case 'i':
      numi=atoi(optarg);
      break;
    case 'x':
      xcol=atoi(optarg);
      break;
    case 'y':
      ycol=atoi(optarg);
      break;
    case 'A':
      minmaxflag=ON;
      maxx=atof(optarg);
      break;
    case 'I':
      minmaxflag=ON;
      minx=atof(optarg);
      break;
    case '@':
      minmaxflag=ON;
      maxy=atof(optarg);
      break;
    case 'j':
      minmaxflag=ON;
      miny=atof(optarg);
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  progname=*argv;

  argc-=optind;
  argv+=optind;

  pi=acos(-1.0);

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }  
  inputfilename  = *argv;
  inputwfilename = *++argv;
  outputfilename = *++argv;

  numstep=0;
  data=(double *)gcemalloc(sizeof(double)*2);
  inputfile1=efopen(inputfilename,"r");
  data=io_scandcoldata2(inputfile1,numi,numcol,xcol,ycol,&numstep,data);
  fclose(inputfile1);

  w=(double *)gcemalloc(sizeof(double)*numstep);
  inputwfile=efopen(inputwfilename,"r");
  for (i=0;i<numstep;++i) {
    fscanf(inputwfile,"%lf",&(w[i]));
  }
  fclose(inputwfile);

  if (minmaxflag==OFF)
    pmf=pmf_2dmap_rew(data,w,numstep,width,&maxx,&maxy,&minx,&miny,&framex,&framey,T);
  else {
    pmf=pmf_2dmap_wmaxmin_rew(data,w,numstep,width,width,maxx,maxy,minx,miny,&framex,&framey,T);

    for (i=0;i<=framex;++i) {
      for (j=0;j<framey;++j) {
	if ( pmf_min == 0.0 && pmf[i*(framey)+j] > 0.0 )
	  pmf_min=pmf[i*(framey)+j];
	if ( pmf_min > pmf[i*(framey)+j] && pmf[i*(framey)+j] > 0.0 )
	  pmf_min=pmf[i*(framey)+j];
	if ( pmf_max == 0.0 && pmf[i*(framey)+j] > 0.0 )
	  pmf_max=pmf[i*(framey)+j];
	if ( pmf_max < pmf[i*(framey)+j] && pmf[i*(framey)+j] > 0.0 )
	  pmf_max=pmf[i*(framey)+j];
      }
    }

    for (i=0;i<=framex;++i) {
      for (j=0;j<framey;++j) {
	if (pmf[i*(framey)+j]!=0.0) {
	  if (udflag==OFF)
	    pmf[i*(framey)+j]=-log(pmf[i*(framey)+j])+log(pmf_min);
	  else
	    pmf[i*(framey)+j]=-log(pmf[i*(framey)+j])+log(pmf_max);
	}
	else 
	  if (udflag==ON)
	    pmf[i*(framey)+j]=-1.0;
      }
    }
  }
  if (T>0) beta=1.0/(k_B*T);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<framex;++i) {
    for (j=0;j<framey;++j)
      fprintf(outputfile,"%e %e %e\n",width*i+minx,width*j+miny,/*1.0/beta**/pmf[i*framey+j]);
    fprintf(outputfile,"\n");
  }
  fclose(outputfile);
  
  return 0;
}

double *io_scandcoldata2(FILE *inputfile,int numi,int numcol,int xcol,int ycol,int *numstep,double *data){
  int i,j,k;
  double f;

  char *line;
  size_t len=0;

  for (i=0;i<numi;++i)
    getline(&line,&len,inputfile);

  for (i=(*numstep);;++i) {
    for (j=0;j<numcol;++j) {
      if (fscanf(inputfile,"%lf",&f)!=-1) {
	if (j==xcol-1) {
	  data=(double *)gcerealloc(data,sizeof(double)*(i+1)*2);
	  data[i*2]=f;
	}
	else if (j==ycol-1) {
	  data[i*2+1]=f;
	}
      }
      else {
	*numstep=i-1;
	return data;
      }
    }
  }
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename inputwfilename  outputfilename\n",progname);
}
