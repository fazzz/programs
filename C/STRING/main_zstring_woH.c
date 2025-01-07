
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "STRING.h"
#include "CSI.h"
#include "PT.h"
#include "FF.h"
#include "EF.h"
#include "IO.h"

#define ON 1
#define OFF 0

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int flag=ON,flagB=ON,flagA=ON,flagD=ON,flagL=ON,flagE=ON;

  int sl=0;

  int outinterval=1;

  int numpoint,numiteration=10000;
  int numatom,numpara,numnb,num14;

  double dt=0.01;

  double scpLJ=1.0,scpE=1.0;

  double maxdelta,*pe_old;
  double criteria=1.0e-3;

  char *line;
  size_t len=0;

  char *progname;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double *path,*path_evoluted;
  double *fe,*pe;

  struct potential e;
  struct force f;
  double *crd;

  char *inifilename,*parmtopname,*outputfilename,*outputfilename2;
  FILE *inifile,*parmtop,*outputfile,*outputfile2,*logfile1,*logfile2;

  progname=argv[0];

  while((c=getopt(argc,argv,"BADLEl:e:hs:t:o:c:m:"))!=-1) {
    switch(c) {
    case 'B':
      flagB=OFF;
      break;
    case 'A':
      flagA=OFF;
      break;
    case 'D':
      flagD=OFF;
      break;
    case 'L':
      flagL=OFF;
      break;
    case 'E':
      flagE=OFF;
      break;
    case 'l':
      scpLJ=atof(optarg);
      break;
    case 'e':
      scpE=atof(optarg);
      break;
    case 's':
      sl=atoi(optarg);
      break;
    case 't':
      dt=atof(optarg);
      break;
    case 'o':
      outinterval=atoi(optarg);
      break;
    case 'c':
      criteria=atof(optarg);
      break;
    case 'm':
      numiteration=atoi(optarg);
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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  numpoint=atoi(*argv);
  inifilename  = *++argv;
  parmtopname  = *++argv;
  outputfilename=*++argv;
  outputfilename2=*++argv;

  parmtop=efopen(parmtopname,"r");
  readParmtop(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;

  path=(double *)gcemalloc(sizeof(double)*numatom*3*numpoint);
  inifile=efopen(inifilename,"r");
  for (i=0;i<sl;++i) getline(&line,&len,inifile);
  io_scantrj2(inifile,numatom,numpoint,path);
  fclose(inifile);

  path_evoluted=(double *)gcemalloc(sizeof(double)*numatom*3*numpoint);
  fe=(double *)gcemalloc(sizeof(double)*numatom*3*numpoint);
  pe=(double *)gcemalloc(sizeof(double)*numpoint);
  pe_old=(double *)gcemalloc(sizeof(double)*numpoint);
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  ff_set_calcffandforce(&e,&f);

  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");
  logfile1=efopen("log1_zswoH","w");
  logfile2=efopen("log2_zswoH","w");
  for (i=0;i<numiteration;++i) {
    for (j=0;j<numpoint*numatom*3;++j) fe[j]=0.0;
    for (j=0;j<numpoint;++j) pe[j]=0.0;

    for (j=0;j<numpoint;++j) {
      for (k=0;k<numatom;++k) {
	if (strncmp(AP.IGRAPH[k],"H",1)!=0 ) {
 	  for (l=0;l<3;++l) {
	    crd[k*3+l]=path[j*numatom*3+k*3+l];
	  }
	}
	else {
	  for (l=0;l<3;++l) {
	    crd[k*3+l]=0.0;
	    path[j*numatom*3+k*3+l]=0.0;
	  }
	}
      }

      ff_calcffandforce_woH(crd,numatom,&e,&f);
	if (flagB==ON) {
	  for (k=0;k<numatom*3;++k) fe[j*numatom*3+k]+=-f.f_b[k];
	  pe[j]+=e.p_b_t;
	}
	if (flagA==ON) {
	  for (k=0;k<numatom*3;++k) fe[j*numatom*3+k]+=/*-*/f.f_a[k];
	  pe[j]+=e.p_a_t;
	}
	if (flagD==ON) {
	  for (k=0;k<numatom*3;++k) fe[j*numatom*3+k]+=f.f_d[k];
	  pe[j]+=e.p_d_t;
	}
	if (flagL==ON) {
	  for (k=0;k<numatom*3;++k)  fe[j*numatom*3+k]+=scpLJ*(f.f_LJ[k]+f.f_LJ_14[k]);
	  pe[j]+=scpLJ*(e.p_LJ_t+e.p_LJ_14_t);
	}
	if (flagE==ON) {
	  for (k=0;k<numatom*3;++k)  fe[j*numatom*3+k]+=scpE*(f.f_e[k]+f.f_e_14[k]);
	  pe[j]+=scpE*(e.p_e_t+e.p_e_14_t);
	}
    }
    z_string_Cartesian(path,path_evoluted,fe,numpoint,numatom,dt);
    printf("%d %lf %lf %lf\n",i,pe[0],pe[(int)(numpoint/2)],pe[numpoint-1]);
    if (i>0) {
      maxdelta=fabs(pe[0]-pe_old[0]);
      for (j=1;j<numpoint;++j) {
	if (maxdelta<abs(pe[j]-pe_old[j])) maxdelta=pe[j]-pe_old[j];
      }
      printf("%lf\n",fabs(maxdelta));
      if (fabs(maxdelta)<criteria) break;
    }
    for (j=0;j<numpoint;++j) pe_old[j]=pe[j];
    if (i%outinterval==0) {
      io_outtrj2(logfile1,numatom,numpoint,path);
      for (j=0;j<numpoint;++j) fprintf(logfile2,"%4d %12.8e \n",j,pe[j]);
    }
  }
  fclose(logfile1);
  fclose(logfile2);
  io_outtrj2(outputfile,numatom,numpoint,path);
  fprintf(outputfile2,"# ene\n");
  for (i=0;i<numpoint;++i) fprintf(outputfile2,"%4d %12.8lf\n",i,pe[i]);
  fclose(outputfile);
  fclose(outputfile2);
}

int USAGE(char *progname) {
  printf("USAGE:%s \n",progname);
  printf("[-B] BOND OFF [-A] ANGLE OFF [-D] DIHED OFF [-L] L-J [-E] elesta \n");
  printf("[-l scl -e sce ]   \n");
  printf("[-s sl] (dissmiss lines)    \n");
  printf("[-t dt]    \n");
  printf("[-o outinterval] \n");
  printf("[-c criteria] \n");
  printf("[-m numiteration] \n");
  printf("[-h] help  \n");
  printf("%s numpoint inifilename(path) parmtopname outputfilename(path) outputfilename2(ene)\n",progname);
}

 
