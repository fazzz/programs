#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "netcdf_mineL.h"

#include "TOPO.h"
#include "LA.h"
#include "mymath.h"

#include "MB.h"
#include "PTL.h"
#include "EF.h"
#include "IO.h"

#define ON 1
#define OFF 0

double CDon(double atom[4][3], double pi, double *theta);

int usage(void);

int main(int argc, char *argv[]) {
  int i,j,k,l,d,dummy,initialstep=0;
  int numatom,numstep,numdihed,numcrd=0;

  int flagout,amberncflag=OFF;
  int **pairs;

  int interval=1;

  double *dihed_target,width;

  double atom[4][3];
  double *theta;
  double pi;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc[MAXATOM][3],*crd;
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id_MD;

  char *inputfilename,*parmfilename,*condfilename;
  char *crdfilenamebase,crdfilename[100],*logfilename;
  char *progname;
  FILE *inputfile,*parmfile,*condfile,*crdfile,*logfile;

  while((c=getopt(argc,argv,"hNk:v:"))!=-1) {
    switch(c) {
    case 'N':
      amberncflag=ON;
      break;
    case 'v':
      interval=atoi(optarg);
      break;
    case 'k':
      initialstep=atoi(optarg);
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

  if (argc < 6) {
    USAGE(progname);
    exit(1);
  }
  inputfilename   = *argv;
  parmfilename    = *++argv;
  condfilename    = *++argv;
  width           = atof(*++argv);
  crdfilenamebase = *++argv;
  logfilename     = *++argv;

  if (amberncflag==OFF) numstep=myncL_get_present_step_MCD(inputfilename,&nc_id_MCD);
  else numstep=myncL_get_present_step_AMBER(inputfilename,&nc_id_MD);

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  pi=acos(-1.0);

  condfile=efopen(condfilename,"r");

  fscanf(condfile,"%d",&numdihed);
  pairs=(int **)gcemalloc(sizeof(int *)*numdihed);
  dihed_target=(double *)gcemalloc(sizeof(double)*numdihed);

  for (i=0;i<numdihed;++i) {
    pairs[i]=(int *)gcemalloc(sizeof(int)*4);
    for (j=0;j<4;++j) {
      fscanf(condfile,"%d",&pairs[i][j]);
      pairs[i][j]-=1;
    }
    fscanf(condfile,"%lf",&dihed_target[i]);
    while (dihed_target[i] > pi) dihed_target[i]-=2.0*pi;
    while (dihed_target[i] < -1.0*pi) dihed_target[i]+=2.0*pi;
  }
  fclose(condfile);

  theta=(double *)gcemalloc(sizeof(double)*numdihed);

  if (initialstep>0) {
    for (i=0;i<initialstep;++i) {
      if (amberncflag==OFF) myncL_open_inq_get_sh_MCD(inputfilename,numatom,0,1,1,&nc_id_MCD,crd_nc);
      else myncL_open_inq_get_sh_AMBER(inputfilename,numatom,0,1,1,&nc_id_MD,crd_nc);
    }
  }

  logfile=efopen(logfilename,"w");

  for (i=initialstep;i<numstep;++i) {
    if (amberncflag==OFF) myncL_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd_nc);
    else myncL_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id_MD,crd_nc);

    if ( (i%interval) == 0 || i==initialstep ) {
      flagout=ON;
      for (j=0;j<numdihed;++j) {
	for (k=0;k<4;++k) 
	  for (l=0;l<3;++l)
	    atom[k][l]=crd_nc[(pairs[j][k])][l];

	CDon(atom,pi,&theta[j]);

	if (theta[j] < dihed_target[j]-width || theta[j] >= dihed_target[j]+width) {
	  flagout=OFF;
	}
      }

      if (flagout==ON) {
	for (j=0;j<numdihed;++j) fprintf(logfile,"%10.8lf ",theta[j]);
	fprintf(logfile,"\n");

	++numcrd;
	sprintf(crdfilename,"%s_%d",crdfilenamebase,numcrd);
	crdfile=efopen(crdfilename,"w");
	for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd[j*3+k]=crd_nc[j][k];

	io_outputconf_Amberform(crdfile,numatom,crd);
	fclose(crdfile);
      }
    }
  }
  fclose(logfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-N] amberncflag \n");
  printf("[-k interval] specify initial \n");
  printf("%s inputfilename outputfilename \n",progname);
}

double CDon(double atom[4][3], double pi, double *theta) {
  int i,j,k,l;
  int ii,jj,kk,ll;

  double m[3],n[3],m_n[3],n_n[3],lm,ln;
  double vij[3],vkj[3],vkl[3];
  double lkj;
  double vijvkj,vklvkj;

  double dihed;

  for (j=0;j<3;++j) {
    vij[j] = atom[1][j]-atom[0][j];
    vkj[j] = atom[1][j]-atom[2][j];
    vkl[j] = atom[3][j]-atom[2][j];
  }
  lkj=sqrt(inprod(vkj,vkj,3));
  
  outprod(vij,vkj,m);
  outprod(vkj,vkl,n);
  lm=sqrt(inprod(m,m,3));
  ln=sqrt(inprod(n,n,3));
  for (j=0;j<3;++j) {
    m_n[j]=m[j]/lm;
    n_n[j]=n[j]/ln;
  }
  
  dihed=inprod(m_n,n_n,3);
  if (dihed>=1.0)
    dihed=0.0;
  else if (dihed<=-1.0)
      dihed=pi;
  else
    dihed=acos(dihed);
  if (inprod(vij,n,3)>0) dihed=-dihed;
  if (dihed<-1.0*pi) dihed=2.0*pi+dihed;
  if (dihed>pi) dihed=-2.0*pi+dihed;
  
  *theta=dihed;
    
  return *theta;
}
