#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

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
  int numatom,numstep,numdihed,numcrd=0,numref;

  int flagout,amberncflag=OFF;
  int **pairs;

  int interval=1;

  double *DA_target,width;

  double *deltad,pi;
  double **refcrd,**refd;

  double atom[4][3];
  double *theta,Thetao,dtheta,a;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc[MAXATOM][3],*crd;
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id_MD;

  char *inputfilename,**refcrdfilename,*parmfilename,*condfilename;
  char *crdfilenamebase,crdfilename[100],*logfilename;
  char *progname;
  FILE *inputfile,**refcrdfile,*parmfile,*condfile,*crdfile,*logfile;

  int opt_idx=1;

  struct option long_opt[] = {
    {"h",0,NULL,'h'},
    {"N",0,NULL,'N'},
    {"K",0,NULL,'K'},
    {"nref",1,NULL,'n'},
    {"int",1,NULL,'v'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hNk:n:v:",long_opt,&opt_idx))!=-1) {
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
    case 'n':
      numref=atoi(optarg); break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  deltad=(double *)gcemalloc(sizeof(double)*numref);
  refcrd=(double **)gcemalloc(sizeof(double *)*numref);
  refd=(double **)gcemalloc(sizeof(double *)*numref);
  refcrdfilename=(char **)gcemalloc(sizeof(char *)*numref);
  for (i=0;i<numref;++i) refcrdfilename[i]=(char *)gcemalloc(sizeof(char )*1000);
  refcrdfile=(FILE **)gcemalloc(sizeof(FILE *)*numref);

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 6+numref) {
    USAGE(progname);
    exit(1);
  }
  inputfilename   = *argv;
  for (i=0;i<numref;++i) refcrdfilename[i] = *++argv;
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

  refcrd=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numref;++i) {
    refcrdfile[i]=efopen(refcrdfilename[i],"r");
    refcrd[i]=(double *)gcemalloc(sizeof(double)*numatom*3);
    getline(&line,&len,refcrdfile[i]);
    fscanf(refcrdfile[i],"%d",&d);
    for (j=0;j<numatom;++j) for (k=0;k<3;++k) fscanf(refcrdfile[i],"%lf",&refcrd[i][j*3+k]);
    fclose(refcrdfile[i]);
  }

  condfile=efopen(condfilename,"r");

  fscanf(condfile,"%d",&numdihed);
  pairs=(int **)gcemalloc(sizeof(int *)*numdihed);
  DA_target=(double *)gcemalloc(sizeof(double)*numref);

  for (i=0;i<numdihed;++i) {
    pairs[i]=(int *)gcemalloc(sizeof(int)*4);
    for (j=0;j<4;++j) {
      fscanf(condfile,"%d",&pairs[i][j]);
      pairs[i][j]-=1;
    }
  }
  for (i=0;i<numref;++i) {
    fscanf(condfile,"%lf",&DA_target[i]);
    while (DA_target[i] > pi) DA_target[i]-=2.0*pi;
    while (DA_target[i] < -1.0*pi) DA_target[i]+=2.0*pi;
  }
  fclose(condfile);

  for (i=0;i<numref;++i) refd[i]=(double *)gcemalloc(sizeof(double)*numdihed);
  
  for (i=0;i<numref;++i) {
    for (j=0;j<numdihed;++j) {
      for (k=0;k<4;++k) 
	for (l=0;l<3;++l) 
	  atom[k][l]=refcrd[i][(pairs[j][k])*3+l];
    
      CDon(atom,pi,&Thetao);
      if (Thetao > 2.0*pi)    Thetao = Thetao -2.0*pi;
      else if (Thetao < 0.0)  Thetao = Thetao +2.0*pi;
      
      refd[i][j]=Thetao;
    }
  }

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
      for (j=0;j<numdihed;++j) {
	for (k=0;k<4;++k) 
	  for (l=0;l<3;++l)
	    atom[k][l]=crd_nc[(pairs[j][k])][l];

	CDon(atom,pi,&theta[j]);

	if (theta[j] > 2.0*pi)  theta[j] = theta[j] -2.0*pi;
	else if (theta[j] < 0.0) theta[j] = theta[j] +2.0*pi;
      }

      flagout=ON;
      for (j=0;j<numref;++j) {
	for (k=0;k<numdihed;++k) {
	  dtheta=theta[k]-refd[j][k];
	  if (dtheta < 0.0) dtheta=-dtheta;
	  a=2.0*pi-dtheta;
	  if (dtheta>a) dtheta=a;

	  deltad[j]+=dtheta;
	}
	deltad[j]=deltad[j]/pi/numdihed;
	if (deltad[j] < DA_target[j]-width || deltad[j] >= DA_target[j]+width) {
	  flagout=OFF;
	}
      }

      //      for (j=0;j<numref;++j) fprintf(logfile,"%10.8lf ",deltad[j]);
      //      fprintf(logfile,"\n");

      if (flagout==ON) {
	for (j=0;j<numref;++j) fprintf(logfile,"%10.8lf ",deltad[j]);
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
