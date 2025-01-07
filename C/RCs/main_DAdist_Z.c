
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "netcdf_mineL.h"

#include "TOPO.h"
#include "LA.h"
#include "mymath.h"

#include "TACCM.h"
#include "PTL.h"
#include "MB.h"
#include "EF.h"

#define ON 0
#define OFF 1

double CDon(double atom[4][3] , double pi);
int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,m,n,nr=1,d;
  int numref=2,numatom,*numdihed,numdihedall,numdihedall2;
  int numstep=1;
  int amberncflag=OFF,flag=1;
  double Theta,*theta,dtheta,a;
  int **pairs;

  double *deltad,pi;
  double **refcrd,**refd;

  double atom[4][3];
  double crd_nc[MAXATOM][3],*crd;
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id_MD;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*inputfilename2,**refcrdfilename,*TACCMfilename,*TACCMfilename2,*parmfilename,*outputfilename;

  FILE *inputfile,*inputfile2,**refcrdfile,*TACCMfile,*TACCMfile2,*parmfile,*outputfile,*logfile;

  char *progname;
  int opt_idx=1;

  pi=acos(-1.0);

  struct option long_opt[] = {
    {"nref",1,NULL,'n'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hn:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'A':
      amberncflag=ON;
      break;
    case 'n':
      numref=atoi(optarg); break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);   exit(1);
    }
  }

  deltad=(double *)gcemalloc(sizeof(double)*numref);
  refcrd=(double **)gcemalloc(sizeof(double *)*numref);
  refd=(double **)gcemalloc(sizeof(double *)*numref);
  refcrdfilename=(char **)gcemalloc(sizeof(char *)*numref);
  for (i=0;i<numref;++i) refcrdfilename[i]=(char *)gcemalloc(sizeof(char )*1000);
  refcrdfile=(FILE **)gcemalloc(sizeof(FILE *)*numref);

  progname=*argv;  argc-=optind;  argv+=optind;

  if (argc < 6+numref) {
    USAGE(progname);
    exit(1);
  }
  inputfilename    = *argv;
  inputfilename2   = *++argv;
  for (i=0;i<numref;++i) refcrdfilename[i] = *++argv;
  TACCMfilename   = *++argv;
  TACCMfilename2  = *++argv;
  parmfilename     = *++argv;
  outputfilename   = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;

  TACCMfile=efopen(TACCMfilename,"r");
  fscanf(TACCMfile,"%d",&numdihedall);
  pairs=(int **)gcemalloc(sizeof(int *)*numdihedall);
  for (i=0;i<numdihedall;++i) pairs[i]=(int *)gcemalloc(sizeof(int)*5);
  for (i=0;i<numdihedall;++i) {
    for (j=0;j<4;++j) fscanf(TACCMfile,"%d",&(pairs[i][j]));
    fscanf(TACCMfile,"%d",&(pairs[i][j]));
  }
  fclose(TACCMfile);

  TACCMfile2=efopen(TACCMfilename2,"r");
  fscanf(TACCMfile2,"%d",&numdihedall2);
  fclose(TACCMfile2);

  refcrd=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numref;++i) {
    refcrdfile[i]=efopen(refcrdfilename[i],"r");
    refcrd[i]=(double *)gcemalloc(sizeof(double)*numatom*3);
    getline(&line,&len,refcrdfile[i]);
    fscanf(refcrdfile[i],"%d",&d);
    for (j=0;j<numatom;++j) for (k=0;k<3;++k) fscanf(refcrdfile[i],"%lf",&refcrd[i][j*3+k]);
    fclose(refcrdfile[i]);
  }

  //  numdihedall=0;for (i=0;i<flag;++i) numdihedall+=numdihed[i];

  theta=(double *)gcemalloc(sizeof(double)*numdihedall);
  for (i=0;i<numref;++i) refd[i]=(double *)gcemalloc(sizeof(double)*numdihedall);
  
  for (n=0;n<numref;++n) {
    m=0;
    for (i=0;i<numdihedall;++i) {
      for (j=0;j<4;++j) for (k=0;k<3;++k) atom[j][k]=refcrd[n][(pairs[i][j]-1)*3+k];
    
      Theta=CDon(atom,pi);
      if (Theta > 2.0*pi)  Theta = Theta -2.0*pi;
      else if (Theta < 0.0)  Theta = Theta +2.0*pi;

      refd[n][m]=Theta;
      ++m;
    }
  }

  if (amberncflag==OFF) numstep=myncL_get_present_step_MCD(inputfilename,&nc_id_MCD);
  else numstep=myncL_get_present_step_AMBER(inputfilename,&nc_id_MD);

  inputfile2=efopen(inputfilename2,"r");
  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    for (j=0;j<numdihedall2;++j) {
      fscanf(inputfile2,"%lf",&theta[j]);
    }
    //    if (amberncflag==OFF) myncL_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd_nc);
    //    else myncL_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id_MD,crd_nc);

    for (nr=0;nr<numref;++nr) {
      deltad[nr]=0.0;
      for (j=0;j<numdihedall;++j) {
	if (theta[j] > 2.0*pi)  theta[j] = theta[j] -2.0*pi;
	else if (theta[j] < 0.0) theta[j] = theta[j] +2.0*pi;

	dtheta=theta[j]-refd[nr][j];
	if (dtheta < 0.0) dtheta=-dtheta;

	a=2.0*pi-dtheta;

	if (dtheta>a) dtheta=a;

	deltad[nr]+=dtheta;
	++n;
      }
      fprintf(outputfile,"%8.4lf  ",deltad[nr]/pi/numdihedall);
    }
    fprintf(outputfile,"\n");
  }
  fclose(inputfile2);
  fclose(outputfile);

}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-A] amberncflag \n");
  printf("[-P] phi psi \n");
  printf("[-O] phi psi omega \n");
  printf("[-K] phi psi omega kai \n");
  printf("[--nref numref] numref \n");
  printf("[-h ] help\n");
  printf("%s inputfilename inputfilename2 refcrdfilename parmfilename outputfilename \n",progname);
}

double CDon(double atom[4][3], double pi) {
  int i,j,k,l;
  int ii,jj,kk,ll;

  double m[3],n[3],m_n[3],n_n[3],lm,ln;
  double vij[3],vkj[3],vkl[3];
  double lkj;
  double vijvkj,vklvkj;

  double dihed,theta;

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
  
  theta=dihed;
    
  return theta;
}
