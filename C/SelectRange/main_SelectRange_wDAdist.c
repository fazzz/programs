
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "netcdf_mineL.h"

#include "PTL.h"
#include "EF.h"

#define ON 0
#define OFF 1

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int numnc=1,numatom;
  int numstep=1;
  int amberncflag=OFF;

  double crd_nc[MAXATOM][3],*crd;
  struct my_netcdf_out_id_MCD *nc_id_MCD_in,nc_id_MCD_out;
  struct my_netcdf_out_id_AMBER *nc_id_MD_in,nc_id_MD_out;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*parmfilename,*outputfilename;
  FILE *parmfile,*logfile;

  char *progname;
  int opt_idx=1;

  pi=acos(-1.0);

  struct option long_opt[] = {
    {"A",0,NULL,'A'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hA",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'A':
      amberncflag=ON;
      break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);   exit(1);
    }
  }

  progname=*argv;  argc-=optind;  argv+=optind;

  if (argc < 7) {
    USAGE(progname);
    exit(1);
  }
  maxx             =  aotf(*argv);
  minx             =  aotf(*++argv);
  maxy             =  aotf(*++argv);
  miny             =  aotf(*++argv);
  inputfilename    = *++argv;
  parmfilename     = *++argv;
  outputfilename   = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;

  if (amberncflag==OFF) numstep=myncL_get_present_step_MCD(inputfilename,&nc_id_MCD);
  else numstep=myncL_get_present_step_AMBER(inputfilename,&nc_id_MD);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    if (amberncflag==OFF) myncL_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd_nc);
    else myncL_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id_MD,crd_nc);

    for (nr=0;nr<numref;++nr) {
      deltad[nr]=0.0;
      n=0;
      for (j=0;j<flag;++j) {
	for (k=0;k<numdihed[j];++k) {
	  for (m=0;m<4;++m) 
	    for (l=0;l<3;++l)
	      atom[m][l]=crd_nc[(adpairs[j][k*4+m])][l];
	
	  theta=CDon(atom,pi);
	  if (theta > 2.0*pi)  theta = theta -2.0*pi;
	  else if (theta < 0.0) theta = theta +2.0*pi;

	  dtheta=theta-refd[nr][n];
	  if (dtheta < 0.0) dtheta=-dtheta;

	  a=2.0*pi-dtheta;

	  if (dtheta>a) dtheta=a;

	  deltad[nr]+=dtheta;
	  ++n;
	}
      }
      fprintf(outputfile,"%8.4lf  ",deltad[nr]/pi/numdihedall);
    }
    fprintf(outputfile,"\n");
  }
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
  printf("%s inputfilename refcrdfilename parmfilename outputfilename \n",progname);
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
