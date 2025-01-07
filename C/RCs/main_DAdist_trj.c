
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
  int numref=2,numatom,*numdihed,numdihedall;
  int numstep=1;
  int amberncflag=OFF,flag=1;
  int normflag=ON;
  double theta,dtheta,a;
  int **adpairs;

  double *deltad,pi;
  double **refcrd,**refd;

  double atom[4][3];
  double crd_nc[MAXATOM][3],*crd;
  //  struct my_netcdf_out_id_MCD nc_id_MCD;
  //  struct my_netcdf_out_id_AMBER nc_id_MD;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,**refcrdfilename,*parmfilename,*outputfilename;

  FILE *inputfile,**refcrdfile,*parmfile,*outputfile,*logfile;

  char *progname;
  int opt_idx=1;

  pi=acos(-1.0);

  struct option long_opt[] = {
    {"A",0,NULL,'A'},
    {"N",0,NULL,'N'},
    {"P",0,NULL,'P'},
    {"O",0,NULL,'O'},
    {"K",0,NULL,'K'},
    {"nref",1,NULL,'n'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hANPOKn:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'A':
      amberncflag=ON;
      break;
    case 'N':
      normflag=OFF;
      break;
    case 'P':
      flag=1;
      break;
    case 'O':
      flag=2;
      break;
    case 'K':
      flag=5;
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

  if (argc < 4+numref) {
    USAGE(progname);
    exit(1);
  }
  numstep          = atoi(*argv);
  inputfilename    = *++argv;
  for (i=0;i<numref;++i) refcrdfilename[i] = *++argv;
  parmfilename     = *++argv;
  outputfilename   = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;

  adpairs=(int **)gcemalloc(sizeof(int *)*5);
  adpairs[0]=(int *)gcemalloc(sizeof(int)*4); // PHI,PSI
  adpairs[1]=(int *)gcemalloc(sizeof(int)*4); // OMEGA
  adpairs[2]=(int *)gcemalloc(sizeof(int)*4); // KI
  adpairs[3]=(int *)gcemalloc(sizeof(int)*8); // ACE
  adpairs[4]=(int *)gcemalloc(sizeof(int)*8); // NME
  numdihed=(int *)gcemalloc(sizeof(int)*5);  
  readdihedpairsL(adpairs,numdihed);

  refcrd=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numref;++i) {
    refcrdfile[i]=efopen(refcrdfilename[i],"r");
    refcrd[i]=(double *)gcemalloc(sizeof(double)*numatom*3);
    getline(&line,&len,refcrdfile[i]);
    fscanf(refcrdfile[i],"%d",&d);
    for (j=0;j<numatom;++j) for (k=0;k<3;++k) fscanf(refcrdfile[i],"%lf",&refcrd[i][j*3+k]);
    fclose(refcrdfile[i]);
  }

  numdihedall=0;for (i=0;i<flag;++i) numdihedall+=numdihed[i];

  for (i=0;i<numref;++i) refd[i]=(double *)gcemalloc(sizeof(double)*numdihedall);
  
  for (n=0;n<numref;++n) {
    m=0;
    for (i=0;i<flag;++i) {
      for (j=0;j<numdihed[i];++j) {
	for (k=0;k<4;++k) for (l=0;l<3;++l) atom[k][l]=refcrd[n][(adpairs[i][j*4+k])*3+l];
    
	theta=CDon(atom,pi);
	if (theta > 2.0*pi)  theta = theta -2.0*pi;
	else if (theta < 0.0)  theta = theta +2.0*pi;

	refd[n][m]=theta;

	++m;
      }
    }
  }

  //  if (amberncflag==OFF) numstep=myncL_get_present_step_MCD(inputfilename,&nc_id_MCD);
  //  else numstep=myncL_get_present_step_AMBER(inputfilename,&nc_id_MD);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    /*****************************************************************************************************/
    /* if (amberncflag==OFF) myncL_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd_nc); */
    /* else myncL_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id_MD,crd_nc);		 */
    /*****************************************************************************************************/

    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k) {
	fscanf(inputfile,"%lf",&crd_nc[j][k]);
      }
    }

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
      if ( normflag==ON )
	fprintf(outputfile,"%8.4lf  ",deltad[nr]/pi/numdihedall);
      else
	fprintf(outputfile,"%8.4lf  ",deltad[nr]/pi);
    }
    fprintf(outputfile,"\n");
  }
  fclose(inputfile);
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
