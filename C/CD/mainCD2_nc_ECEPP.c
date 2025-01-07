#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "netcdf_mine.h"

#include "IO.h"
#include "ECEPE.h"
#include "PDB.h"
#include "TOPO.h"
#include "PT.h"
#include "EF.h"

#define ON 1
#define OFF 0

int usage(void);

int main(int argc, char *argv[]) {
  int i,j,k,l,d,m,nt,dummy,num,n;
  int flagcn='n',flagrst=OFF;
  int numatom,numstep;
  int interval=1;

  double atom[4][3];
  double *protcoord;
  double theta,*old_theta;
  double pi;

  double *co;
  double *dihed,*delta_dihed;
  double dihed_check;
  double *co_dummy,*dihed_dummy,*ene;

  int *bp;
  int **bp_f;
  int *numb;
  
  struct ECEPE_parms ECEPE_p;
  struct pnb nb_p;
  struct ECEPE_pote p;
  struct ECEPE_force f;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc[MAXATOM][3],*crd;

  char *line;
  size_t len=0;

  char *inputfilename,*outputfilename;
  char *preofilename,*bd8filename,*coofilename_for_sflag;
  char *progname;
  FILE *inputfile,*outputfile;
  FILE *preofile,*bd8file,*coofile_for_sflag;

  while((c=getopt(argc,argv,"hcnsv:"))!=-1) {
    switch(c) {
    case 'v':
      interval=atoi(optarg);
      break;
    case 's':
      flagrst=ON;
      break;
    case 'c':
      flagcn='c';
      break;
    case 'n':
      flagcn='n';
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

  if (argc < 5) {
    USAGE(progname);
    exit(1);
  }
  inputfilename  = *argv;
  preofilename   = *++argv;
  bd8filename    = *++argv;
  coofilename_for_sflag = *++argv;
  outputfilename = *++argv;

  read_ECEPE_parm(preofilename,bd8filename,&ECEPE_p,&nb_p);
  numatom=ECEPE_p.NUMATM;

  pi=acos(-1.0);

  bp=(int *)gcemalloc(sizeof(int)*(ECEPE_p.NUMATM-1)*2);
  bp_f=(int **)gcemalloc(sizeof(int *)*ECEPE_p.NUMATM);

  numb=(int *)gcemalloc(sizeof(int)*ECEPE_p.NUMATM);
  make_bd_pair_list(ECEPE_p,bp,bp_f,numb);
  make_dihed_pairs_list(ECEPE_p,bp_f,numb);

  co=(double *)gcemalloc(sizeof(double)*(ECEPE_p.NUMATM)*3);
  dihed=(double *)gcemalloc(sizeof(double)*(ECEPE_p.NUMVAR));
  delta_dihed=(double *)gcemalloc(sizeof(double)*(ECEPE_p.NUMVAR));
  for (i=0;i<ECEPE_p.NUMVAR;++i) delta_dihed[i]=0.0;

  co_dummy=(double *)gcemalloc(sizeof(double)*40*3);
  dihed_dummy=(double *)gcemalloc(sizeof(double)*10);
  ene=(double *)gcemalloc(sizeof(double)*6);

  coofile_for_sflag=efopen(coofilename_for_sflag,"r");
  read_ECEPE_detail_coo_cyc(coofile_for_sflag,co_dummy,dihed_dummy,ene);
  for (i=0;i<ECEPE_p.NUMVAR;++i) {
    dihed_dummy[i]=dihed_dummy[i]*pi/180.0;
    if (dihed_dummy[i]<-pi)
      dihed_dummy[i]+=2.0*pi;
    else if (dihed_dummy[i]>pi)
      dihed_dummy[i]-=2.0*pi;
  }
  k=0;
  for (i=0;i<ECEPE_p.NUMATM;++i) {
    for (j=0;j<3;++j) {
      co[i*3+j]=co_dummy[k];
      ++k;
    }
  }    
  calc_TORS_for_get_sabun(ECEPE_p.NUMVAR,co,ECEPE_p,dihed_dummy,delta_dihed);
  fclose(coofile_for_sflag);
  
  old_theta=(int *)gcemalloc(sizeof(int)*ECEPE_p.NUMVAR);

  inputfile=efopen(inputfilename,"r");
  if (flagrst==ON) fscanf(inputfile,"%d",&n);
  outputfile=efopen(outputfilename,"w");
  i=0;
  numstep=0;
  crd=(double *)gcemalloc(sizeof(double )*numatom*3);
  d=0;
  while (d!=-1) {
    d=io_scanconfwj(inputfile,numatom,crd,'x');
    /***************************************/
    /* for (j=0;j<numatom;++j)		   */
    /*   for (k=0;k<3;++k)		   */
    /* 	crd_nc[j][k]=crd[j*3+k];	   */
    /***************************************/
    for (j=0;j<ECEPE_p.NUMATM;++j)
      for (k=0;k<3;++k) 
	crd_nc[j][k]=crd[(ECEPE_p.atom[j].katom-1)*3+k];
    ++numstep;

    if ( (numstep%interval) == 0 ) {
      for (k=0;k<ECEPE_p.NUMVAR;++k) {
	for (m=0;m<4;++m) 
	  for (l=0;l<3;++l)
	    atom[m][l]=crd_nc[(ECEPE_p.dihed[k].dpairs[m])][l];

	if (flagcn=='n' || i==0 ) {
	  theta=dih(atom[0],atom[1],atom[2],atom[3]);
	  theta+=delta_dihed[k];
	  if (theta > pi)  theta-=2.0*pi;
	  else if (theta < -pi)  theta+=2.0*pi;
	}
	else if (flagcn=='c') {
	  theta=dih(atom[0],atom[1],atom[2],atom[3]);
	  theta+=delta_dihed[k];
	  if (theta > pi)  theta-=2.0*pi;
	  else if (theta < -pi)  theta+=2.0*pi;
	  if (fabs(theta-old_theta[k]) < fabs((theta-2.0*pi)-old_theta[k])) {
	    if (fabs(theta-old_theta[k]) < fabs((theta+2.0*pi)-old_theta[k])) {
	      ;
	    }
	    else {
	      theta=theta+2.0*pi;
	    }
	  }
	  else {
	    if (fabs((theta-2.0*pi)-old_theta[k]) < fabs((theta+2.0*pi)-old_theta[k])) {
	      theta=theta-2.0*pi;
	    }
	    else {
	      theta=theta+2.0*pi;
	    }
	  }
	}
	
	fprintf(outputfile,"%e ",theta*180/pi);
	old_theta[k]=theta;
      }
      fprintf(outputfile,"\n");
      ++i;
    }
  }
  fclose(inputfile);
  fclose(outputfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h ] help \n");
  printf("[-s ] flag for rst \n");
  printf("[-k interval] specify interval \n");
  printf("%s inputfilename preofilename bd8filename coofilename_for_sflag outputfilename \n",progname);
}


 
