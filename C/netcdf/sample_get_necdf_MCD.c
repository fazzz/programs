
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "EF.h"
#include "IO.h"

#include "TOPO.h"
#include "RAND.h"
#include "BOXMULL.h"

#include "PT.h"

#include "netcdf_mine.h"

#define ON 1
#define OFF 0

#define kb 1.98723e-3

struct drestparameters {
  int numdrest;
  double k;
  double *dihed_equ;

  int *dihed_index;
};

void USAGE(char *progname);
void USAGEopt(char *progname);
double calc_direst(double *crd,int numatom,struct drestparameters drestparm, double *cv, double *cv_x);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int d;
  double f;
  int drestflag=OFF;
  int numarg;
  int numstep=10000,intervalt=1,intervale=1,intervald=1;
  int numatom,numpara,numnb,num14;

  double pi;

  double beta=300.0,delta;
  double dx=0.1;
  double *crd,*crd_trial;
  double *cv,*cv_trial,*cv_x,*cv_x_trial;
  double drest,drest_trial;
  struct potential ene,ene_trial;
  struct drestparameters drestparm;

  char *progname;
  char *inifilename,*parmtopname,*outfilename,*drestfilename;
  FILE *inifile,*parmtop,*drestfile;

  double **crd_nc,*energy;
  struct my_netcdf_out_id_MCD nc_id_MCD;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  pi=acos(-1.0);

  progname=argv[0];
  while((c=getopt(argc,argv,"chan:i:j:k:b:d:"))!=-1) {
    switch(c) {
    case 'h':
      USAGE(progname);
      exit(1);
    case 'c':
      USAGEopt(progname);
      exit(1);
    case 'n':
      numstep=atoi(optarg);
      break;
    case 'i':
      intervalt=atoi(optarg);
      break;
    case 'j':
      intervale=atoi(optarg);
      break;
    case 'k':
      intervald=atoi(optarg);
      break;
    case 'b':
      beta=atof(optarg);
      break;
    case 'd':
      dx=atof(optarg);
      break;
    case 'a':
      drestflag=ON;
      break;
    default:
      USAGE(progname);
      exit(1);
    }
  }

  argc-=optind;
  argv+=optind;

  numarg=3;
  if (drestflag==ON)
    ++numarg;

  if (argc < numarg) {
    USAGE(progname);
    exit(1);
  }
  outfilename  = *argv;

  readParmtop(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;
  numpara=AP.NTYPES*(AP.NTYPES+1)/2;

  crd=(double *)gcemalloc(sizeof(double )*numatom*3);
  crd_trial=(double *)gcemalloc(sizeof(double )*numatom*3);

  crd_nc=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) crd_nc[i]=(double *)gcemalloc(sizeof(double)*3);
  energy=(double *)gcemalloc(sizeof(double)*(NDIMS_ENERGY-1));

  inifile=efopen(inifilename,"r");
  io_scanconf_Amber_ini(inifile,numatom,crd);
  fclose(inifile);

  if (drestflag==ON) {
    drestfile=efopen(drestfilename,"r");
    fscanf(drestfile,"%lf",&drestparm.k);
    fscanf(drestfile,"%d",&drestparm.numdrest);
    cv=(double *)gcemalloc(sizeof(double)*drestparm.numdrest);
    cv_trial=(double *)gcemalloc(sizeof(double)*drestparm.numdrest);
    cv_x=(double *)gcemalloc(sizeof(double)*drestparm.numdrest);
    cv_x_trial=(double *)gcemalloc(sizeof(double)*drestparm.numdrest);
    drestparm.dihed_index=(int *)gcemalloc(sizeof(int)*drestparm.numdrest*4);
    for (i=0;i<drestparm.numdrest;++i) {
      for (j=0;j<4;++j) {
	fscanf(drestfile,"%d",&d);
	drestparm.dihed_index[i*4+j]=d-1;
      }
      fscanf(drestfile,"%lf",&f);
      cv[i]=f/180*pi;
      if (cv[i]<-pi)
	cv[i]+=2.0*pi;
      else if (cv[i]>pi)
	cv[i]-=2.0*pi;
    }
    fclose(drestfile);
  }

  mync_create_def_MCD(outfilename,numatom,&nc_id_MCD);
  
  ff_set_calcffsp(&ene);
  ff_set_calcffsp(&ene_trial);

  ff_calcff(crd,numatom,&ene);
  if (drestflag == ON)
    drest=calc_direst(crd,numatom,drestparm,cv,cv_x);
  for (i=0;i<numatom*3;++i)
    crd_trial[i]=crd[i];
  
  l=0;
  for (i=0;i<numstep;++i) {
    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k)
	crd_trial[j*3+k]=crd[j*3+k]+dx*Box_Muller(i,0.0,1.0);
      ff_calcff(crd_trial,numatom,&ene_trial);
      delta=ene_trial.p_t-ene.p_t;
      if (drestflag == ON) {
	drest_trial=calc_direst(crd_trial,numatom,drestparm,cv,cv_x_trial);
	delta+=drest_trial-drest;
      }
      if((c=Metropolis(beta*delta))==1) {
	for (k=0;k<numatom*3;++k)
	  crd[k]=crd_trial[k];
	memcpy(&ene,&ene_trial,sizeof(ene_trial));
	if (drestflag == ON) {
	  drest=drest_trial;
	  for (k=0;k<drestparm.numdrest;++k)
	    cv_x[k]=cv_x_trial[k];
	}
      }
    }

    if (i%intervalt==0) {
      for (j=0;j<numatom;++j)
	for (k=0;k<3;++k) 
	  crd_nc[j][k]=crd[j*3+k];

      mync_put_crd_ene_MCD(nc_id_MCD,l,crd_nc,ene);
      ++l;
    }
  }
    
  nc_close((nc_id_MCD.ncid));

  return 0;
}
 
void USAGE(char *progname) {
  printf("[-h] -- help\n");
  printf("[-c] -- help for option\n");
  printf("[-a drestfilename ] -- drestflag\n");
  printf("[-n numstep ]\n");
  printf("[-i intervalt ]\n");
  printf("[-j intervale ]\n");
  printf("[-k intervald ]\n");
  printf("[-b beta ]\n");
  printf("[-d dx ]\n");
  printf("USAGE: %s inifilename parmtopname outfilename\n", progname);
}

void USAGEopt(char *progname) {
  printf("drestfile\n");
  printf("force constsnt for dihed restraint\n");
  printf("num of dihed restraint \n");
  printf("atom1 atom2 atom3 atom4 dihed \n");
}

double calc_direst(double *crd,int numatom,struct drestparameters drestparm, double *cv, double *cv_x){
  int i,j,k;
  int numdrest;
  double atom[4][3],dihedang,dang;
  double drest=0.0,pi,f;

  pi=acos(-1.0);
  numdrest=drestparm.numdrest;

  for (i=0;i<numdrest;++i) {
    for (j=0;j<4;++j)
      for (k=0;k<3;++k)
	atom[j][k]=crd[drestparm.dihed_index[i*4+j]*3+k];
    f = dih(atom[0],atom[1],atom[2],atom[3]);
    if (f>pi)
      f-=2.0*pi;
    else if (f<-1.0*pi)
      f+=2.0*pi;
    cv_x[i]=f;

    if ((dang=cv_x[i]-cv[i])>pi)
      dang-=2.0*pi;
    else if ((dang=cv_x[i]-cv[i])<-1.0*pi)
      dang+=2.0*pi;

    drest+=0.5*(drestparm.k)*dang*dang;
  }

  return drest;
}


