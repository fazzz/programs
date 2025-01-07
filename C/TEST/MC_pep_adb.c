
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"
#include "IO.h"

#include "RAND.h"
#include "BOXMULL.h"

#include "MB.h"
#include "PT.h"
#include "FF.h"

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
double calc_direst(double *crd,int numatom,struct drestparameters drestparm, double *cv, double *cv_x);
double calc_ene_cv(double *cv,double *cv_ave,struct drestparameters drestparm);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int d;
  double f;
  int drestflag=OFF;
  int numarg;
  int numstep=10000,interval=10,intervalt=1,intervale=1,intervald=1;
  int numatom,numpara,numnb,num14;

  double pi;

  double beta=300.0,beta2=300.0,delta;
  double dx=0.01,dcv=0.01;
  double *crd,*crd_trial;
  double *cv,*cv_trial,*cv_x,*cv_x_trial,*cv_ave;
  double drest,drest_trial;
  double drest_cv,drest_trial_cv;
  struct potential ene,ene_trial;
  struct drestparameters drestparm;

  char *progname;
  char *inifilename,*parmtopname,*indexfilename,*trjfilename,*enefilename,*enedfilename,*drestfilename,*dtrjfilename,*cvtrjfilename;
  FILE *inifile,*parmtop,*indexfile,*trjfile,*enefile,*enedfile,*drestfile,*dtrjfile,*cvtrjfile;

  char *line,*dummy;
  size_t len=0;

  int c,c2;
  extern char *optarg;
  extern int optind,opterr,optopt;

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  pi=acos(-1.0);

  progname=argv[0];
  while((c=getopt(argc,argv,"hn:i:j:k:b:c:d:"))!=-1) {
    switch(c) {
    case 'h':
      USAGE(progname);
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
    case 'c':
      beta2=atof(optarg);
      break;
    case 'd':
      dx=atof(optarg);
      break;
    default:
      USAGE(progname);
      exit(1);
    }
  }

  argc-=optind;
  argv+=optind;

  numarg=11;

  if (argc < numarg) {
    USAGE(progname);
    exit(1);
  }
  numnb = atoi(*argv);
  num14 = atoi(*++argv);
  inifilename  = *++argv;
  parmtopname  = *++argv;
  indexfilename = *++argv;
  drestfilename = *++argv;
  enefilename  = *++argv;
  enedfilename = *++argv;
  trjfilename  = *++argv;
  dtrjfilename = *++argv;
  cvtrjfilename = *++argv;

  beta=1.0/(kb*beta);
  beta2=1.0/(kb*beta2);

  parmtop=efopen(parmtopname,"r");
  readParmtop(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;
  numpara=AP.NTYPES*(AP.NTYPES+1)/2;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  crd_trial=(double *)gcemalloc(sizeof(double)*numatom*3);
  inifile=efopen(inifilename,"r");
  io_scanconf_Amber_ini(inifile,numatom,crd);
  fclose(inifile);

  drestfile=efopen(drestfilename,"r");
  fscanf(drestfile,"%lf",&drestparm.k);
  fscanf(drestfile,"%d",&drestparm.numdrest);
  cv=(double *)gcemalloc(sizeof(double)*drestparm.numdrest);
  cv_trial=(double *)gcemalloc(sizeof(double)*drestparm.numdrest);
  cv_x=(double *)gcemalloc(sizeof(double)*drestparm.numdrest);
  cv_x_trial=(double *)gcemalloc(sizeof(double)*drestparm.numdrest);
  cv_ave=(double *)gcemalloc(sizeof(double)*drestparm.numdrest);
  drestparm.dihed_index=(int *)gcemalloc(sizeof(int)*drestparm.numdrest*4);
  for (i=0;i<drestparm.numdrest;++i) {
    for (j=0;j<4;++j) {
      fscanf(drestfile,"%d",&d);
      drestparm.dihed_index[i*4+j]=d-1;
    }
    fscanf(drestfile,"%lf",&f);
    cv[i]=wraped_angle(f/180*pi,pi);
  }
  fclose(drestfile);

  trjfile=efopen(trjfilename,"w");
  enefile=efopen(enefilename,"w");
  enedfile=efopen(enedfilename,"w");
  dtrjfile=efopen(dtrjfilename,"w");
  cvtrjfile=efopen(cvtrjfilename,"w");

  indexfile=efopen(indexfilename,"r");
  ff_set_calcff(numnb,num14,indexfile,&ene);
  fclose(indexfile);
  indexfile=efopen(indexfilename,"r");
  ff_set_calcff(numnb,num14,indexfile,&ene_trial);
  fclose(indexfile);
  ff_calcff(crd,numatom,&ene);
  drest=calc_direst(crd,numatom,drestparm,cv,cv_x);
  for (i=0;i<numatom*3;++i)
    crd_trial[i]=crd[i];
  drest_cv=drest;
  for (i=0;i<drestparm.numdrest;++i)
    cv_ave[i]=0.0;
  l=0;

  for (i=0;i<numstep;++i) {
    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k)
	crd_trial[j*3+k]=crd[j*3+k]+dx*Box_Muller(i,0.0,1.0);
      ff_calcff(crd_trial,numatom,&ene_trial);
      drest_trial=calc_direst(crd_trial,numatom,drestparm,cv,cv_x_trial);
      delta=ene_trial.p_t-ene.p_t+drest_trial-drest;
      if((c=Metropolis(beta*delta))==1) {
	for (k=0;k<numatom*3;++k)
	  crd[k]=crd_trial[k];
	memcpy(&ene,&ene_trial,sizeof(ene_trial));
	drest=drest_trial;
	for (k=0;k<drestparm.numdrest;++k)
	  cv_x[k]=cv_x_trial[k];
      }
      for (k=0;k<drestparm.numdrest;++k)
	cv_ave[k]=(l*cv_ave[k]+cv_x[k])/(l+1);
      ++l;
    }

    if (i%interval==0) {
      for (j=0;j<drestparm.numdrest;++j) {
	cv_trial[j]=cv[j]+dcv*Box_Muller(i,0.0,1.0);
	drest_trial_cv=calc_ene_cv(cv_trial,cv_ave,drestparm);
	drest_trial=calc_direst(crd,numatom,drestparm,cv,cv_x_trial);
	delta=drest_trial_cv-drest_cv;
	if((c2=Metropolis(beta2*delta))==1) {
	  for (k=0;k<drestparm.numdrest;++k)
	    cv[k]=crd_trial[k];
	  drest_cv=drest_trial_cv;
	  drest=drest_trial;
	}
	for (k=0;k<drestparm.numdrest;++k)
	  cv_ave[k]=0.0;
	l=0;
      }
    }

    if (i%intervalt==0) {
      io_outputconf(trjfile,numatom,crd,'x');
      fprintf(dtrjfile,"%d ",i+1);
      for (j=0;j<drestparm.numdrest;++j)
	fprintf(dtrjfile,"%12.8lf ",cv_x[j]);
      fprintf(dtrjfile,"\n");
      fprintf(cvtrjfile,"%d ",i+1);
      for (j=0;j<drestparm.numdrest;++j)
	fprintf(cvtrjfile,"%12.8lf ",cv[j]);
      fprintf(cvtrjfile,"\n");
    }
    if (i%intervald==0) {
      fprintf(enedfile,"%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",ene.p_t,ene.p_e_t,ene.p_LJ_t,ene.p_e_14_t,ene.p_LJ_14_t,ene.p_d_t,ene.p_a_t,ene.p_b_t,drest);
    }
    if (i%intervale==0) {
      fprintf(enefile,"%d %12.8lf %12.8lf %12.8lf %d\n",i+1,ene.p_t+drest,ene.p_t,drest,c);
    }
  }
    
  fclose(trjfile);
  fclose(dtrjfile);
  fclose(cvtrjfile);
  fclose(enedfile);
  fclose(enefile);
  
  return 0;
}
 
void USAGE(char *progname) {
  printf("USAGE: %s numnb num14 inifilename parmtopname indexfilename enefilename enedfilename logfilename \n", progname);
}

double calc_direst(double *crd,int numatom,struct drestparameters drestparm, double *cv, double *cv_x){
  int i,j,k;
  int numdrest;
  double atom[4][3],dihedang,dang;
  double drest=0.0,pi;

  pi=acos(-1.0);
  numdrest=drestparm.numdrest;

  for (i=0;i<numdrest;++i) {
    for (j=0;j<4;++j)
      for (k=0;k<3;++k)
      atom[j][k]=crd[drestparm.dihed_index[i*4+j]*3+k];
    cv_x[i] = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);

    if ((dang=cv_x[i]-/*drestparm.dihed_equ*/cv[i])>pi)
      dang-=pi;
    else if (dang<-1.0*pi)
      dang+=pi;

    drest+=0.5*(drestparm.k)*dang*dang;
  }

  return drest;
}

double calc_ene_cv(double *cv,double *cv_ave,struct drestparameters drestparm){
  int i,j,k;
  int numdrest;
  double dang;
  double drest=0.0,pi;

  pi=acos(-1.0);
  numdrest=drestparm.numdrest;

  for (i=0;i<numdrest;++i) {
    if ((dang=cv[i]-cv_ave[i])>pi)
      dang-=pi;
    else if (dang<-1.0*pi)
      dang+=pi;

    drest+=0.5*(drestparm.k)*dang*dang;
  }
  return drest;
}
