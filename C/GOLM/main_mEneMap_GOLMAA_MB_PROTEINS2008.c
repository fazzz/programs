#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "ABAb.h"
#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"
#include "GOLMAA_MB_PROTEINS2008.h"

#include "PTL.h"
#include "FFL.h"
#include "EF.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,dummyd;
  int numatom,numheavyatom,numres,numclut;
  double dihed,da=0.1;
  double dummy;

  CLTb *clt;

  double Tobj=300,KEobj,KBT,BE;
  double k_B=1.98723e-3;
  double UNITT=418.4070;
  int DOF;

  double *qrot;

  double ep=0.3;
  int NCmode=3,nibnum=3,criteria=6.5;

  double *crd,*refcrd1,*refcrd2,*mass;
  double d=1.0,de=1.0,d2;

  double Kai;
  double *Kai_s,*p_t_s;

  double min,max,frame,width=0.05,sum;
  double *hist;

  struct potential e;
  struct force f;
  struct potential_GOLMAA_MB_PROTEINS2008 e_GOLM;
  double p_t=0.0;

  char *line;
  size_t len=0;

  double pi;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*reffilename1,*reffilename2,*parmfilename,*clustfilename;
  char *outputfilename,*outputfilename2;

  FILE *inputfile,*reffile1,*reffile2,*velfile,*parmfile,*clustfile;
  FILE *outputfile,*outputfile2;

  FILE *logfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"ep",1,NULL,'e'},
    {"da",1,NULL,'@'},
    {"width",1,NULL,'w'},
    {"temp",1,NULL,'t'},
    {"cutoff",1,NULL,'c'},
    {"de",1,NULL,'d'},
    {"dh",1,NULL,'2'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"he:w:@:t:c:d:2:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case '@':
      da=atof(optarg);
      break;
    case 'w':
      width=atof(optarg);
      break;
    case 'e':
      ep=atof(optarg);
      break;
    case 't':
      Tobj=atof(optarg);
      break;
    case 'x':
      da=atof(optarg);
      break;
    case 'c':
      criteria=atof(optarg);
      break;
    case 'd':
      de=atof(optarg);
      break;
    case '2':
      d=atof(optarg);
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  pi=acos(-1.0);

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 7) {
    USAGE(progname);
    exit(1);
  }
  inputfilename   = *argv;
  reffilename1    = *++argv;
  reffilename2    = *++argv;
  clustfilename   = *++argv;
  parmfilename    = *++argv;
  outputfilename  = *++argv;
  outputfilename2 = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;

  j=0;
  for (i=0;i<numatom;++i) {
    if (strncmp(AP.IGRAPH[i],"H",1)==0) {
      ++j;
    }
  }
  numheavyatom=numatom-j;
  numres=AP.NRES;
  mass=(double *)gcemalloc(sizeof(double)*numatom);
  for (i=0;i<numatom;++i) mass[i]=AP.AMASS[i];
  
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrd1=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrd2=(double *)gcemalloc(sizeof(double)*numatom*3);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&dummyd);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);

  reffile1=efopen(reffilename1,"r");
  getline(&line,&len,reffile1);
  fscanf(reffile1,"%d",&dummyd);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(reffile1,"%lf",&refcrd1[i*3+j]);
  fclose(reffile1);

  reffile2=efopen(reffilename2,"r");
  getline(&line,&len,reffile2);
  fscanf(reffile2,"%d",&dummyd);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(reffile2,"%lf",&refcrd2[i*3+j]);
  fclose(reffile2);

  clustfile=efopen(clustfilename,"r");
  clt=ABAbp_clustscan(clustfile,&numclut);
  fclose(clustfile);

  ffL_set_calcffandforce(&e,&f);
  GOLMAA_MB_PROTEINS2008_ff_calcff_set(&e_GOLM,refcrd1,refcrd2,numatom,numres,e.parm.indexnb,e.parm.numnb,ep,nibnum,criteria);

  qrot=(double *)gcemalloc(sizeof(double)*numclut);

  DOF=(numclut-1); DOF+=6; KEobj=0.5*DOF*k_B*Tobj; KBT=k_B*Tobj;

  de=de*KBT;  d=d*KBT;  d2=d*d;
  GOLMAA_MB_PROTEINS2008_ff_calcff_wobaimp(crd,numatom,de,d2,&e_GOLM);

  Kai_s=(double *)gcemalloc(sizeof(double)*1);
  p_t_s=(double *)gcemalloc(sizeof(double)*1);

  outputfile=efopen(outputfilename,"w");
  l=0;
  for (i=1;i<numclut;++i) {
    for (j=0;j<numclut;++j)  qrot[j]=0.0;
    for (dihed=0.0;dihed<2.0*pi;dihed+=da) {

      qrot[i]=dihed;

      ABAb_update(clt,crd,qrot,numclut,numatom);
    
      GOLMAA_MB_PROTEINS2008_ff_calcff_wobaimp(crd,numatom,de,d2,&e_GOLM);
      Kai=GOLMAA_MB_PROTEINS2008_Kai(crd,numatom,de,d,d2,&e_GOLM);
    
      p_t=e_GOLM.p_MB;

      Kai_s=(double *)gcerealloc(Kai_s,sizeof(double)*(l+1));
      p_t_s=(double *)gcerealloc(p_t_s,sizeof(double)*(l+1));

      Kai_s[l]=Kai;
      p_t_s[l]=p_t;

      if (l==0) {
	min=Kai;
	max=Kai;
      }
      else {
	if (min>Kai && Kai > -50) min=Kai;
	if (max<Kai && Kai <  50) max=Kai;
      }

      ++l;
      fprintf(outputfile,"%e %e \n",p_t,Kai);
    }
    fprintf(outputfile,"\n");
  }
  fclose(outputfile);

  frame=(int)round((max-min)/width)+1;
  hist=(double *)gcemalloc(sizeof(double)*(frame));

  for (i=0;i<l;++i) {
    if ( Kai_s[i] >= -50 && Kai_s[i] <= 50 ) {
      BE=exp(-1.0/KBT*p_t_s[i]);
      hist[((int)round((Kai_s[i]-min)/width))]+=1.0*BE;
    }
  }

  sum=0.0;  for (i=0;i<frame;++i)  sum+=hist[i];

  //  for (i=0;i<frame;++i) hist[i]/=sum;

  outputfile2=efopen(outputfilename2,"w");
  for (i=0;i<frame;++i) {
    fprintf(outputfile2,"%e %e \n",min+width*i,hist[i]);
  }
  fclose(outputfile2);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename refcrdfilename clustfilename parmfilename outputfilename1 \n",progname);
}
