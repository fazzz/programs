
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>
#include <netcdf.h>
#include <getopt.h>

#include "netcdf_mineL.h"

#include "NC.h"
#include "NCcount.h"

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"
#include "GOLMAA_MB_PROTEINS2008.h"

//#include "GOLM_Clementi_set_wmCutOff.h"
//#include "GOLM_Clementi.h"
//#include "GOLM_Clementi_MB.h"
//#include "GOLM_Clementi_MB_wmCutOff.h"

#include "PTL.h"
#include "EF.h"


#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,a,dummy;
  int numstep;
  double Kai,aveKai=0.0,varKai=0.0;
  double de=1.0,d=1.0,d2,ep=0.3/*ep_natatt_hybrid*/;

  double *crd,*crdref1,*crdref2/*,*refcrdAA1,*refcrdAA2*/,vec[3];
  int numatom,numheavyatom,/*numCAatom,*/numres;

  int AMBERMODEflag=OFF,MODE=1,AAMODE=OFF,crdMODE=OFF,wobaimpflag=OFF;
  int nibnum=3;

  double criteria=criteria_NC;

  int DOF=1;
  double Tobj=300,KEobj,KBT;
  double k_B=1.98723e-3;
  double UNITT=418.4070;

  double length;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  struct potential e;
  struct force f;
  struct potential_GOLMAA_MB_PROTEINS2008 e_GOLM;
  //  struct potential_GOLM_Clementi_MB e_GOLM;
  double p_t=0.0,E_t;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id;

  char *line;
  size_t len=0;

  char *progname;
  char *inputfilename,*parmtopname,*crdreffilename1,*crdreffilename2;
  char *outputfilename;

  FILE *crdfile,*parmtop,*crdreffile1,*crdreffile2;
  FILE *outputfile;

  int opt_idx=1;

  struct option long_opt[] = {
    {"Amber",0,NULL,'A'},
    {"wobaimp",0,NULL,'w'},
    {"crd",0,NULL,'d'},
    {"DOF",1,NULL,'D'},
    {"temp",1,NULL,'t'},
    {"ep",1,NULL,'e'},
    {"de",1,NULL,'3'},
    {"dh",1,NULL,'4'},
    {"cutoff",1,NULL,'c'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hdAx12w:D:t:e:c:b:3:4:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'A':
      AMBERMODEflag=ON;
      break;
    case 'd':
      crdMODE=ON;
      break;
    case 'w':
      wobaimpflag=ON;
      break;
    case 'x':
      AAMODE=ON;
      break;
    case 't':
      Tobj=atof(optarg);
      break;
    case 'D':
      DOF=atoi(optarg);
      break;
    case 'b':
      nibnum=atoi(optarg);
      break;
    case 'c':
      criteria=atof(optarg);
      break;
    case 'e':
      ep=atof(optarg);
      break;
    case '3':
      de=atof(optarg);
      break;
    case '4':
      d=atof(optarg);
      break;
    default:
      USAGE(progname);
      exit(1);
    case 'h':
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
  inputfilename   = *argv;
  parmtopname     = *++argv;
  crdreffilename1 = *++argv;
  crdreffilename2 = *++argv;
  outputfilename  = *++argv;

  parmtop=efopen(parmtopname,"r");
  readParmtopL(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;
  j=0;
  for (i=0;i<numatom;++i) {
    //    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
    if (strncmp(AP.IGRAPH[i],"H",1)==0) {
      ++j;
    }
  }
  numheavyatom=numatom-j;
  //  numCAatom=j;
  numres=AP.NRES;

  crd=(double *)gcemalloc(sizeof(double)*/*numCAatom*/numatom*3);
  crdref1=(double *)gcemalloc(sizeof(double)*/*numCAatom*/numatom*3);
  crdref2=(double *)gcemalloc(sizeof(double)*/*numCAatom*/numatom*3);
  //  refcrdAA1=(double *)gcemalloc(sizeof(double)*numatom*3);
  //  refcrdAA2=(double *)gcemalloc(sizeof(double)*numatom*3);

  crdreffile1=efopen(crdreffilename1,"r");
  getline(&line,&len,crdreffile1);
  fscanf(crdreffile1,"%d",&dummy);
  //j=0;
  for (i=0;i<numatom;++i) {
    //    for (k=0;k<3;++k) fscanf(crdreffile1,"%lf",&refcrdAA1[i*3+k]);
    /********************************************************/
    /* if (strncmp(AP.IGRAPH[i],"CA",2)==0) {		    */
    /*   for (k=0;k<3;++k) crdref1[j*3+k]=refcrdAA1[i*3+k]; */
    /*   ++j;						    */
    /* }						    */
    /********************************************************/
    for (k=0;k<3;++k) fscanf(crdreffile1,"%lf",&crdref1[i*3+k]);
  }
  fclose(crdreffile1);

  crdreffile2=efopen(crdreffilename2,"r");
  getline(&line,&len,crdreffile2);
  fscanf(crdreffile2,"%d",&dummy);
  //j=0;
  for (i=0;i<numatom;++i) {
    //    for (k=0;k<3;++k) fscanf(crdreffile2,"%lf",&refcrdAA2[i*3+k]);
    /********************************************************/
    /* if (strncmp(AP.IGRAPH[i],"CA",2)==0) {		    */
    /*   for (k=0;k<3;++k) crdref2[j*3+k]=refcrdAA2[i*3+k]; */
    /*   ++j;						    */
    /* }						    */
    /********************************************************/
    for (k=0;k<3;++k) fscanf(crdreffile2,"%lf",&crdref2[i*3+k]);
  }
  fclose(crdreffile2);

  if (AMBERMODEflag==ON) numstep=myncL_get_present_step_AMBER(inputfilename,&nc_id);
  else numstep=myncL_get_present_step_MCD(inputfilename,&nc_id_MCD);

  ffL_set_calcffandforce(&e,&f);
  GOLMAA_MB_PROTEINS2008_ff_calcff_set(&e_GOLM,crdref1,crdref2,numatom,numres,e.parm.indexnb,e.parm.numnb,ep,nibnum,criteria);

  //  GOLM_Clementi_MB_ff_set_calcff_wmCutOff(&e_GOLM,crdref1,crdref2,refcrdAA1,refcrdAA2,/*numCAatom*/numatom,numatom,ep,criteria);

  KEobj=0.5*DOF*k_B*Tobj;
  KBT=k_B*Tobj;

  de=de*KBT;
  d=d*KBT;

  d2=d*d;

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    if (AMBERMODEflag==ON) myncL_open_inq_get_sh_AMBER(inputfilename,/*numCAatom*/numatom,i,1,i+1,&nc_id,crd_nc);
    else myncL_open_inq_get_sh_MCD(inputfilename,/*numCAatom*/numatom,i,1,i+1,&nc_id_MCD,crd_nc);

    for (j=0;j</*numCAatom*/numatom;++j)for (k=0;k<3;++k)crd[j*3+k]=crd_nc[j][k];

    //    Kai=GOLM_Clementi_MB_Kai(crd,/*numCAatom*/numatom,de,d,d2,&e_GOLM);

    Kai=GOLMAA_MB_PROTEINS2008_Kai(crd,numatom,de,d,d2,&e_GOLM);
    
    aveKai=(i*aveKai+Kai)/(i+1);
    varKai=(i*varKai+Kai*Kai)/(i+1);
    fprintf(outputfile,"%e\n",Kai);
  }
  fclose(outputfile);

  varKai=sqrt(varKai-aveKai*aveKai);

  printf("%s \n",inputfilename);
  printf("%s \n",outputfilename);
  printf("ave= %e var= %e\n",aveKai,varKai);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("[--Amber] \n");
  printf("[--nibnum] \n");
  printf("[--criteria] \n");
  printf("%s [-h] inputfilename crdreffilename1 crdreffilename2 parmfilename outputfilename \n",progname);
}
