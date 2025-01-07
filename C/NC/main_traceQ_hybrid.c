
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "NC.h"
#include "PTL.h"
#include "EF.h"
#include "TOPO.h"
#include "MB.h"
#include "GOLMAA_hybrid_set.h"

#include "netcdf_mine.h"

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int na,resi,resj;
  int numstep,HMODE=EXC;
  double Q;

  double *crd,*crdref,vec[3];
  int numatom,numres;

  int AMBERMODEflag=OFF;
  int nibnum=1;

  double criteria=criteria_NC;
  int numncaa,numncres;
  int  *index_natatt,**ncmap,**ncmapres;
  double atom1[3],atom2[3];

  double *cradii_natatt;
  double length;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id;

  char *progname;
  char *inputfilename,*parmtopname,*crdreffilename;
  char *outputfilename;

  FILE *crdfile,*parmtop,*crdreffile;
  FILE *outputfile;

  while((c=getopt(argc,argv,"hAc:b:k:b:c:"))!=-1) {
    switch(c) {
    case 'A':
      AMBERMODEflag=ON;
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    case 'b':
      nibnum=atoi(optarg);
      break;
    case 'c':
      criteria=atof(optarg);
      break;
    default:
      USAGE(progname);
      exit(1);
    }
  }

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  inputfilename = *argv;
  parmtopname = *++argv;
  crdreffilename = *++argv;
  outputfilename = *++argv;

  parmtop=efopen(parmtopname,"r");
  readParmtopL(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;
  numres=AP.NRES;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  crdref=(double *)gcemalloc(sizeof(double)*numatom*3);

  crdreffile=efopen(crdreffilename,"r");
  io_scanconf_Amber_ini(crdreffile,numatom,crdref);
  fclose(crdreffile);

  if (AMBERMODEflag==ON)
    numstep=mync_get_present_step_AMBER(inputfilename,&nc_id);
  else
    numstep=mync_get_present_step_MCD(inputfilename,&nc_id_MCD);

  //  ncmap=GOLMAA_hybrid_ff_set_make_native_contact(crdref,criteria,&numncaa,numatom,numres);
  ncmap=GOLMAA_hybrid_ff_set_make_native_contact_6(crdref,criteria,&numncaa,numatom,numres,nibnum);
  cradii_natatt=(double *)gcemalloc(sizeof(double)*numncaa);

  ncmapres=(int **)gcemalloc(sizeof(int *)*numres);
  for (i=0;i<numres;++i) ncmapres[i]=(int *)gcemalloc(sizeof(int)*numres);

  for (i=0;i<numres;++i) for (j=0;j<numres;++j) ncmapres[i][j]=-1;
      
  na=0;
  numncres=0;
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (ncmap[i][j]==0 ) {
	length = 0.0;
	for(k=0;k<3;++k){
	  vec[k] = crdref[i*3+k]-crdref[j*3+k];
	  length += vec[k]*vec[k];
	}
	length = sqrt(length);
	cradii_natatt[na]=length;
	++na;
	resi=PTL_resnum(i)-1;
	resj=PTL_resnum(j)-1;
	if (ncmapres[resi][resj]==-1) ncmapres[resi][resj]=1;
	else ++numncres;
      }
    }
  }

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    //      mync_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id,crd_nc);
    //    mync_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id,crd_nc);
    if (AMBERMODEflag==ON)
      mync_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id,crd_nc);
    else
      mync_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd_nc);
    for (j=0;j<numatom;++j)for (k=0;k<3;++k)crd[j*3+k]=crd_nc[j][k];
    Q=count_native_contact_hybrid(numncaa,crd,numatom,numres,ncmap,cradii_natatt);
    
    fprintf(outputfile,"%e\n",Q);
  }
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename parmfilename outputfilename \n",progname);
}
