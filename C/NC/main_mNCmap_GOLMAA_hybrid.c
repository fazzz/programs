
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "GOLMAA_hybrid_set.h"
#include "NC.h"
#include "PTL.h"
#include "EF.h"
#include "TOPO.h"

#include "netcdf_mine.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d,na,ii,jj;
  int resi,resj,resk;

  double *crdref;
  int numatom,numres;
  int dummy;

  int max;
  double sc;

  int flag=3,nibnum=1;

  double criteria=criteria_NC;
  int numncaa,numncres;
  int  **ncmap,**ncmapres;

  double vec[3];
  double length;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *progname;
  char *crdfilename,*parmtopname,*crdreffilename;
  char *outputfilename;

  FILE *crdfile,*parmtop,*crdreffile;
  FILE *outputfile;

  while((c=getopt(argc,argv,"hc:n:b:"))!=-1) {
    switch(c) {
    case 'b':
      nibnum=atoi(optarg);
      break;
    case 'n':
      flag=atoi(optarg);
      break;
    case 'c':
      criteria=atof(optarg);
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

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  crdreffilename = *argv;
  parmtopname = *++argv;
  outputfilename = *++argv;

  parmtop=efopen(parmtopname,"r");
  readParmtopL(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;
  numres=AP.NRES;

  crdref=(double *)gcemalloc(sizeof(double)*numatom*3);

  crdreffile=efopen(crdreffilename,"r");
  io_scanconf_Amber_ini(crdreffile,numatom,crdref);
  fclose(crdreffile);

  ncmapres=(int **)gcemalloc(sizeof(int *)*numres);
  for (i=0;i<numres;++i) ncmapres[i]=(int *)gcemalloc(sizeof(int)*numres);
  for (i=0;i<numres;++i) for (j=0;j<numres;++j) ncmapres[i][j]=/*-10*/0;

  if (flag==1)
    ncmap=GOLMAA_hybrid_ff_set_make_native_contact(crdref,criteria,&numncaa,numatom,numres);
  else if (flag==2)
    ncmap=GOLMAA_hybrid_ff_set_make_native_contact_2(crdref,criteria,&numncaa,numatom,numres);
  else if (flag==3)
    ncmap=GOLMAA_hybrid_ff_set_make_native_contact_3(crdref,criteria,&numncaa,numatom,numres);
  else if (flag==4)
    ncmap=GOLMAA_hybrid_ff_set_make_native_contact_4(crdref,criteria,&numncaa,numatom,numres);
  else if (flag==5)
    ncmap=GOLMAA_hybrid_ff_set_make_native_contact_5(crdref,criteria,&numncaa,numatom,numres,nibnum);
  else if (flag==6)
    ncmap=GOLMAA_hybrid_ff_set_make_native_contact_6(crdref,criteria,&numncaa,numatom,numres,nibnum);

  numncres=0;
  na=0;
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (ncmap[i][j]==0 ) {
	resi=PTL_resnum(i)-1;
	resj=PTL_resnum(j)-1;
	//	if (ncmapres[resi][resj]==-100) ncmapres[resi][resj]=1;
	/*else*/ ncmapres[resi][resj]+=1;
	if (ncmapres[resi][resj]==1) ++numncres;
      }
    }
  }

  max=0;
  for (i=0;i<numres;++i) 
    for (j=i+1;j<numres;++j) 
      if (max < ncmapres[i][j]) max=ncmapres[i][j];

  if (max!=0)
    sc=10000/max;
  else
    sc=0;

  for (i=0;i<numres;++i) 
    for (j=i+1;j<numres;++j) 
      ncmapres[i][j]=(int)(sc*ncmapres[i][j]);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numres;++i) {
    for (j=0;j<numres;++j) {
      fprintf(outputfile,"%2d ",ncmapres[i][j]);
    }
  }
  fclose(outputfile);

  printf("%d %d\n",numncaa,numncres);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] reffilename parmfilename \n",progname);
}

