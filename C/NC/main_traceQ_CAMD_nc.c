
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "GOLM_Clementi_set.h"
#include "PTL.h"
#include "TOPO.h"
#include "EF.h"

#include "netcdf_mine.h"

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d,na,ii,jj,nc;
  int numstep;
  double Q;
  double length;

  double *crd,*crdref,*crdrefCA;
  int numatom,numCAatom;

  double criteria=6.5;
  int numnc;
  int  *index_nc,**ncmap;
  double atom1[3],atom2[3];

  double *cradii_ca;
  double vec[3];

  char *line;
  size_t le=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id;

  char *progname;
  char *inputfilename,*parmtopname,*crdreffilename;
  char *outputfilename;

  FILE *crdfile,*parmtop,*crdreffile;
  FILE *outputfile;

  while((c=getopt(argc,argv,"hc:"))!=-1) {
    switch(c) {
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
  j=0;
  for (i=0;i<numatom;++i) if (strncmp(AP.IGRAPH[i],"CA",2)==0) ++j;
  numCAatom=j;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  crdref=(double *)gcemalloc(sizeof(double)*numatom*3);
  crdrefCA=(double *)gcemalloc(sizeof(double)*numCAatom*3);

  crdreffile=efopen(crdreffilename,"r");
  getline(&line,&le,crdreffile);
  fscanf(crdreffile,"%d",&d);
  j=0;
  for (i=0;i<numatom;++i) {
    for (k=0;k<3;++k) fscanf(crdreffile,"%lf",&crdref[i*3+k]);
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
      for (k=0;k<3;++k) crdrefCA[j*3+k]=crdref[i*3+k];
      ++j;
    }
  }
  fclose(crdreffile);

  numstep=mync_get_present_step_MCD(inputfilename,&nc_id);

  ncmap=GOLM_Clementi_make_native_contact(crdref,criteria,&numnc,numatom,numCAatom);
  cradii_ca=(double *)gcemalloc(sizeof(double)*numnc);
  index_nc=(int *)gcemalloc(sizeof(int)*numnc*2);

  nc=0;
  for (i=0;i<numCAatom;++i) {
    for (j=i+1;j<numCAatom;++j) {
      if (ncmap[i][j]==0) {
	index_nc[nc*2+0]=i;
	index_nc[nc*2+1]=j;
	++nc;
      }
    }
  }

  for (i=0;i<numnc;++i) {
    for (j=0;j<3;++j) {
      atom1[j]=crdrefCA[(index_nc[i*2+0])*3+j];
      atom2[j]=crdrefCA[(index_nc[i*2+1])*3+j];
    }
    length=len(atom1,atom2);
    cradii_ca[i]=length;
  }

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    //      mync_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id,crd_nc);
    mync_open_inq_get_sh_MCD(inputfilename,numCAatom,i,1,i+1,&nc_id,crd_nc);
    /*********************************************************/
    /* k=0;						       */
    /* for (j=0;j<numatom;++j) {			       */
    /* 	if (strncmp(AP.IGRAPH[j],"CA",2)==0) {	       */
    /* 	  for (l=0;l<3;++l) crd[k*3+l]=crd_nc[j][l];   */
    /* 	  ++k;					       */
    /* 	}					       */
    /* }						       */
    /*********************************************************/
    for (j=0;j<numCAatom;++j) for (k=0;k<3;++k) crd[j*3+k]=crd_nc[j][k];

    Q=0.0;
    for (j=0;j<numnc;++j) {
      for (k=0;k<3;++k) {
	atom1[k]=crd[(index_nc[j*2+0])*3+k];
	atom2[k]=crd[(index_nc[j*2+1])*3+k];
      }
      length=len(atom1,atom2);
      if (length < 1.2*cradii_ca[j]) Q+=1.0;
    }
    Q=Q/numnc;
    
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

