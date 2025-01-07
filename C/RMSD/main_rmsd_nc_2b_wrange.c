#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "rmsd.h"
#include "PTL.h"
#include "EF.h"
#include "PDB.h"

#include "netcdf_mine.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int readdata(char *inputfilename,double *crd,double **crdA,PDBF PDB,
	     int numatom,int numiniatom,int numfinatom, int numatomwrange,
	     double *mass,int pdbflag,int MODE);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int numatom,numatomp,numinires,numfinres,numiniatom,numfinatom,numatomwrange;
  int numatomref,numatompref,numiniresref,numfinresref,numiniatomref,numfinatomref,numatomwrangeref;

  double rmsd,*rot;
  double *mass;

  PDBF PDB,PDBref;

  int MODE=AA,pdbflag=OFF,normflag=OFF;
  
  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double *crd,*crdref;
  double **crdA,**crdB;
  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;

  char *inputfilename,*parmfilename,*parmfilenameref,*reffilename;
  char *progname;
  FILE *inputfile,*parmfile,*reffile;

  int opt_idx=1;

  struct option long_opt[] = {
    {"ca",0,NULL,'c'},
    {"n",0,NULL,'n'},
    {"ha",0,NULL,'h'},
    {"aa",0,NULL,'a'},
    {"pdb",0,NULL,'p'},
    {"H",0,NULL,'H'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"chapnH",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'a':
      MODE=AA;
      break;
    case 'h':
      MODE=HV;
      break;
    case 'c':
      MODE=CA;
      break;
    case 'p':
      pdbflag=ON;
      break;
    case 'n':
      normflag=ON;
      break;
    case 'H':
      USAGE(progname);
      exit(1);
    }
  }

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 8) {
    USAGE(progname);
    exit(1);
  }
  inputfilename = *argv;
  reffilename = *++argv;
  parmfilename = *++argv;
  parmfilenameref = *++argv;
  numinires = atoi(*++argv);
  numfinres = atoi(*++argv);
  numiniresref = atoi(*++argv);
  numfinresref = atoi(*++argv);

  parmfile=efopen(parmfilename,"r");  
  readParmtopL(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  if (pdbflag==ON) {
    PDB.PDBa=(PDBA *)gcemalloc(sizeof(PDBA)*numatom);
  }
  numiniatom=AP.IPRES[numinires-1]-1;
  numfinatom=AP.IPRES[numfinres]-2;
  numatomwrange=(numfinatom-numiniatom+1);
  crdA=(double **)gcemalloc(sizeof(double *)*3);
  for (j=0;j<3;++j) crdA[j]=(double *)gcemalloc(sizeof(double)*numatomwrange);
  rot=(double *)gcemalloc(sizeof(double)*9);
  mass=(double *)gcemalloc(sizeof(double)*numatomwrange);

  numatomp=readdata(inputfilename,crd,crdA,PDB,
		    numatom,numiniatom,numfinatom,numatomwrange,
		    mass,pdbflag,MODE);

  parmfile=efopen(parmfilenameref,"r");
  readParmtopL(parmfile);
  fclose(parmfile);
  numatomref=AP.NATOM;
  crdref=(double *)gcemalloc(sizeof(double)*numatomref*3);
  if (pdbflag==ON) {
    PDBref.PDBa=(PDBA *)gcemalloc(sizeof(PDBA)*numatomref);
  }
  numiniatomref=AP.IPRES[numiniresref-1]-1;
  numfinatomref=AP.IPRES[numfinresref]-2;
  numatomwrangeref=(numfinatomref-numiniatomref+1);
  crdB=(double **)gcemalloc(sizeof(double *)*3);
  for (j=0;j<3;++j) crdB[j]=(double *)gcemalloc(sizeof(double)*numatomwrangeref);
  rot=(double *)gcemalloc(sizeof(double)*9);
  mass=(double *)gcemalloc(sizeof(double)*numatomwrangeref);

  numatompref=readdata(reffilename,crdref,crdB,PDBref,
		       numatomref,numiniatomref,numfinatomref,numatomwrangeref,
		       mass,pdbflag,MODE);

  if (numatomp!=numatompref || numatomp<0 || numatompref<0 ) {
    printf("error of res num\n");
  }
  
  rmsd=CalcRMSDRotationalMatrix(crdB,crdA,numatomp,rot,mass);
  printf("%10.8e",rmsd);
  if (normflag==ON)  printf(" %10.8e",rmsd/(numfinres-numinires+1));
  printf("\n");
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-ca]\n");
  printf("[-ha]\n");
  printf("[-aa]\n");
  printf("[-pdb]\n");
  printf("[-help]\n");
  printf("%s inputfilename reffilename parmfilename numinires numfinres numiniresref numfinresref\n",progname);
}

int readdata(char *inputfilename,double *crd,double **crdA,PDBF PDB,
	     int numatom,int numiniatom,int numfinatom, int numatomwrange,
	     double *mass,int pdbflag,int MODE) {
  int i,j,k,l;
  int numatomp;
  FILE *inputfile;


  inputfile=efopen(inputfilename,"r");
  if (pdbflag==OFF) {
    io_scanconf_Amber_ini(inputfile,numatom,crd);
    fclose(inputfile);
  }
  else {
    readPDB(inputfile,PDB,numatom);
    for (i=0;i<numatom;++i) {
      for (j=0;j<3;++j) {
	crd[i*3+j]=PDB.PDBa[i].coord[j];
      }
    }
  }

  l=0;
  for (j=numiniatom;j<numfinatom;++j) {
    if (MODE==CA) {
      if (strncmp(AP.IGRAPH[j],"CA",2)==0) {
	mass[l]=AP.AMASS[j];
	for (k=0;k<3;++k) {
	  crdA[k][l]=crd[j*3+k];
	}
	++l;
      }
    }
    else if (MODE==HV) {
      if (strncmp(AP.IGRAPH[j],"H",1)!=0) {
	mass[l]=AP.AMASS[j];
	for (k=0;k<3;++k) {
	  crdA[k][l]=crd[j*3+k];
	}
	++l;
      }
    }
    else {
      for (k=0;k<3;++k) {
	crdA[k][l]=crd[j*3+k];
      }
      ++l;
    }
  }
  numatomp=l;

  return numatomp;
}

