
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "NC.h"
#include "PDB.h"
#include "PT.h"
#include "netcdf_mine.h"
#include "EF.h"
#include "IO.h"

#define ON 1
#define OFF 0

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,t;
  int num;
  int pdbflag=OFF,AAMODEflag=OFF;

  int numnc,*indexncb;
  int **ncmap;
  double *crd;
  int numatom,numres;
  double criteria=criteria_NC;

  PDBF PDB;

  char *progname;
  char *crdfilename,*parmtopname;
  char *ncmapfilename,*nclistfilename;

  FILE *crdfile,*parmtop;
  FILE *ncmapfile,*nclistfile;

  int ncid,fene_varid;
  size_t start[1],count[1];

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"hapc:"))!=-1) {
    switch(c) {
    case 'a':
      AAMODEflag=ON;
      break;
    case 'c':
      criteria=atof(optarg);
      break;
    case 'p':
      pdbflag=ON;
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  argc-=optind;
  argv+=optind;

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  crdfilename=*argv;
  parmtopname=*++argv;
  ncmapfilename=*++argv;
  nclistfilename=*++argv;

  parmtop=efopen(parmtopname,"r");
  readParmtopL(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;
  numres=AP.NRES;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  if (pdbflag==ON) PDB.PDBa=(PDBA *)gcemalloc(sizeof(PDBA)*numatom);

  crdfile=efopen(crdfilename,"r");
  if (pdbflag==OFF) {
    io_scanconf_Amber_ini(crdfile,numatom,crd);
    fclose(crdfile);
  }
  else {
    readPDB(crdfile,PDB,numatom);
    for (i=0;i<numatom;++i) 
      for (j=0;j<3;++j) 
	crd[i*3+j]=PDB.PDBa[i].coord[j];
  }

  if (AAMODEflag==OFF) num=numres;
  else num=numatom;

  ncmap=(int **)gcemalloc(sizeof(int *)*numatom);
  for (i=0;i<numatom++i) ncmap[i]=(int *)gcemalloc(sizeof(int)*numatom);

  if (AAMODEflag==OFF)
    indexncb=make_native_contact_list(&numnc,crd,numatom,numres,criteria,ncmap);
  else
    indexncb=make_native_contact_list_aa(&numnc,crd,numatom,criteria,ncmap,OFF);

  ncmapfile=efopen(ncmapfilename,"w");
  for (i=0;i<num;++i) {
    for (j=0;j<num;++j)
      fprintf(ncmapfile,"%2d ",ncmap[i*2+j]);
    fprintf(ncmapfile,"\n");
  }
  fclose(ncmapfile);

  nclistfile=efopen(nclistfilename,"w");
  for (i=0;i<numnc;++i)
      fprintf(nclistfile,"%d-%d\n",indexncb[i*2],indexncb[i*2+1]);
  fclose(nclistfile);

  return 0;

}

void USAGE(char *progname) {
  printf("[-h] -- help\n");
  printf("USAGE: %s  crdfilename parmtopname ncmapname nclistmap \n", progname);
}
