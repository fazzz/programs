#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <getopt.h>

#include "PDBL.h"

#include "TOPO.h"
#include "PT.h"
#include "EF.h"

#define ON 1
#define OFF 0
#define SUM 1
#define SERIES 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,m,d;
  int numstep,numatom,interval=1;

  int IHflag=OFF,MODE=AA,INMODE=AA;
  int OUTMODE;

  PDBLF PDBL;

  char *line;
  size_t len=0;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double *crd;

  char *trjfilename,*outputfilename,*parmfilename,*progname,*logfilename;
  char outputfilename2[10000];
  FILE *trjfile,*pdbfile,*parmfile,*outputfile,*logfile;

  int opt_idx=1;

  struct option long_opt[] = {
    {"int",1,NULL,'i'},
    {"PDBs",1,NULL,'P'},
    {"CA",0,NULL,'C'},
    {"IMCA",0,NULL,'I'},
    {"IMHV",0,NULL,'J'},
    {"H",0,NULL,'H'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hCHIJHi:P:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'i':
      interval=atoi(optarg);
      break;
    case 'C':
      MODE=CA;
      break;
    case 'I':
      INMODE=CA;
      break;
    case 'J':
      INMODE=HV;
      break;
    case 'H':
      MODE=HV;
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    case 'P':
      logfilename=optarg;
      OUTMODE=SERIES;
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
  numstep = atoi(*argv);
  trjfilename  = *++argv;
  parmfilename = *++argv;
  outputfilename  = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;
  if (INMODE==CA) {
    j=0;
    for (i=0;i<numatom;++i) {
      if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
	++j;
      }
    }
    numatom=j;
  }
  else if (INMODE==HV) {
    j=0;
    for (i=0;i<numatom;++i) {
      if (strncmp(AP.IGRAPH[i],"H",1)!=0) {
	++j;
      }
    }
    numatom=j;
  }

  PDBL.numatom=numatom;
  PDBL.PDBLa=(PDBLA *)gcemalloc(sizeof(PDBLA)*numatom);
  readPDBLdatafromParmtop(PDBL);

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  if (OUTMODE!=SERIES) {
    sprintf(outputfilename2,"%s.pdb",outputfilename);
    outputfile=efopen(outputfilename2,"w");
  }
  else if (OUTMODE==SERIES) logfile=efopen(logfilename,"w");

  trjfile=efopen(trjfilename,"r");
  getline(&line,&len,trjfile);

  for (i=0;i<numstep;++i) {
    for (j=0;j<numatom;++j) for (k=0;k<3;++k) fscanf(trjfile,"%lf",&crd[j*3+k]);

    if (i%interval == 0) {
      if (INMODE!=CA) {
	for (j=0;j<AP.NATOM;++j) for (k=0;k < 3; ++k) PDBL.PDBLa[j].coord[k]=crd[j*3+k];
      }
      else if (INMODE==CA) {
	m=0;
	for (j=0;j<AP.NATOM;++j) {
	  if (strncmp(PDBL.PDBLa[j].name," CA",3)==0) {
	    for (k=0;k < 3; ++k) PDBL.PDBLa[j].coord[k]=crd[m*3+k];
	    ++m;
	  }
	  else {
	    for (k=0;k < 3; ++k) PDBL.PDBLa[j].coord[k]=0.0;
	  }
	}
      }
      else if (INMODE==HV) {
	m=0;
	for (j=0;j<AP.NATOM;++j) {
	  if (strncmp(PDBL.PDBLa[j].name," H",2)!=0) {
	    for (k=0;k < 3; ++k) PDBL.PDBLa[j].coord[k]=crd[m*3+k];
	    ++m;
	  }
	  else {
	    for (k=0;k < 3; ++k) PDBL.PDBLa[j].coord[k]=0.0;
	  }
	}
      }

      if (OUTMODE!=SERIES) {
	writPDBL_wopt_series(outputfile,PDBL,MODE);
      }
      else if (OUTMODE==SERIES) {
	++l;
	sprintf(outputfilename2,"%s_%d.pdb",outputfilename,l);
	outputfile=efopen(outputfilename2,"w");
	writPDBL_wopt(outputfile,PDBL,MODE);
	fclose(outputfile);
	fprintf(logfile,"%s\n",outputfilename2);
      }
    }
  }
  if (OUTMODE!=SERIES) fclose(outputfile);
  else fclose(logfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-H] (include Hydrogen) \n");
  printf("[--CA] (include Hydrogen) \n");
  printf("[--PDBs] logfilename \n");
  printf("[--int] interval \n");
  printf("[-h] help \n");
  printf("%s [-h] numstep trjfilename parmfilename pdbfilename \n",progname);
}
