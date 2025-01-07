#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PDBL.h"

#include "TOPO.h"
#include "PT.h"
#include "EF.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0
#define SUM 1
#define SERIES 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,m;
  int numstep,numatom,interval=1;

  int IHflag=OFF,AMBncflag=OFF,MODE=AA,INMODE=AA;
  int OUTMODE;

  PDBLF PDBL;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double *crd;
  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id;

  char *inputfilename,*outputfilename,*parmfilename,*progname,*logfilename;
  char outputfilename2[10000];
  FILE *outputfile,*parmfile,*logfile;

  int opt_idx=1;

  struct option long_opt[] = {
    {"int",1,NULL,'i'},
    {"PDBs",1,NULL,'P'},
    {"CA",0,NULL,'C'},
    {"IMCA",0,NULL,'I'},
    {"IMHV",0,NULL,'J'},
    {"H",0,NULL,'H'},
    {"Amber",0,NULL,'K'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hCIJHKi:P:",long_opt,&opt_idx))!=-1) {
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
    case 'P':
      logfilename=optarg;
      OUTMODE=SERIES;
      break;
    case 'K':
      AMBncflag=ON;
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
  inputfilename = *argv;
  parmfilename = *++argv;
  outputfilename = *++argv;

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

  PDBL.numatom=AP.NATOM;
  PDBL.PDBLa=(PDBLA *)gcemalloc(sizeof(PDBLA)*AP.NATOM);
  readPDBLdatafromParmtop(PDBL);

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  if (AMBncflag==ON) numstep=myncL_get_present_step_AMBER(inputfilename,&nc_id);
  else numstep=myncL_get_present_step_MCD(inputfilename,&nc_id_MCD);

  if (OUTMODE!=SERIES) {
    sprintf(outputfilename2,"%s.pdb",outputfilename);
    outputfile=efopen(outputfilename2,"w");
  }
  else if (OUTMODE==SERIES) logfile=efopen(logfilename,"w");

  for (i=0;i<numstep;++i) {
    if (AMBncflag==ON) myncL_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id,crd_nc);
    else myncL_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd_nc);

    if (i%interval == 0) {
      if (INMODE!=CA) {
	for (j=0;j<AP.NATOM;++j) for (k=0;k < 3; ++k) PDBL.PDBLa[j].coord[k]=crd_nc[j][k];
      }
      else if (INMODE==CA) {
	m=0;
	for (j=0;j<AP.NATOM;++j) {
	  if (strncmp(PDBL.PDBLa[j].name," CA",3)==0) {
	    for (k=0;k < 3; ++k) PDBL.PDBLa[j].coord[k]=crd_nc[m][k];
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
	    for (k=0;k < 3; ++k) PDBL.PDBLa[j].coord[k]=crd_nc[m][k];
	    ++m;
	  }
	  else {
	    for (k=0;k < 3; ++k) PDBL.PDBLa[j].coord[k]=0.0;
	  }
	}
      }

      if (OUTMODE!=SERIES) {
	//	fprintf(outputfile,"MODEL \n");
	//    writPDBL(outputfile,PDBL);
	writPDBL_wopt_series(outputfile,PDBL,MODE);
	//	fprintf(outputfile,"ENDMOD\n");
	//	fprintf(outputfile,"MODEL \n");
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
  printf("[--Amber] (amber) \n");
  printf("[--int] interval \n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename parmfilename outputfilename \n",progname);
}
