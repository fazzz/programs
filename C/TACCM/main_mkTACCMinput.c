#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "TACCM.h"
#include "ABAb.h"

#include "d2rc.h"

#include "PTL.h"
#include "FFL.h"
#include "EF.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,dummyd;
  int joinflag;
  int numatom,numheavyatom,numres;
  int specifymove=OFF;
  int numclut,*numdihed;

  int Omegaflag=OFF,Kaiflag=OFF;

  double *Z;
  int numZ;

  int *resid;
  int **adpairs;
  int *inpindex,inpnumH,inpnumA,*indexclut;
  int *numclutparent,*terminal,*origin;
  CLTb *clt;

  double pi;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*TACCMfilename,*parmfilename,*clustfilename,*logfilename="mkTACCMinput.log";
  FILE *inputfile,*TACCMfile,*parmfile,*clustfile,*logfile;

  char *progname;
  int opt_idx=1;

  struct option long_opt[] = {
    {"Omega",0,NULL,'O'},
    {"Kai",0,NULL,'K'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hOK",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'O':
      Omegaflag=ON;
      break;
    case 'K':
      Kaiflag=ON;
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

  pi=acos(-1.0);

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  inputfilename  = *argv;
  clustfilename  = *++argv;
  parmfilename   = *++argv;
  TACCMfilename  = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;

  clustfile=efopen(clustfilename,"r");
  clt=ABAbp_clustscan(clustfile,&numclut);
  fclose(clustfile);

  inputfile=efopen(inputfilename,"r");
  resid=readd2rcinput(inputfile,&numres,specifymove);
  fclose(inputfile);

  numclutparent=(int *)gcemalloc(sizeof(int)*numclut);
  terminal=(int *)gcemalloc(sizeof(int)*numclut);
  origin=(int *)gcemalloc(sizeof(int)*numclut);
  for (i=0;i<numclut;++i) {
    numclutparent[i]=clt[i].nNumClutOfParent;
    terminal[i]=clt[i].terminal_atom_a[0];
    origin[i]=clt[i].origin_atom_a;
  }

  adpairs=(int **)gcemalloc(sizeof(int *)*5);
  adpairs[0]=(int *)gcemalloc(sizeof(int)*4);
  adpairs[1]=(int *)gcemalloc(sizeof(int)*4);
  adpairs[2]=(int *)gcemalloc(sizeof(int)*4);
  adpairs[3]=(int *)gcemalloc(sizeof(int)*8);
  adpairs[4]=(int *)gcemalloc(sizeof(int)*8);
  numdihed=(int *)gcemalloc(sizeof(int)*5);  
  readdihedpairsL(adpairs,numdihed);

  indexclut=(int *)gcemalloc(sizeof(int)*(AP.NPHIH+AP.MPHIA));
  inpindex=ffL_make_inpindex(&inpnumH,&inpnumA,indexclut,numclut,numclutparent,terminal,origin);

  TACCMfile=efopen(TACCMfilename,"w");
  numZ=0;
  for (i=0;i<numdihed[0];++i) {
    if (      PTL_which_include(PTL_resnum2(adpairs[0][i*4+1]),resid,numres)==0 
	   && PTL_which_include(PTL_resnum2(adpairs[0][i*4+2]),resid,numres)==0) {
      ++numZ;
    }
  }
  if (Omegaflag==ON) {
    for (i=0;i<numdihed[1];++i) {
      if (      PTL_which_include(PTL_resnum2(adpairs[1][i*4+1]),resid,numres)==0 
		&& PTL_which_include(PTL_resnum2(adpairs[1][i*4+2]),resid,numres)==0) {
	++numZ;
      }
    }
  }
  if (Kaiflag==ON) {
    for (i=0;i<numdihed[2];++i) {
      if (      PTL_which_include(PTL_resnum2(adpairs[2][i*4+1]),resid,numres)==0 
		&& PTL_which_include(PTL_resnum2(adpairs[2][i*4+2]),resid,numres)==0) {
	++numZ;
      }
    }
  }

  fprintf(TACCMfile,"%4d\n",numZ);
  for (i=0;i<numdihed[0];++i) {
    if (      PTL_which_include(PTL_resnum2(adpairs[0][i*4+1]),resid,numres)==0 
	   && PTL_which_include(PTL_resnum2(adpairs[0][i*4+2]),resid,numres)==0) {
      for (j=0;j<AP.MPHIA;++j) {
	if (adpairs[0][i*4+0]==abs(AP.PA[j][0])/3 && adpairs[0][i*4+1]==abs(AP.PA[j][1])/3  &&
	    adpairs[0][i*4+2]==abs(AP.PA[j][2])/3 && adpairs[0][i*4+3]==abs(AP.PA[j][3])/3  )
	  k=indexclut[j+AP.NPHIH];
      }
      for (j=0;j<4;++j) {
	fprintf(TACCMfile,"%4d ",adpairs[0][i*4+j]+1);
      }
      fprintf(TACCMfile,"%4d\n",k-1);
    }
  }
  if (Omegaflag==ON) {
    for (i=0;i<numdihed[1];++i) {
      if (      PTL_which_include(PTL_resnum2(adpairs[1][i*4+1]),resid,numres)==0 
		&& PTL_which_include(PTL_resnum2(adpairs[1][i*4+2]),resid,numres)==0) {
	for (j=0;j<AP.MPHIA;++j) {
	  if (adpairs[1][i*4+0]==abs(AP.PA[j][0])/3 && adpairs[1][i*4+1]==abs(AP.PA[j][1])/3  &&
	      adpairs[1][i*4+2]==abs(AP.PA[j][2])/3 && adpairs[1][i*4+3]==abs(AP.PA[j][3])/3  )
	    k=indexclut[j+AP.NPHIH];
	}
	for (j=0;j<4;++j) {
	  fprintf(TACCMfile,"%4d ",adpairs[1][i*4+j]+1);
	}
	fprintf(TACCMfile,"%4d\n",k-1);
      }
    }
  }
  if (Kaiflag==ON) {
    for (i=0;i<numdihed[2];++i) {
      if (      PTL_which_include(PTL_resnum2(adpairs[2][i*4+1]),resid,numres)==0 
		&& PTL_which_include(PTL_resnum2(adpairs[2][i*4+2]),resid,numres)==0) {
	for (j=0;j<AP.MPHIA;++j) {
	  if (adpairs[2][i*4+0]==abs(AP.PA[j][0])/3 && adpairs[2][i*4+1]==abs(AP.PA[j][1])/3  &&
	      adpairs[2][i*4+2]==abs(AP.PA[j][2])/3 && adpairs[2][i*4+3]==abs(AP.PA[j][3])/3  )
	    k=indexclut[j+AP.NPHIH];
	}
	for (j=0;j<4;++j) {
	  fprintf(TACCMfile,"%4d ",adpairs[2][i*4+j]+1);
	}
	fprintf(TACCMfile,"%4d\n",k-1);
      }
    }
  }
  fclose(TACCMfile);

  return 0;

}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-O] Omegaflag \n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename clustfilename parmfilename TACCMfilename \n",progname);
}
