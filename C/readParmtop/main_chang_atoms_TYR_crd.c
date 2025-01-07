#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//#include "PTL.h"
#include "PT.h"
#include "EF.h"

#define OHO 14
#define HOO 15
#define C1O 16
#define H1O 17
#define C2O 18
#define H2O 19

#define OHN 18
#define HON 19
#define C1N 14
#define H1N 15
#define C2N 16
#define H2N 17

#define NATOMCH 6

int USAGE(char *progname);

int TYRAOLD[NATOMCH]={14,15,16,17,18,19};
int TYRANEW[NATOMCH]={18,19,14,15,16,17};

int CHON1[NATOMCH]={16,17,18,19,14,15};
int CHON2[NATOMCH]={14,15,16,17,18,19};

int NATEX_new[NATOMCH][10]={{15,16,17,18,19},{16,17,18,19},{17,18},{0},{0},{0}};
int NUMEX_new[NATOMCH]={5,4,2,1,1,1};

int main(int argc, char *argv[]) {
  int i,j,k,numatomini,numatomfin;
  int num,num_new;

  char IGRAPH_temp[NATOMCH][4];
  double CHRG_temp[NATOMCH];
  double AMASS_temp[NATOMCH];
  int IAC_temp[NATOMCH];
  int NUMEX_temp[NATOMCH];
  int *NATEX_temp;

  char ISYMBL_temp[NATOMCH][4];
  char ITREE_temp[NATOMCH][4];
  int JOIN_temp[NATOMCH];
  int IROTAT_temp[NATOMCH];

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *parmfilenamein,*parmfilenameout,*progname;
  FILE *parmfilein,*parmfileout;

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  parmfilenamein  = *argv;
  parmfilenameout = *++argv;

  parmfilein=efopen(parmfilenamein,"r");
  //  readParmtopL(parmfilein);
  readParmtop(parmfilein);
  fclose(parmfilein);

  for (i=0;i<AP.NRES;++i) {
    if (strncmp(AP.LABERES[i],"TYR",3)==0) {
      numatomini=AP.IPRES[i]-1;
      numatomfin=AP.IPRES[i+1];
      for (j=0;j<NATOMCH;++j) {
	strncpy(IGRAPH_temp[j],AP.IGRAPH[numatomini+CHON1[j]-1],4);
	CHRG_temp[j]=AP.CHRG[numatomini+CHON1[j]-1];
	AMASS_temp[j]=AP.AMASS[numatomini+CHON1[j]-1];
	IAC_temp[j]=AP.IAC[numatomini+CHON1[j]-1];
	strncpy(ISYMBL_temp[j],AP.ISYMBL[numatomini+CHON1[j]-1],4);
	strncpy(ITREE_temp[j],AP.ITREE[numatomini+CHON1[j]-1],4);
	JOIN_temp[j]=AP.JOIN[numatomini+CHON1[j]-1];
	IROTAT_temp[j]=AP.IROTAT[numatomini+CHON1[j]-1];
      }

      for (j=0;j<NATOMCH;++j) {
	strncpy(AP.IGRAPH[numatomini+CHON2[j]-1],IGRAPH_temp[j],4);
	AP.CHRG[numatomini+CHON2[j]-1]=CHRG_temp[j];
	AP.AMASS[numatomini+CHON2[j]-1]=AMASS_temp[j];
	AP.IAC[numatomini+CHON2[j]-1]=IAC_temp[j];
	strncpy(AP.ISYMBL[numatomini+CHON2[j]-1],ISYMBL_temp[j],4);
	strncpy(AP.ITREE[numatomini+CHON2[j]-1],ITREE_temp[j],4);
	AP.JOIN[numatomini+CHON2[j]-1]=JOIN_temp[j];
	AP.IROTAT[numatomini+CHON2[j]-1]=IROTAT_temp[j];
      }
    }
  }

  for (i=0;i<AP.NBONH;++i) {
    for (j=0;j<2;++j) {
      if (AP.BH[i][j]==((OHO+numatomini-1)-1)*3) AP.BH[i][j]=((OHN+numatomini-1)-1)*3;
      else if (AP.BH[i][j]==((HOO+numatomini-1)-1)*3) AP.BH[i][j]=((HON+numatomini-1)-1)*3;
      else if (AP.BH[i][j]==((C1O+numatomini-1)-1)*3) AP.BH[i][j]=((C1N+numatomini-1)-1)*3;
      else if (AP.BH[i][j]==((H1O+numatomini-1)-1)*3) AP.BH[i][j]=((H1N+numatomini-1)-1)*3;
      else if (AP.BH[i][j]==((C2O+numatomini-1)-1)*3) AP.BH[i][j]=((C2N+numatomini-1)-1)*3;
      else if (AP.BH[i][j]==((H2O+numatomini-1)-1)*3) AP.BH[i][j]=((H2N+numatomini-1)-1)*3;
    }
  }
  
  for (i=0;i<AP.MBONA;++i) {
    for (j=0;j<2;++j) {
      if (AP.BA[i][j]==((OHO+numatomini-1)-1)*3) AP.BA[i][j]=((OHN+numatomini-1)-1)*3;
      else if (AP.BA[i][j]==((HOO+numatomini-1)-1)*3) AP.BA[i][j]=((HON+numatomini-1)-1)*3;
      else if (AP.BA[i][j]==((C1O+numatomini-1)-1)*3) AP.BA[i][j]=((C1N+numatomini-1)-1)*3;
      else if (AP.BA[i][j]==((H1O+numatomini-1)-1)*3) AP.BA[i][j]=((H1N+numatomini-1)-1)*3;
      else if (AP.BA[i][j]==((C2O+numatomini-1)-1)*3) AP.BA[i][j]=((C2N+numatomini-1)-1)*3;
      else if (AP.BA[i][j]==((H2O+numatomini-1)-1)*3) AP.BA[i][j]=((H2N+numatomini-1)-1)*3;
    }
  }

  for (i=0;i<AP.NTHETH;++i) {
    for (j=0;j<3;++j) {
      if (AP.TH[i][j]==((OHO+numatomini-1)-1)*3) AP.TH[i][j]=((OHN+numatomini-1)-1)*3;
      else if (AP.TH[i][j]==((HOO+numatomini-1)-1)*3) AP.TH[i][j]=((HON+numatomini-1)-1)*3;
      else if (AP.TH[i][j]==((C1O+numatomini-1)-1)*3) AP.TH[i][j]=((C1N+numatomini-1)-1)*3;
      else if (AP.TH[i][j]==((H1O+numatomini-1)-1)*3) AP.TH[i][j]=((H1N+numatomini-1)-1)*3;
      else if (AP.TH[i][j]==((C2O+numatomini-1)-1)*3) AP.TH[i][j]=((C2N+numatomini-1)-1)*3;
      else if (AP.TH[i][j]==((H2O+numatomini-1)-1)*3) AP.TH[i][j]=((H2N+numatomini-1)-1)*3;
    }
  }

  for (i=0;i<AP.MTHETA;++i) {
    for (j=0;j<3;++j) {
      if (AP.TA[i][j]==((OHO+numatomini-1)-1)*3) AP.TA[i][j]=((OHN+numatomini-1)-1)*3;
      else if (AP.TA[i][j]==((HOO+numatomini-1)-1)*3) AP.TA[i][j]=((HON+numatomini-1)-1)*3;
      else if (AP.TA[i][j]==((C1O+numatomini-1)-1)*3) AP.TA[i][j]=((C1N+numatomini-1)-1)*3;
      else if (AP.TA[i][j]==((H1O+numatomini-1)-1)*3) AP.TA[i][j]=((H1N+numatomini-1)-1)*3;
      else if (AP.TA[i][j]==((C2O+numatomini-1)-1)*3) AP.TA[i][j]=((C2N+numatomini-1)-1)*3;
      else if (AP.TA[i][j]==((H2O+numatomini-1)-1)*3) AP.TA[i][j]=((H2N+numatomini-1)-1)*3;
    }
  }

  for (i=0;i<AP.NPHIH;++i) {
    for (j=0;j<4;++j) {
      if (AP.PH[i][j]==((OHO+numatomini-1)-1)*3) AP.PH[i][j]=((OHN+numatomini-1)-1)*3;
      else if (AP.PH[i][j]==((HOO+numatomini-1)-1)*3) AP.PH[i][j]=((HON+numatomini-1)-1)*3;
      else if (AP.PH[i][j]==((C1O+numatomini-1)-1)*3) AP.PH[i][j]=((C1N+numatomini-1)-1)*3;
      else if (AP.PH[i][j]==((H1O+numatomini-1)-1)*3) AP.PH[i][j]=((H1N+numatomini-1)-1)*3;
      else if (AP.PH[i][j]==((C2O+numatomini-1)-1)*3) AP.PH[i][j]=((C2N+numatomini-1)-1)*3;
      else if (AP.PH[i][j]==((H2O+numatomini-1)-1)*3) AP.PH[i][j]=((H2N+numatomini-1)-1)*3;
    }
  }

  for (i=0;i<AP.MPHIA;++i) {
    for (j=0;j<4;++j) {
      if (AP.PA[i][j]==((OHO+numatomini-1)-1)*3) AP.PA[i][j]=((OHN+numatomini-1)-1)*3;
      else if (AP.PA[i][j]==((HOO+numatomini-1)-1)*3) AP.PA[i][j]=((HON+numatomini-1)-1)*3;
      else if (AP.PA[i][j]==((C1O+numatomini-1)-1)*3) AP.PA[i][j]=((C1N+numatomini-1)-1)*3;
      else if (AP.PA[i][j]==((H1O+numatomini-1)-1)*3) AP.PA[i][j]=((H1N+numatomini-1)-1)*3;
      else if (AP.PA[i][j]==((C2O+numatomini-1)-1)*3) AP.PA[i][j]=((C2N+numatomini-1)-1)*3;
      else if (AP.PA[i][j]==((H2O+numatomini-1)-1)*3) AP.PA[i][j]=((H2N+numatomini-1)-1)*3;
    }
  }

  for (i=0;i<AP.NNB;++i) {
    if (AP.NATEX[i]==OHO+numatomini) 
      AP.NATEX[i]=OHN+numatomini;
    else if (AP.NATEX[i]==HOO+numatomini) 
      AP.NATEX[i]=HON+numatomini;
    else if (AP.NATEX[i]==C1O+numatomini) 
      AP.NATEX[i]=C1N+numatomini;
    else if (AP.NATEX[i]==H1O+numatomini) 
      AP.NATEX[i]=H1N+numatomini;
    else if (AP.NATEX[i]==C2O+numatomini) 
      AP.NATEX[i]=C2N+numatomini;
    else if (AP.NATEX[i]==H2O+numatomini) 
      AP.NATEX[i]=H2N+numatomini;
  }

  NATEX_temp=(int *)gcemalloc(sizeof(int)*AP.NNB+2);

  num=0;
  num_new=0;
  for (i=0;i<AP.NRES;++i) {
    if (strncmp(AP.LABERES[i],"TYR",3)==0) {
      numatomini=AP.IPRES[i]-1;
      numatomfin=AP.IPRES[i+1];
      for (j=numatomini;j<numatomfin;++j) {
	if (j==CHON2[0]+numatomini-1) {
	  for (k=0;k<NUMEX_new[0];++k) {
	    NATEX_temp[num_new+k]=NATEX_new[0][k]+numatomini;
	  }
	  NUMEX_temp[j]=NUMEX_new[0];
	  num+=AP.NUMEX[j];
	  num_new+=NUMEX_temp[j];
	}
	else if (j==CHON2[1]+numatomini-1) {
	  for (k=0;k<NUMEX_new[1];++k) {
	    NATEX_temp[num_new+k]=NATEX_new[1][k]+numatomini;
	  }
	  NUMEX_temp[j]=NUMEX_new[1];
	  num+=AP.NUMEX[j];
	  num_new+=NUMEX_temp[j];
	}
	else if (j==CHON2[2]+numatomini-1) {
	  for (k=0;k<NUMEX_new[2];++k) {
	    NATEX_temp[num_new+k]=NATEX_new[2][k]+numatomini;
	  }
	  NUMEX_temp[j]=NUMEX_new[2];
	  num+=AP.NUMEX[j];
	  num_new+=NUMEX_temp[j];
	}
	else if (j==CHON2[3]+numatomini-1) {
	  for (k=0;k<NUMEX_new[3];++k) {
	    NATEX_temp[num_new+k]=NATEX_new[3][k]+numatomini;
	  }
	  NUMEX_temp[j]=NUMEX_new[3];
	  num+=AP.NUMEX[j];
	  num_new+=NUMEX_temp[j];
	}
	else if (j==CHON2[4]+numatomini-1) {
	  for (k=0;k<NUMEX_new[4];++k) {
	    NATEX_temp[num_new+k]=NATEX_new[4][k]/*+numatomini*/;
	  }
 	  NUMEX_temp[j]=NUMEX_new[4];
	  num+=AP.NUMEX[j];
	  num_new+=NUMEX_temp[j];
	}
	else if (j==CHON2[5]+numatomini-1) {
	  for (k=0;k<NUMEX_new[5];++k) {
	    NATEX_temp[num_new+k]=NATEX_new[5][k]/*+numatomini*/;
	  }
	  NUMEX_temp[j]=NUMEX_new[5];
	  num+=AP.NUMEX[j];
	  num_new+=NUMEX_temp[j];
	}
	else {
	  for (k=0;k<AP.NUMEX[j];++k) {
	    NATEX_temp[num_new+k]=AP.NATEX[num+k];
	  }
	  NUMEX_temp[j]=AP.NUMEX[j];
	  num+=AP.NUMEX[j];
	  num_new+=AP.NUMEX[j];
	}
      }
    }
    else {
      numatomini=AP.IPRES[i]-1;
      if (i!=AP.NRES-1) {
	numatomfin=AP.IPRES[i+1];
      }
      else {
	numatomfin=AP.NATOM;
      }
      for (j=numatomini;j<numatomfin;++j) {
	for (k=0;k<AP.NUMEX[j];++k) {
	  NATEX_temp[num_new+k]=AP.NATEX[num+k];
	}
	NUMEX_temp[j]=AP.NUMEX[j];
	num+=AP.NUMEX[j];
	num_new+=AP.NUMEX[j];
      }
    }
  }

  AP.NNB+=2;
  AP.NEXT+=2;
  //  AP.NATEX=(int *)gcrealloc(AP.NATEX,sizeof(int)*AP.NNB);
  for (i=0;i<AP.NATOM;++i) {
    AP.NUMEX[i]=NUMEX_temp[i];
  }
  for (i=0;i<AP.NNB;++i) {
    AP.NATEX[i]=NATEX_temp[i];
  }

  parmfileout=efopen(parmfilenameout,"w");

  writeParmtop(parmfileout);

  fclose(parmfileout);

}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] crdin crdout \n",progname);
}
