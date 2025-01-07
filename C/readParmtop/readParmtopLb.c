#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>

#include "PTLb.h"
#include "EF.h"

int readParmtopLb(FILE *parmfile,struct AmberParmL*ap){
  int i,j;
  char *line;
  size_t len=0;

  getline(&line,&len,parmfile);
  getline(&line,&len,parmfile);
  getline(&line,&len,parmfile);
  getline(&line,&len,parmfile);
  //  fscanf(parmfile,"%4s",ap.ITITLE));

  for (i=0;i<2;++i){
    getline(&line,&len,parmfile);
  }
  fscanf(parmfile,"%8d",&(ap->NATOM));
  fscanf(parmfile,"%8d",&(ap->NTYPES));
  fscanf(parmfile,"%8d",&(ap->NBONH));
  fscanf(parmfile,"%8d",&(ap->MBONA));
  fscanf(parmfile,"%8d",&(ap->NTHETH));
  fscanf(parmfile,"%8d",&(ap->MTHETA));
  fscanf(parmfile,"%8d",&(ap->NPHIH));
  fscanf(parmfile,"%8d",&(ap->MPHIA));
  fscanf(parmfile,"%8d",&(ap->NHPARM));
  fscanf(parmfile,"%8d",&(ap->NPARM));
  fscanf(parmfile,"%8d",&(ap->NNB));
  (*ap).NEXT = (*ap).NNB;
  fscanf(parmfile,"%8d",&(ap->NRES));
  fscanf(parmfile,"%8d",&(ap->NBONA));
  fscanf(parmfile,"%8d",&(ap->NTHETA));
  fscanf(parmfile,"%8d",&(ap->NPHIA));
  fscanf(parmfile,"%8d",&(ap->NUMBND));
  fscanf(parmfile,"%8d",&(ap->NUMANG));
  fscanf(parmfile,"%8d",&(ap->NPTRA));
  fscanf(parmfile,"%8d",&(ap->NATYP));
  fscanf(parmfile,"%8d",&(ap->NPHB));
  fscanf(parmfile,"%8d",&(ap->IFPERT));
  fscanf(parmfile,"%8d",&(ap->NBPER));
  fscanf(parmfile,"%8d",&(ap->NGPER));
  fscanf(parmfile,"%8d",&(ap->NDPER));
  fscanf(parmfile,"%8d",&(ap->MBPER));
  fscanf(parmfile,"%8d",&(ap->MGPER));
  fscanf(parmfile,"%8d",&(ap->MDPER));
  fscanf(parmfile,"%8d",&(ap->IFBOX));
  fscanf(parmfile,"%8d",&(ap->NMXPS));
  fscanf(parmfile,"%8d",&(ap->IFCAP));
  fscanf(parmfile,"%8d",&(ap->NEXTRA));

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<(*ap).NATOM;++i){
    fscanf(parmfile,"%4s",&(ap->IGRAPH[i]));
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).CHRG=(double *)gcemalloc(sizeof(double)*(*ap).NATOM);
  for (i=0;i<(*ap).NATOM;++i){
    fscanf(parmfile,"%lf",&(ap->CHRG[i]));
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).AMASS=(double *)gcemalloc(sizeof(double)*(*ap).NATOM);
  for (i=0;i<(*ap).NATOM;++i){
    fscanf(parmfile,"%lf",&(ap->AMASS[i]));
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).IAC=(int/*double*/ *)gcemalloc(sizeof(int/*double*/)*(*ap).NATOM); // 0811
  for (i=0;i<(*ap).NATOM;++i){
    fscanf(parmfile,"%8d",&(ap->IAC[i]));
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).NUMEX=(int/*double*/ *)gcemalloc(sizeof(int/*double*/)*(*ap).NATOM); // 0811
  for (i=0;i<(*ap).NATOM;++i){
    fscanf(parmfile,"%8d",&(ap->NUMEX[i]));
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).ICO=(int/*double*/ *)gcemalloc(sizeof(int/*double*/)*((*ap).NTYPES)*((*ap).NTYPES)); // 0811
  for (i=0;i<((*ap).NTYPES)*((*ap).NTYPES);++i){
    fscanf(parmfile,"%8d",&(ap->ICO[i]));
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<(*ap).NRES;++i){
    fscanf(parmfile,"%4s",&(ap->LABERES[i]));
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).IPRES=(int *)gcemalloc(sizeof(int)*(*ap).NRES);
  for (i=0;i<(*ap).NRES;++i){
    fscanf(parmfile,"%8d",&(ap->IPRES[i]));
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).RK=(double *)gcemalloc(sizeof(double)*(*ap).NUMBND);
  for (i=0;i<(*ap).NUMBND;++i){
    fscanf(parmfile,"%lf",&(ap->RK[i]));
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).REQ=(double *)gcemalloc(sizeof(double)*(*ap).NUMBND);
  for (i=0;i<(*ap).NUMBND;++i){
    fscanf(parmfile,"%lf",&(ap->REQ[i]));
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).TK=(double *)gcemalloc(sizeof(double)*(*ap).NUMANG);
  for (i=0;i<(*ap).NUMANG;++i){
    fscanf(parmfile,"%lf",&(ap->TK[i]));
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).TEQ=(double *)gcemalloc(sizeof(double)*(*ap).NUMANG);
  for (i=0;i<(*ap).NUMANG;++i){
    fscanf(parmfile,"%lf",&(ap->TEQ[i]));
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).PK=(double *)gcemalloc(sizeof(double)*(*ap).NPTRA);
  for (i=0;i<(*ap).NPTRA;++i){
    fscanf(parmfile,"%lf",&(ap->PK[i]));
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).PN=(double *)gcemalloc(sizeof(double)*(*ap).NPTRA);
  for (i=0;i<(*ap).NPTRA;++i){
    fscanf(parmfile,"%lf",&(ap->PN[i]));
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).PHASE=(double *)gcemalloc(sizeof(double)*(*ap).NPTRA);
  for (i=0;i<(*ap).NPTRA;++i){
    fscanf(parmfile,"%lf",&(ap->PHASE[i]));
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).SOLTY=(double *)gcemalloc(sizeof(double)*(*ap).NATYP);
  for (i=0;i<(*ap).NATYP;++i){
    fscanf(parmfile,"%lf",&(ap->SOLTY[i]));
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).CN1=(double *)gcemalloc(sizeof(double)*((*ap).NTYPES)*((*ap).NTYPES+1)/2);
  for (i=0;i<((*ap).NTYPES)*((*ap).NTYPES+1)/2;++i){
    fscanf(parmfile,"%lf",&(ap->CN1[i]));
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).CN2=(double *)gcemalloc(sizeof(double)*((*ap).NTYPES)*((*ap).NTYPES+1)/2);
  for (i=0;i<((*ap).NTYPES)*((*ap).NTYPES+1)/2;++i){
    fscanf(parmfile,"%lf",&(ap->CN2[i]));
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).BH=(int/*double*/ **)gcemalloc(sizeof(int/*double*/ *)*(*ap).NBONH); // 0811
  for (i=0;i<(*ap).NBONH;++i) (*ap).BH[i]=(int/*double*/ *)gcemalloc(sizeof(int/*double*/)*3);
  for (i=0;i<(*ap).NBONH;++i){
    for (j=0;j<3;++j){
      fscanf(parmfile,"%8d",&(ap->BH[i][j]));
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).BA=(int/*double*/ **)gcemalloc(sizeof(int/*double*/ *)*(*ap).NBONA);
  for (i=0;i<(*ap).NBONA;++i) (*ap).BA[i]=(int/*double*/ *)gcemalloc(sizeof(int/*double*/)*3);
  for (i=0;i<(*ap).NBONA;++i){
    for (j=0;j<3;++j){
      fscanf(parmfile,"%d",&(ap->BA[i][j]));
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).TH=(int/*double*/ **)gcemalloc(sizeof(int/*double*/ *)*(*ap).NTHETH);
  for (i=0;i<(*ap).NTHETH;++i) (*ap).TH[i]=(int/*double*/ *)gcemalloc(sizeof(int/*double*/)*4);
  for (i=0;i<(*ap).NTHETH;++i){
    for (j=0;j<4;++j){
      fscanf(parmfile,"%8d",&(ap->TH[i][j]));
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).TA=(int/*double*/ **)gcemalloc(sizeof(int/*double*/ *)*(*ap).NTHETA);
  for (i=0;i<(*ap).NTHETA;++i) (*ap).TA[i]=(int/*double*/ *)gcemalloc(sizeof(int/*double*/)*4);
  for (i=0;i<(*ap).NTHETA;++i){
    for (j=0;j<4;++j){
      fscanf(parmfile,"%d",&(ap->TA[i][j]));
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).PH=(int/*double*/ **)gcemalloc(sizeof(int/*double*/ *)*(*ap).NPHIH);
  for (i=0;i<(*ap).NPHIH;++i) (*ap).PH[i]=(int/*double*/ *)gcemalloc(sizeof(int/*double*/)*5);
  for (i=0;i<(*ap).NPHIH;++i){
    for (j=0;j<5;++j){
      fscanf(parmfile,"%8d",&(ap->PH[i][j]));
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).PA=(int/*double*/ **)gcemalloc(sizeof(int/*double*/ *)*(*ap).NPHIA);
  for (i=0;i<(*ap).NPHIA;++i) (*ap).PA[i]=(int/*double*/ *)gcemalloc(sizeof(int/*double*/)*5);
  for (i=0;i<(*ap).NPHIA;++i){
    for (j=0;j<5;++j){
      fscanf(parmfile,"%8d",&(ap->PA[i][j]));
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).NATEX=(int *)gcemalloc(sizeof(int)*(*ap).NEXT);
  for (i=0;i<(*ap).NEXT;++i){
    fscanf(parmfile,"%8d",&(ap->NATEX[i]));
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).ASOL=(double *)gcemalloc(sizeof(double)*(*ap).NPHB);
  for (i=0;i<(*ap).NPHB;++i){
    fscanf(parmfile,"%lf",&(ap->ASOL[i]));
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).BSOL=(double *)gcemalloc(sizeof(double)*(*ap).NPHB);
  for (i=0;i<(*ap).NPHB;++i){
    fscanf(parmfile,"%lf",&(ap->BSOL[i]));
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).HBCUT=(double *)gcemalloc(sizeof(double)*(*ap).NPHB);
  for (i=0;i<(*ap).NPHB;++i){
    fscanf(parmfile,"%lf",&(ap->HBCUT[i]));
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<(*ap).NATOM;++i){
    fscanf(parmfile,"%4s",&(ap->ISYMBL[i]));
  }
  
  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<(*ap).NATOM;++i){
    fscanf(parmfile,"%4s",&(ap->ITREE[i]));
  }
  
  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).JOIN=(int *)gcemalloc(sizeof(int)*(*ap).NATOM);
  for (i=0;i<(*ap).NATOM;++i){
    fscanf(parmfile,"%8d",&(ap->JOIN[i]));
  }
  
  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  (*ap).IROTAT=(int *)gcemalloc(sizeof(int)*(*ap).NATOM);
  for (i=0;i<(*ap).NATOM;++i){
    fscanf(parmfile,"%8d",&(ap->IROTAT[i]));
  }
  
  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  if ((*ap).IFBOX > 0) {
    fscanf(parmfile,"%d",&(ap->IPTRES));
    
    fscanf(parmfile,"%d",&(ap->NSPM));
    
    fscanf(parmfile,"%d",&(ap->NSPSOL));
    
    (*ap).NSP=(int *)gcemalloc(sizeof(int)*(*ap).NSPM);
    for (i=0;i<(*ap).NSPM;++i){
      fscanf(parmfile,"%d",&(ap->NSP[i]));
    }
    
    fscanf(parmfile,"%lf",&(ap->BETA));
    
    for (i=0;i<3;++i){
      fscanf(parmfile,"%12d",&(ap->BOX[i]));
    }
  }
  
  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  if ((*ap).IFCAP > 0) {
    fscanf(parmfile,"%d",&(ap->NATCAP));
    
    fscanf(parmfile,"%lf",&(ap->CUTCAP));
    
    fscanf(parmfile,"%lf",&(ap->XCAP));
    
    fscanf(parmfile,"%lf",&(ap->YCAP));
    
    fscanf(parmfile,"%lf",&(ap->ZCAP));
  }
  
  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  if ((*ap).IFPERT > 0) {
    (*ap).BPER=(int **)gcemalloc(sizeof(int *)*(*ap).NBPER);
    for (i=0;i<(*ap).NBPER;++i) (*ap).BPER[i]=(int *)gcemalloc(sizeof(int)*2);
    for (i=0;i<(*ap).NBPER;++i){
      for (j=0;j<2;++j){
	fscanf(parmfile,"%d",&(ap->BPER[i][j]));
      }
    }
    
    (*ap).ICBPER=(int *)gcemalloc(sizeof(int)*(*ap).NBPER*1);
    for (i=0;i<(*ap).NBPER*2;++i){
      fscanf(parmfile,"%d",&(ap->ICBPER[i]));
    }
    
    (*ap).TPER=(int **)gcemalloc(sizeof(int *)*(*ap).NGPER);
    for (i=0;i<(*ap).NGPER;++i) (*ap).TPER[i]=(int *)gcemalloc(sizeof(int)*3);
    for (i=0;i<(*ap).NGPER;++i){
      for (j=0;j<3;++j){
	fscanf(parmfile,"%d",&(ap->TPER[i][j]));
      }
    }
    
    (*ap).ICTPER=(int *)gcemalloc(sizeof(int)*(*ap).NGPER*2);
    for (i=0;i<(*ap).NGPER*2;++i){
      fscanf(parmfile,"%d",&(ap->ICTPER[i]));
    }

    (*ap).PPER=(int **)gcemalloc(sizeof(int *)*(*ap).NDPER);
    for (i=0;i<(*ap).NDPER;++i)   (*ap).PPER=(int *)gcemalloc(sizeof(int)*4);
    for (i=0;i<(*ap).NDPER;++i) {
      for (j=0;j<4;++j){
	fscanf(parmfile,"%d",&(ap->PPER[i][j]));
      }
    }

    (*ap).ICTPER=(int *)gcemalloc(sizeof(int)*(*ap).NDPER*2);
    for (i=0;i<(*ap).NDPER*2;++i){
      fscanf(parmfile,"%d",&(ap->ICPPER[i]));
    }

    for (i=0;i<(*ap).NRES;++i){
      fscanf(parmfile,"%4s",&(ap->LABERES[i]));
    }

    for (i=0;i<(*ap).NATOM;++i){
      fscanf(parmfile,"%4s",&(ap->IGRPER[i]));
    }

    for (i=0;i<(*ap).NATOM;++i){
      fscanf(parmfile,"%4s",&(ap->ISMPER[i]));
    }

    for (i=0;i<(*ap).NATOM;++i){
      fscanf(parmfile,"%lf",&(ap->ALMPER[i]));
    }

    (*ap).IAPER=(int *)gcemalloc(sizeof(int)*(*ap).NATOM);
    for (i=0;i<(*ap).NATOM;++i){
      fscanf(parmfile,"%lf",&(ap->IAPER[i]));
    }

    (*ap).IACPER=(int *)gcemalloc(sizeof(int)*(*ap).NATOM);
    for (i=0;i<(*ap).NATOM;++i){
      fscanf(parmfile,"%d",&(ap->IACPER[i]));
    }

    (*ap).CGPER=(int *)gcemalloc(sizeof(int)*(*ap).NATOM);
    for (i=0;i<(*ap).NATOM;++i){
      fscanf(parmfile,"%lf",&(ap->CGPER[i]));
    }
} 

  /* if (ap.IPOL == 1) {
    fscanf(parmfile,"%18.8e",&ap.ATPOL[i]);

    fscanf(parmfile,"%18.8e",&ap.ATPOL1[i]);
    }*/

if ((*ap).NPARM == 1) {
    fscanf(parmfile,"%d",&(ap->NLES_NTYP));

    (*ap).LES_TYPE=(int/*double*/ *)gcemalloc(sizeof(int/*double*/)*(*ap).NATOM);
    for (i=0;i<(*ap).NATOM;++i){
      fscanf(parmfile,"%lf",&(ap->LES_TYPE[i]));
    }

    (*ap).LES_FAC=(double *)gcemalloc(sizeof(double)*(*ap).NATOM);
    for (i=0;i<((*ap).NLES_NTYP)*((*ap).NLES_NTYP);++i){
      fscanf(parmfile,"%lf",&(ap->LES_FAC[i]));
    }

    (*ap).LES_CNUM=(double *)gcemalloc(sizeof(double)*(*ap).NATOM);
    for (i=0;i<(*ap).NATOM;++i){
      fscanf(parmfile,"%lf",&(ap->LES_CNUM[i]));
    }

    (*ap).LES_ID=(double *)gcemalloc(sizeof(double)*(*ap).NATOM);
    for (i=0;i<(*ap).NATOM;++i){
      fscanf(parmfile,"%lf",&(ap->LES_ID[i]));
    }

  }

}

int readdihedpairsLb(int **atomdihedpairs, int *num,struct AmberParmL ap) {
  int i,j;
  int phi[4],psi[4],omega[4],ipsi[4],iomega[4],fomega[4],fphi[4],sumnum;
  int kai[10][4];
  int PHIFLAG=-1,PSIFLAG=-1,OMEGAFLAG=-1,ACEFLAG,NMEFLAG;
  int numphsi,numomega,numkai,numres;
  char RESNAME[4];

  numphsi=0;numomega=0;numkai=0;
  
  for (i=0;i<ap.NATOM;++i) {
    if (strncmp(ap.ITREE[i],"M",1)==0) {
      if (strncmp(ap.IGRAPH[i],"CH3",3)==0){
	;
      }
      else if (strncmp(ap.IGRAPH[i],"CA",2)==0) {
	if (PHIFLAG==1) {
	  PHIFLAG=2;
	  phi[PHIFLAG]=i;
	}
	if (PSIFLAG==0) {
	  PSIFLAG=1;
	  psi[PSIFLAG]=i;
	}
	if (OMEGAFLAG==2) {
	  ++numomega;
	  atomdihedpairs[1]=(int *)gcerealloc(atomdihedpairs[1],sizeof(int)*numomega*4);
	  atomdihedpairs[1][numomega*4-4]=omega[0];
	  atomdihedpairs[1][numomega*4-3]=omega[1];
	  atomdihedpairs[1][numomega*4-2]=omega[2];
	  atomdihedpairs[1][numomega*4-1]=i;
	}
	OMEGAFLAG=0;
	omega[0]=i;
      }
      else if (strncmp(ap.IGRAPH[i],"C",1)==0){
	if (PSIFLAG==1) {
	  PSIFLAG=2;
	  psi[PSIFLAG]=i;
	}
	if (OMEGAFLAG==0) {
	  OMEGAFLAG=1;
	  omega[OMEGAFLAG]=i;
	}
	if (PHIFLAG==2) {
	  phi[3]=i;
	  ++numphsi;
	  atomdihedpairs[0]=(int *)gcerealloc(atomdihedpairs[0],sizeof(int)*numphsi*4);
	  atomdihedpairs[0][numphsi*4-4]=phi[0];
	  atomdihedpairs[0][numphsi*4-3]=phi[1];
	  atomdihedpairs[0][numphsi*4-2]=phi[2];
	  atomdihedpairs[0][numphsi*4-1]=i;
	}
	PHIFLAG=0;
	phi[0]=i;
      }
      else if (strncmp(ap.IGRAPH[i],"N",1)==0){
	if (PHIFLAG==0) {
	  PHIFLAG=1;
	  phi[PHIFLAG]=i;
	}
	if (OMEGAFLAG==1) {
	  OMEGAFLAG=2;
	  omega[OMEGAFLAG]=i;
	}
	if (PSIFLAG==2) {
	  psi[3]=i;
 	  ++numphsi;
	  atomdihedpairs[0]=(int *)gcerealloc(atomdihedpairs[0],sizeof(int)*numphsi*4);
	  atomdihedpairs[0][numphsi*4-4]=psi[0];
	  atomdihedpairs[0][numphsi*4-3]=psi[1];
	  atomdihedpairs[0][numphsi*4-2]=psi[2];
	  atomdihedpairs[0][numphsi*4-1]=i;
	}
	PSIFLAG=0;
	psi[0]=i;
      }
    }
    numres=PTL_joinatomtores(i,RESNAME);
    if (strncmp(RESNAME,"ALA",3)==0) {
      if (strncmp(ap.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
      }
      if (strncmp(ap.IGRAPH[i],"HB1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"ACE",3)==0) {
      ACEFLAG=1;
      if (strncmp(ap.IGRAPH[i],"HH31",4)==0){
    	ipsi[0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CH3",3)==0){
    	ipsi[1]=i;
	iomega[0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"C",1)==0){
    	ipsi[2]=i;
	iomega[1]=i;
      }
    }
    else if (strncmp(RESNAME,"NME",3)==0) {
      NMEFLAG=1;
      if (strncmp(ap.IGRAPH[i],"N",1)==0){
    	fomega[0]=i;
	fphi[0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CH3",3)==0){
	atomdihedpairs[4][0]=phi[2];
	atomdihedpairs[4][1]=phi[3];
	atomdihedpairs[4][2]=fomega[0];
	atomdihedpairs[4][3]=i;
	fphi[1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"HH31",4)==0){
	atomdihedpairs[4][4]=phi[3];
	atomdihedpairs[4][5]=fphi[0];
	atomdihedpairs[4][6]=fphi[1];
	atomdihedpairs[4][7]=i;
      }
    }
    /*   else if (strncmp(ap.IGRAPH[i],"N",1)==0){	    */
    /* 	;						    */
    /*   }						    */
    /*   else if (strncmp(ap.IGRAPH[i],"CB",1)==0){	    */
    /* 	;						    */
    /*   }						    */
    /*   else if (strncmp(ap.IGRAPH[i],"HA1",3)==0){	    */
    /* 	;						    */
    /*   }						    */
    /* }						    */

    else if (strncmp(RESNAME,"ASP",3)==0) {
      if (strncmp(ap.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
	kai[2][1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"OD1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"GLU",3)==0) {
      if (strncmp(ap.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
	kai[2][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
	kai[2][1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CD",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[2][2]=i;
      }
      if (strncmp(ap.IGRAPH[i],"OE1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[2][0];
	atomdihedpairs[2][numkai*4-3]=kai[2][1];
	atomdihedpairs[2][numkai*4-2]=kai[2][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"LEU",3)==0) {
      if (strncmp(ap.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
	kai[2][0]=i;
	kai[3][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
	kai[2][1]=i;
	kai[3][1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CD1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[2][2]=i;
      }
      if (strncmp(ap.IGRAPH[i],"HD1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[2][0];
	atomdihedpairs[2][numkai*4-3]=kai[2][1];
	atomdihedpairs[2][numkai*4-2]=kai[2][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CD2",3)==0){
	kai[3][2]=i;
      }
      if (strncmp(ap.IGRAPH[i],"HD2",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[3][0];
	atomdihedpairs[2][numkai*4-3]=kai[3][1];
	atomdihedpairs[2][numkai*4-2]=kai[3][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"ILE",3)==0) {
      if (strncmp(ap.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
	kai[2][0]=i;
	kai[3][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CG2",3)==0){
	kai[1][2]=i;
      }
      if (strncmp(ap.IGRAPH[i],"HG21",4)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"GG1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[2][2]=i;
	kai[3][1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CD1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[2][0];
	atomdihedpairs[2][numkai*4-3]=kai[2][1];
	atomdihedpairs[2][numkai*4-2]=kai[2][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[3][2]=i;
      }
      if (strncmp(ap.IGRAPH[i],"HD1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[3][0];
	atomdihedpairs[2][numkai*4-3]=kai[3][1];
	atomdihedpairs[2][numkai*4-2]=kai[3][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }

    else if (strncmp(RESNAME,"ASN",3)==0) {
      if (strncmp(ap.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
	kai[2][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
	kai[2][1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"ND2",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[2][2]=i;
      }
      if (strncmp(ap.IGRAPH[i],"HD21",4)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[2][0];
	atomdihedpairs[2][numkai*4-3]=kai[2][1];
	atomdihedpairs[2][numkai*4-2]=kai[2][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }

    else if (strncmp(RESNAME,"GLN",3)==0) {
      if (strncmp(ap.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
	kai[2][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
	kai[2][1]=i;
	kai[3][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CD",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[2][2]=i;
	//	kai[2][1]=i;
	kai[3][1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"NE2",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[2][0];
	atomdihedpairs[2][numkai*4-3]=kai[2][1];
	atomdihedpairs[2][numkai*4-2]=kai[2][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[3][2]=i;
      }
      if (strncmp(ap.IGRAPH[i],"HE21",4)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[3][0];
	atomdihedpairs[2][numkai*4-3]=kai[3][1];
	atomdihedpairs[2][numkai*4-2]=kai[3][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }

    else if (strncmp(RESNAME,"VAL",3)==0) {
      if (strncmp(ap.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
	kai[2][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CG1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
      }
      if (strncmp(ap.IGRAPH[i],"HG1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CG2",3)==0){
	kai[2][2]=i;
      }
      if (strncmp(ap.IGRAPH[i],"HG21",4)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[2][0];
	atomdihedpairs[2][numkai*4-3]=kai[2][1];
	atomdihedpairs[2][numkai*4-2]=kai[2][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"SER",3)==0) {
      if (strncmp(ap.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"OG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
      }
      if (strncmp(ap.IGRAPH[i],"HG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"THR",3)==0) {
      if (strncmp(ap.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
	kai[2][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
	kai[2][1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CG2",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
      }
      if (strncmp(ap.IGRAPH[i],"HG21",4)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"OG1",3)==0){
	kai[2][2]=i;
      }
      if (strncmp(ap.IGRAPH[i],"HG1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[2][0];
	atomdihedpairs[2][numkai*4-3]=kai[2][1];
	atomdihedpairs[2][numkai*4-2]=kai[2][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"CYX",3)==0) {
      if (strncmp(ap.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
      }
      if (strncmp(ap.IGRAPH[i],"SG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"CYS",3)==0) {
      if (strncmp(ap.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"SG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
      }
      if (strncmp(ap.IGRAPH[i],"HG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"HID",3)==0) {
      ;
    }
    else if (strncmp(RESNAME,"HIE",3)==0) {
      if (strncmp(ap.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
      }
      if (strncmp(ap.IGRAPH[i],"ND1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"HIP",3)==0) {
      if (strncmp(ap.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
      }
      if (strncmp(ap.IGRAPH[i],"ND1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"MET",3)==0) {
      if (strncmp(ap.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
	kai[2][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
	kai[2][1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"SD",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[2][2]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CE",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[2][0];
	atomdihedpairs[2][numkai*4-3]=kai[2][1];
	atomdihedpairs[2][numkai*4-2]=kai[2][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    /*******************************************/
    /* else if (strncmp(RESNAME,"PRO",3)==0) { */
    /*   ;				       */
    /* }				       */
    /*******************************************/
    else if (strncmp(RESNAME,"ARG",3)==0) {
      if (strncmp(ap.IGRAPH[i],"N\0",2)==0){
    	kai[0][0]=i;
      }
      else if (strncmp(ap.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      else if (strncmp(ap.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
	kai[2][0]=i;
      }
      else if (strncmp(ap.IGRAPH[i],"CG",2)==0){
	++numkai;
	kai[0][3]=i;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
	kai[2][1]=i;
	kai[3][0]=i;
      }
      else if (strncmp(ap.IGRAPH[i],"CD",2)==0){
	++numkai;
	kai[1][3]=i;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[2][2]=i;
	kai[3][1]=i;
	kai[4][0]=i;
      }
      else if (strncmp(ap.IGRAPH[i],"NE",2)==0){
	++numkai;
	kai[2][3]=i;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[2][0];
	atomdihedpairs[2][numkai*4-3]=kai[2][1];
	atomdihedpairs[2][numkai*4-2]=kai[2][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[3][2]=i;
	kai[4][1]=i;
	kai[5][0]=i;
	kai[6][0]=i;
      }
      else if (strncmp(ap.IGRAPH[i],"CZ",2)==0){
	++numkai;
	kai[3][3]=i;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[3][0];
	atomdihedpairs[2][numkai*4-3]=kai[3][1];
	atomdihedpairs[2][numkai*4-2]=kai[3][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[4][2]=i;
	kai[5][1]=i;
	kai[6][1]=i;
      }
      else if (strncmp(ap.IGRAPH[i],"NH1",3)==0){
	++numkai;
	kai[4][3]=i;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[4][0];
	atomdihedpairs[2][numkai*4-3]=kai[4][1];
	atomdihedpairs[2][numkai*4-2]=kai[4][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[5][2]=i;
	kai[6][2]=i;
      }
      else if (strncmp(ap.IGRAPH[i],"HH11",4)==0){
	++numkai;
	kai[5][3]=i;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[5][0];
	atomdihedpairs[2][numkai*4-3]=kai[5][1];
	atomdihedpairs[2][numkai*4-2]=kai[5][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
      else if (strncmp(ap.IGRAPH[i],"HH21",4)==0){
	++numkai;
	kai[6][3]=i;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[6][0];
	atomdihedpairs[2][numkai*4-3]=kai[6][1];
	atomdihedpairs[2][numkai*4-2]=kai[6][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"LYS",3)==0) {
      if (strncmp(ap.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
	kai[2][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
	kai[2][1]=i;
	kai[3][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CD",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[2][2]=i;
	kai[3][1]=i;
	kai[4][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CE",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[2][0];
	atomdihedpairs[2][numkai*4-3]=kai[2][1];
	atomdihedpairs[2][numkai*4-2]=kai[2][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[4][1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"NZ",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[3][0];
	atomdihedpairs[2][numkai*4-3]=kai[3][1];
	atomdihedpairs[2][numkai*4-2]=kai[3][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[4][2]=i;
      }
      if (strncmp(ap.IGRAPH[i],"HZ1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[4][0];
	atomdihedpairs[2][numkai*4-3]=kai[4][1];
	atomdihedpairs[2][numkai*4-2]=kai[4][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"PHE",3)==0) {
      if (strncmp(ap.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CD1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"TYR",3)==0) {
      if (strncmp(ap.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CD1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CE1",3)==0){
    	kai[2][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CZ",2)==0){
	kai[2][1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"OH",2)==0){
	kai[2][2]=i;
      }
      if (strncmp(ap.IGRAPH[i],"HH",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[2][0];
	atomdihedpairs[2][numkai*4-3]=kai[2][1];
	atomdihedpairs[2][numkai*4-2]=kai[2][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"TRP",3)==0) {
      if (strncmp(ap.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
      }
      if (strncmp(ap.IGRAPH[i],"CD1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }

    if (numres==1 && ACEFLAG==1) {
      if (strncmp(ap.IGRAPH[i],"N",1)==0){
	iomega[2]=i;
	atomdihedpairs[3][0]=ipsi[0];
	atomdihedpairs[3][1]=ipsi[1];
	atomdihedpairs[3][2]=ipsi[2];
	atomdihedpairs[3][3]=i;
      }
      else if (strncmp(ap.IGRAPH[i],"C",1)==0){
	iomega[3]=i;
	atomdihedpairs[3][4]=iomega[0];
	atomdihedpairs[3][5]=iomega[1];
	atomdihedpairs[3][6]=iomega[2];
	atomdihedpairs[3][7]=i;
      }
    }
  }

  num[0]=numphsi;num[1]=numomega;num[2]=numkai;
  if (ACEFLAG==1)
    num[3]=2;
  if (NMEFLAG==1)
    num[4]=2;

  sumnum=0;
  for (i=0;i<4;++i)
    sumnum+=num[i];
  
  return sumnum;
}

int PTL_resb_ca(int numres,struct AmberParmL ap) {
  int i;

  for (i=ap.IPRES[numres]-1;i<ap.IPRES[numres+1]-1;++i) {
    if (strncmp(ap.IGRAPH[i],"CA",2)==0) {
      return i;
      break;
    }
  }

  if (i==ap.IPRES[numres+1]-1)
    return -1;
}

/**********************************************************************/
/* int *PTL_ca_res_set(void) {					      */
/*   int i,ii;							      */
/*   int numatom;						      */
/*   int *index_ca_res;						      */
/* 								      */
/*   numatom=ap.NATOM;						      */
/*   index_ca_res=(int *)gcemalloc(sizeof(int)*1);		      */
/* 								      */
/*   ii=0;							      */
/*   for (i=0;i<numatom;++i) {					      */
/*     if (strncmp(ap.IGRAPH[i],"CA",2)==0) {			      */
/*       index_ca_res[ii]=i;					      */
/*       ++ii;							      */
/*       index_ca_res=(int *)gcerealloc(index_ca_res,sizeof(int)*ii); */
/*     }							      */
/*   }								      */
/* 								      */
/*   return index_ca_res;					      */
/* }								      */
/* 								      */
/* int PTL_ca_res(int* index_ca_res,int numatom) {		      */
/*   int i;							      */
/* 								      */
/*   for (i=ap.IPRES[numres]-1;i<ap.IPRES[numres+1]-1;++i) {	      */
/*     if (strncmp(ap.IGRAPH[i],"CA",1)==0) {			      */
/*       return i;						      */
/*       break;							      */
/*     }							      */
/*   }								      */
/* 								      */
/*   if (i==ap.IPRES[numres+1]-1)				      */
/*     return -1;						      */
/* }								      */
/**********************************************************************/

int PTL_joinatomtoresb(int numatom, char LABERES[4],struct AmberParmL ap) {
  int i,j;

  for (i=0;i<ap.NRES;++i) {
    if (numatom >= ap.IPRES[i]-1 ) {
      if (i==(ap.NRES-1)){
	for (j=0;j<3;++j) {
	  LABERES[j]=ap.LABERES[i][j];
	}
	LABERES[3]='\0';
      }
      else {
	if (numatom < ap.IPRES[i+1]-1) {
	  for (j=0;j<3;++j) {
	    LABERES[j]=ap.LABERES[i][j];
	  }
	  LABERES[3]='\0';
	  break;
	}
      }
    }
  }
  
  return i;
}

int PTL_resnumb(int numatom,struct AmberParmL ap) {
  int i;

  if (numatom>ap.NATOM) return -1;

  for (i=0;i<ap.NRES;++i) {
    if (ap.IPRES[i]-1 > numatom)
      return i;
  }
  return ap.NRES;
}

int PTL_resnum2b(int numatom,struct AmberParmL ap) {
  int i;

  if (numatom>ap.NATOM) return -1;

  for (i=0;i<ap.NRES;++i) {
    //    if (ap.IPRES[i]-1 > numatom)
    if (ap.IPRES[i]-1 > numatom) //12-02-20
      return i;
  }
  return ap.NRES;
}

int PTL_canum_fromresnumb(int numres,struct AmberParmL ap) {
  int i;
  int fin;

  if (numres==ap.NRES-1) 
    fin=ap.NATOM;
  else 
    fin=ap.IPRES[numres+1];

  for (i=ap.IPRES[numres]-1;i<fin;++i) {
    if (strncmp(ap.IGRAPH[i],"CA",2)==0) {
      return i;
      break;
    }
  }

  return -1;
}

int PTL_which_includeb(int numres,int *listres,int numlistres,struct AmberParmL ap){
  int i,j,k;
  int flag=1;

  for (i=0;i<numlistres;++i) {
    if (numres==listres[i]) {
      flag=0;
      break;
    }
  }

  return flag;
}
