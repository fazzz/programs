#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>

#include "PTL.h"
#include "EF.h"

int readParmtopL(FILE *parmfile){
  int i,j;
  char *line;
  size_t len=0;

  double *a; // 2014-07-15
  //  struct AmberParmL APP; // 2014-07-15
  getline(&line,&len,parmfile);
  getline(&line,&len,parmfile);
  getline(&line,&len,parmfile);
  getline(&line,&len,parmfile);
  //  fscanf(parmfile,"%4s",&AP.ITITLE);

  for (i=0;i<2;++i){
    getline(&line,&len,parmfile);
  }
  fscanf(parmfile,"%8d",&AP.NATOM);
  fscanf(parmfile,"%8d",&AP.NTYPES);
  fscanf(parmfile,"%8d",&AP.NBONH);
  fscanf(parmfile,"%8d",&AP.MBONA);
  fscanf(parmfile,"%8d",&AP.NTHETH);
  fscanf(parmfile,"%8d",&AP.MTHETA);
  fscanf(parmfile,"%8d",&AP.NPHIH);
  fscanf(parmfile,"%8d",&AP.MPHIA);
  fscanf(parmfile,"%8d",&AP.NHPARM);
  fscanf(parmfile,"%8d",&AP.NPARM);
  fscanf(parmfile,"%8d",&AP.NNB);
  AP.NEXT = AP.NNB;
  fscanf(parmfile,"%8d",&AP.NRES);
  fscanf(parmfile,"%8d",&AP.NBONA);
  fscanf(parmfile,"%8d",&AP.NTHETA);
  fscanf(parmfile,"%8d",&AP.NPHIA);
  fscanf(parmfile,"%8d",&AP.NUMBND);
  fscanf(parmfile,"%8d",&AP.NUMANG);
  fscanf(parmfile,"%8d",&AP.NPTRA);
  fscanf(parmfile,"%8d",&AP.NATYP);
  fscanf(parmfile,"%8d",&AP.NPHB);
  fscanf(parmfile,"%8d",&AP.IFPERT);
  fscanf(parmfile,"%8d",&AP.NBPER);
  fscanf(parmfile,"%8d",&AP.NGPER);
  fscanf(parmfile,"%8d",&AP.NDPER);
  fscanf(parmfile,"%8d",&AP.MBPER);
  fscanf(parmfile,"%8d",&AP.MGPER);
  fscanf(parmfile,"%8d",&AP.MDPER);
  fscanf(parmfile,"%8d",&AP.IFBOX);
  fscanf(parmfile,"%8d",&AP.NMXPS);
  fscanf(parmfile,"%8d",&AP.IFCAP);
  fscanf(parmfile,"%8d",&AP.NEXTRA);

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NATOM;++i){
    // fscanf(parmfile,"%4s",&AP.IGRAPH[i]); // 2014-07-04
    for (j=0;j<4;++j) {
      fscanf(parmfile,"%c",&AP.IGRAPH[i][j]);
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  
  a=(double *)malloc(sizeof(double)*2); // 2014-07-17
  free(a);                              // 2014-07-17
  //  a=(double *)GC_malloc(sizeof(double)*2); // 2014-07-15 
  //  a=(double *)gcemalloc(sizeof(double)*2); // 2014-07-15
  //AP_test.CHRG=(double *)gcemalloc(sizeof(double)*10); // 2014-07-15
  //  a=(double *)gcemalloc(sizeof(double)*AP.NATOM); // 2014-07-17
  //  AP.CHRG=(double *)gcemalloc(sizeof(double)*AP.NATOM);
  AP.CHRG=(double *)malloc(sizeof(double)*AP.NATOM); // 2014-07-17
  //  AP.CHRG=(double *)GC_malloc(sizeof(double)*AP.NATOM); // 2014-07-15
  for (i=0;i<AP.NATOM;++i){
    fscanf(parmfile,"%lf",&AP.CHRG[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.AMASS=(double *)gcemalloc(sizeof(double)*AP.NATOM);
  for (i=0;i<AP.NATOM;++i){
    fscanf(parmfile,"%lf",&AP.AMASS[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.IAC=(int/*double*/ *)gcemalloc(sizeof(int/*double*/)*AP.NATOM); // 0811
  for (i=0;i<AP.NATOM;++i){
    fscanf(parmfile,"%8d",&AP.IAC[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.NUMEX=(int/*double*/ *)gcemalloc(sizeof(int/*double*/)*AP.NATOM); // 0811
  for (i=0;i<AP.NATOM;++i){
    fscanf(parmfile,"%8d",&AP.NUMEX[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.ICO=(int/*double*/ *)gcemalloc(sizeof(int/*double*/)*(AP.NTYPES)*(AP.NTYPES)); // 0811
  for (i=0;i<(AP.NTYPES)*(AP.NTYPES);++i){
    fscanf(parmfile,"%8d",&AP.ICO[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NRES;++i){
    //    fscanf(parmfile,"%4s",&AP.LABERES[i]);
    for (j=0;j<4;++j) {
      fscanf(parmfile,"%c",&AP.LABERES[i][j]);
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.IPRES=(int *)gcemalloc(sizeof(int)*AP.NRES);
  for (i=0;i<AP.NRES;++i){
    fscanf(parmfile,"%8d",&AP.IPRES[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.RK=(double *)gcemalloc(sizeof(double)*AP.NUMBND);
  for (i=0;i<AP.NUMBND;++i){
    fscanf(parmfile,"%lf",&AP.RK[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.REQ=(double *)gcemalloc(sizeof(double)*AP.NUMBND);
  for (i=0;i<AP.NUMBND;++i){
    fscanf(parmfile,"%lf",&AP.REQ[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.TK=(double *)gcemalloc(sizeof(double)*AP.NUMANG);
  for (i=0;i<AP.NUMANG;++i){
    fscanf(parmfile,"%lf",&AP.TK[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.TEQ=(double *)gcemalloc(sizeof(double)*AP.NUMANG);
  for (i=0;i<AP.NUMANG;++i){
    fscanf(parmfile,"%lf",&AP.TEQ[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.PK=(double *)gcemalloc(sizeof(double)*AP.NPTRA);
  for (i=0;i<AP.NPTRA;++i){
    fscanf(parmfile,"%lf",&AP.PK[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.PN=(double *)gcemalloc(sizeof(double)*AP.NPTRA);
  for (i=0;i<AP.NPTRA;++i){
    fscanf(parmfile,"%lf",&AP.PN[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.PHASE=(double *)gcemalloc(sizeof(double)*AP.NPTRA);
  for (i=0;i<AP.NPTRA;++i){
    fscanf(parmfile,"%lf",&AP.PHASE[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.SOLTY=(double *)gcemalloc(sizeof(double)*AP.NATYP);
  for (i=0;i<AP.NATYP;++i){
    fscanf(parmfile,"%lf",&AP.SOLTY[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.CN1=(double *)gcemalloc(sizeof(double)*(AP.NTYPES)*(AP.NTYPES+1)/2);
  for (i=0;i<(AP.NTYPES)*(AP.NTYPES+1)/2;++i){
    fscanf(parmfile,"%lf",&AP.CN1[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.CN2=(double *)gcemalloc(sizeof(double)*(AP.NTYPES)*(AP.NTYPES+1)/2);
  for (i=0;i<(AP.NTYPES)*(AP.NTYPES+1)/2;++i){
    fscanf(parmfile,"%lf",&AP.CN2[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.BH=(int/*double*/ **)gcemalloc(sizeof(int/*double*/ *)*AP.NBONH); // 0811
  for (i=0;i<AP.NBONH;++i) AP.BH[i]=(int/*double*/ *)gcemalloc(sizeof(int/*double*/)*3);
  for (i=0;i<AP.NBONH;++i){
    for (j=0;j<3;++j){
      fscanf(parmfile,"%8d",&AP.BH[i][j]);
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.BA=(int/*double*/ **)gcemalloc(sizeof(int/*double*/ *)*AP.NBONA);
  for (i=0;i<AP.NBONA;++i) AP.BA[i]=(int/*double*/ *)gcemalloc(sizeof(int/*double*/)*3);
  for (i=0;i<AP.NBONA;++i){
    for (j=0;j<3;++j){
      fscanf(parmfile,"%d",&AP.BA[i][j]);
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.TH=(int/*double*/ **)gcemalloc(sizeof(int/*double*/ *)*AP.NTHETH);
  for (i=0;i<AP.NTHETH;++i) AP.TH[i]=(int/*double*/ *)gcemalloc(sizeof(int/*double*/)*4);
  for (i=0;i<AP.NTHETH;++i){
    for (j=0;j<4;++j){
      fscanf(parmfile,"%8d",&AP.TH[i][j]);
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.TA=(int/*double*/ **)gcemalloc(sizeof(int/*double*/ *)*AP.NTHETA);
  for (i=0;i<AP.NTHETA;++i) AP.TA[i]=(int/*double*/ *)gcemalloc(sizeof(int/*double*/)*4);
  for (i=0;i<AP.NTHETA;++i){
    for (j=0;j<4;++j){
      fscanf(parmfile,"%d",&AP.TA[i][j]);
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.PH=(int/*double*/ **)gcemalloc(sizeof(int/*double*/ *)*AP.NPHIH);
  for (i=0;i<AP.NPHIH;++i) AP.PH[i]=(int/*double*/ *)gcemalloc(sizeof(int/*double*/)*5);
  for (i=0;i<AP.NPHIH;++i){
    for (j=0;j<5;++j){
      fscanf(parmfile,"%8d",&AP.PH[i][j]);
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.PA=(int/*double*/ **)gcemalloc(sizeof(int/*double*/ *)*AP.NPHIA);
  for (i=0;i<AP.NPHIA;++i) AP.PA[i]=(int/*double*/ *)gcemalloc(sizeof(int/*double*/)*5);
  for (i=0;i<AP.NPHIA;++i){
    for (j=0;j<5;++j){
      fscanf(parmfile,"%8d",&AP.PA[i][j]);
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.NATEX=(int *)gcemalloc(sizeof(int)*AP.NEXT);
  for (i=0;i<AP.NEXT;++i){
    fscanf(parmfile,"%8d",&AP.NATEX[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.ASOL=(double *)gcemalloc(sizeof(double)*AP.NPHB);
  for (i=0;i<AP.NPHB;++i){
    fscanf(parmfile,"%lf",&AP.ASOL[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.BSOL=(double *)gcemalloc(sizeof(double)*AP.NPHB);
  for (i=0;i<AP.NPHB;++i){
    fscanf(parmfile,"%lf",&AP.BSOL[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.HBCUT=(double *)gcemalloc(sizeof(double)*AP.NPHB);
  for (i=0;i<AP.NPHB;++i){
    fscanf(parmfile,"%lf",&AP.HBCUT[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NATOM;++i){
    //    fscanf(parmfile,"%4s",&AP.ISYMBL[i]); // 2014-07-04
    for (j=0;j<4;++j) {                         // 2014-07-04
      fscanf(parmfile,"%c",&AP.ISYMBL[i][j]);   // 2014-07-04
    }                                           // 2014-07-04
  }
  
  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NATOM;++i){
    //    fscanf(parmfile,"%4s",&AP.ITREE[i]); // 2014-07-04
    for (j=0;j<4;++j) {                        // 2014-07-04
      fscanf(parmfile,"%c",&AP.ITREE[i][j]);   // 2014-07-04
    }                                          // 2014-07-04
  }
  
  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.JOIN=(int *)gcemalloc(sizeof(int)*AP.NATOM);
  for (i=0;i<AP.NATOM;++i){
    fscanf(parmfile,"%8d",&AP.JOIN[i]);
  }
  
  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.IROTAT=(int *)gcemalloc(sizeof(int)*AP.NATOM);
  for (i=0;i<AP.NATOM;++i){
    fscanf(parmfile,"%8d",&AP.IROTAT[i]);
  }
  
  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  if (AP.IFBOX > 0) {
    fscanf(parmfile,"%d",&AP.IPTRES);
    
    fscanf(parmfile,"%d",&AP.NSPM);
    
    fscanf(parmfile,"%d",&AP.NSPSOL);
    
    AP.NSP=(int *)gcemalloc(sizeof(int)*AP.NSPM);
    for (i=0;i<AP.NSPM;++i){
      fscanf(parmfile,"%d",&AP.NSP[i]);
    }
    
    fscanf(parmfile,"%lf",&AP.BETA);
    
    for (i=0;i<3;++i){
      //      fscanf(parmfile,"%12d",&AP.BOX[i]); // 2014-07-04
      fscanf(parmfile,"%lf",&AP.BOX[i]); // 2014-07-04
    }
  }
  
  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  if (AP.IFCAP > 0) {
    fscanf(parmfile,"%d",&AP.NATCAP);
    
    fscanf(parmfile,"%lf",&AP.CUTCAP);
    
    fscanf(parmfile,"%lf",&AP.XCAP);
    
    fscanf(parmfile,"%lf",&AP.YCAP);
    
    fscanf(parmfile,"%lf",&AP.ZCAP);
  }
  
  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  if (AP.IFPERT > 0) {
    AP.BPER=(int **)gcemalloc(sizeof(int *)*AP.NBPER);
    for (i=0;i<AP.NBPER;++i) AP.BPER[i]=(int *)gcemalloc(sizeof(int)*2);
    for (i=0;i<AP.NBPER;++i){
      for (j=0;j<2;++j){
	fscanf(parmfile,"%d",&AP.BPER[i][j]);
      }
    }
    
    AP.ICBPER=(int *)gcemalloc(sizeof(int)*AP.NBPER*1);
    for (i=0;i<AP.NBPER*2;++i){
      fscanf(parmfile,"%d",&AP.ICBPER[i]);
    }
    
    AP.TPER=(int **)gcemalloc(sizeof(int *)*AP.NGPER);
    for (i=0;i<AP.NGPER;++i) AP.TPER[i]=(int *)gcemalloc(sizeof(int)*3);
    for (i=0;i<AP.NGPER;++i){
      for (j=0;j<3;++j){
	fscanf(parmfile,"%d",&AP.TPER[i][j]);
      }
    }
    
    AP.ICTPER=(int *)gcemalloc(sizeof(int)*AP.NGPER*2);
    for (i=0;i<AP.NGPER*2;++i){
      fscanf(parmfile,"%d",&AP.ICTPER[i]);
    }

    AP.PPER=(int **)gcemalloc(sizeof(int *)*AP.NDPER);
    //    for (i=0;i<AP.NDPER;++i)   AP.PPER=(int *)gcemalloc(sizeof(int)*4);
    for (i=0;i<AP.NDPER;++i)   AP.PPER=(int **)gcemalloc(sizeof(int *)*4); // 2014-06-17
    for (i=0;i<AP.NDPER;++i) {
      for (j=0;j<4;++j){
	fscanf(parmfile,"%d",&AP.PPER[i][j]);
      }
    }

    AP.ICTPER=(int *)gcemalloc(sizeof(int)*AP.NDPER*2);
    for (i=0;i<AP.NDPER*2;++i){
      fscanf(parmfile,"%d",&AP.ICPPER[i]);
    }

    for (i=0;i<AP.NRES;++i){
      //      fscanf(parmfile,"%4s",&AP.LABERES[i]); // 2014-07-04
      for (j=0;j<4;++j) {                            // 2014-07-04
	fscanf(parmfile,"%c",&AP.LABERES[i][j]);     // 2014-07-04
      }                                              // 2014-07-04
    }

    for (i=0;i<AP.NATOM;++i){
      //      fscanf(parmfile,"%4s",&AP.IGRPER[i]);  // 2014-07-04
      for (j=0;j<4;++j) {                            // 2014-07-04
	fscanf(parmfile,"%c",&AP.IGRPER[i][j]);      // 2014-07-04
      }                                              // 2014-07-04
    }

    for (i=0;i<AP.NATOM;++i){
      //      fscanf(parmfile,"%4s",&AP.ISMPER[i]);  // 2014-07-04
      for (j=0;j<4;++j) {                            // 2014-07-04
	fscanf(parmfile,"%c",&AP.ISMPER[i][j]);      // 2014-07-04
      }                                              // 2014-07-04
    }

    for (i=0;i<AP.NATOM;++i){
      fscanf(parmfile,"%lf",&AP.ALMPER[i]);
    }

    //    AP.IAPER=(int *)gcemalloc(sizeof(int)*AP.NATOM);
    AP.IAPER=(double *)gcemalloc(sizeof(double)*AP.NATOM); // 2014-06-17
    for (i=0;i<AP.NATOM;++i){
      fscanf(parmfile,"%lf",&AP.IAPER[i]);
    }

    AP.IACPER=(int *)gcemalloc(sizeof(int)*AP.NATOM);
    for (i=0;i<AP.NATOM;++i){
      fscanf(parmfile,"%d",&AP.IACPER[i]);
    }

    //    AP.CGPER=(int *)gcemalloc(sizeof(int)*AP.NATOM);
    AP.CGPER=(double *)gcemalloc(sizeof(double)*AP.NATOM); // 2014-06-17
    for (i=0;i<AP.NATOM;++i){
      fscanf(parmfile,"%lf",&AP.CGPER[i]);
    }
} 

  /* if (AP.IPOL == 1) {
    fscanf(parmfile,"%18.8e",&AP.ATPOL[i]);

    fscanf(parmfile,"%18.8e",&AP.ATPOL1[i]);
    }*/

if (AP.NPARM == 1) {
    fscanf(parmfile,"%d",&AP.NLES_NTYP);

    AP.LES_TYPE=(int/*double*/ *)gcemalloc(sizeof(int/*double*/)*AP.NATOM);
    for (i=0;i<AP.NATOM;++i){
      //      fscanf(parmfile,"%lf",&AP.LES_TYPE[i]);
      fscanf(parmfile,"%d",&AP.LES_TYPE[i]);
    }

    AP.LES_FAC=(double *)gcemalloc(sizeof(double)*AP.NATOM);
    for (i=0;i<(AP.NLES_NTYP)*(AP.NLES_NTYP);++i){
      fscanf(parmfile,"%lf",&AP.LES_FAC[i]);
    }

    AP.LES_CNUM=(double *)gcemalloc(sizeof(double)*AP.NATOM);
    for (i=0;i<AP.NATOM;++i){
      fscanf(parmfile,"%lf",&AP.LES_CNUM[i]);
    }

    AP.LES_ID=(double *)gcemalloc(sizeof(double)*AP.NATOM);
    for (i=0;i<AP.NATOM;++i){
      fscanf(parmfile,"%lf",&AP.LES_ID[i]);
    }

  }

 return 0; // 2014-07-04
}

int readdihedpairsL(int **atomdihedpairs, int *num) {
  int i,j;
  int phi[4],psi[4],omega[4],ipsi[4],iomega[4],fomega[4],fphi[4],sumnum;
  int kai[10][4];
  int PHIFLAG=-1,PSIFLAG=-1,OMEGAFLAG=-1,ACEFLAG,NMEFLAG/*2014-03-23*/,ACEFLAG2=1/*2014-03-23*/;
  int numphsi,numomega,numkai,numres;
  char RESNAME[4];

  numphsi=0;numomega=0;numkai=0;
  
  for (i=0;i<AP.NATOM;++i) {
    if (strncmp(AP.ITREE[i],"M",1)==0) {
      if (strncmp(AP.IGRAPH[i],"CH3",3)==0){
	;
      }
      else if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
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
      else if (strncmp(AP.IGRAPH[i],"C",1)==0){
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
      else if (strncmp(AP.IGRAPH[i],"N",1)==0){
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
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"HB1",3)==0){
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
      if (strncmp(AP.IGRAPH[i],"HH31",4)==0){
    	ipsi[0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CH3",3)==0){
    	ipsi[1]=i;
	iomega[0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"C",1)==0){
    	ipsi[2]=i;
	iomega[1]=i;
      }
    }
    else if (strncmp(RESNAME,"NME",3)==0) {
      NMEFLAG=1;
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
    	fomega[0]=i;
	fphi[0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CH3",3)==0){
	atomdihedpairs[4][0]=phi[2];
	atomdihedpairs[4][1]=phi[3];
	atomdihedpairs[4][2]=fomega[0];
	atomdihedpairs[4][3]=i;
	fphi[1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"HH31",4)==0){
	atomdihedpairs[4][4]=phi[3];
	atomdihedpairs[4][5]=fphi[0];
	atomdihedpairs[4][6]=fphi[1];
	atomdihedpairs[4][7]=i;
      }
    }
    /*   else if (strncmp(AP.IGRAPH[i],"N",1)==0){	    */
    /* 	;						    */
    /*   }						    */
    /*   else if (strncmp(AP.IGRAPH[i],"CB",1)==0){	    */
    /* 	;						    */
    /*   }						    */
    /*   else if (strncmp(AP.IGRAPH[i],"HA1",3)==0){	    */
    /* 	;						    */
    /*   }						    */
    /* }						    */

    else if (strncmp(RESNAME,"ASP",3)==0) {
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
	kai[2][1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"OD1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"GLU",3)==0) {
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
	kai[2][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
	kai[2][1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CD",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[2][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"OE1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[2][0];
	atomdihedpairs[2][numkai*4-3]=kai[2][1];
	atomdihedpairs[2][numkai*4-2]=kai[2][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"LEU",3)==0) {
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
	kai[2][0]=i;
	kai[3][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CG",2)==0){
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
      if (strncmp(AP.IGRAPH[i],"CD1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[2][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"HD1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[2][0];
	atomdihedpairs[2][numkai*4-3]=kai[2][1];
	atomdihedpairs[2][numkai*4-2]=kai[2][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CD2",3)==0){
	kai[3][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"HD2",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[3][0];
	atomdihedpairs[2][numkai*4-3]=kai[3][1];
	atomdihedpairs[2][numkai*4-2]=kai[3][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"ILE",3)==0) {
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
	kai[2][0]=i;
	kai[3][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CG2",3)==0){
	kai[1][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"HG21",4)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"GG1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[2][2]=i;
	kai[3][1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CD1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[2][0];
	atomdihedpairs[2][numkai*4-3]=kai[2][1];
	atomdihedpairs[2][numkai*4-2]=kai[2][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[3][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"HD1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[3][0];
	atomdihedpairs[2][numkai*4-3]=kai[3][1];
	atomdihedpairs[2][numkai*4-2]=kai[3][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }

    else if (strncmp(RESNAME,"ASN",3)==0) {
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
	kai[2][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
	kai[2][1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"ND2",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[2][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"HD21",4)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[2][0];
	atomdihedpairs[2][numkai*4-3]=kai[2][1];
	atomdihedpairs[2][numkai*4-2]=kai[2][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }

    else if (strncmp(RESNAME,"GLN",3)==0) {
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
	kai[2][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CG",2)==0){
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
      if (strncmp(AP.IGRAPH[i],"CD",2)==0){
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
      if (strncmp(AP.IGRAPH[i],"NE2",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[2][0];
	atomdihedpairs[2][numkai*4-3]=kai[2][1];
	atomdihedpairs[2][numkai*4-2]=kai[2][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[3][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"HE21",4)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[3][0];
	atomdihedpairs[2][numkai*4-3]=kai[3][1];
	atomdihedpairs[2][numkai*4-2]=kai[3][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }

    else if (strncmp(RESNAME,"VAL",3)==0) {
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
	kai[2][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CG1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"HG1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CG2",3)==0){
	kai[2][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"HG21",4)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[2][0];
	atomdihedpairs[2][numkai*4-3]=kai[2][1];
	atomdihedpairs[2][numkai*4-2]=kai[2][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"SER",3)==0) {
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"OG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"HG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"THR",3)==0) {
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
	kai[2][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
	kai[2][1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CG2",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"HG21",4)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"OG1",3)==0){
	kai[2][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"HG1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[2][0];
	atomdihedpairs[2][numkai*4-3]=kai[2][1];
	atomdihedpairs[2][numkai*4-2]=kai[2][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"CYX",3)==0) {
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"SG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"CYS",3)==0) {
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"SG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"HG",2)==0){
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
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"ND1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"HIP",3)==0) {
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"ND1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"MET",3)==0) {
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
	kai[2][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
	kai[2][1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"SD",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[2][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CE",2)==0){
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
      if (strncmp(AP.IGRAPH[i],"N\0",2)==0){
    	kai[0][0]=i;
      }
      else if (strncmp(AP.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      else if (strncmp(AP.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
	kai[2][0]=i;
      }
      else if (strncmp(AP.IGRAPH[i],"CG",2)==0){
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
      else if (strncmp(AP.IGRAPH[i],"CD",2)==0){
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
      else if (strncmp(AP.IGRAPH[i],"NE",2)==0){
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
      else if (strncmp(AP.IGRAPH[i],"CZ",2)==0){
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
      else if (strncmp(AP.IGRAPH[i],"NH1",3)==0){
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
      else if (strncmp(AP.IGRAPH[i],"HH11",4)==0){
	++numkai;
	kai[5][3]=i;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[5][0];
	atomdihedpairs[2][numkai*4-3]=kai[5][1];
	atomdihedpairs[2][numkai*4-2]=kai[5][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
      else if (strncmp(AP.IGRAPH[i],"HH21",4)==0){
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
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
	kai[2][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CG",2)==0){
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
      if (strncmp(AP.IGRAPH[i],"CD",2)==0){
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
      if (strncmp(AP.IGRAPH[i],"CE",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[2][0];
	atomdihedpairs[2][numkai*4-3]=kai[2][1];
	atomdihedpairs[2][numkai*4-2]=kai[2][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[4][1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"NZ",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[3][0];
	atomdihedpairs[2][numkai*4-3]=kai[3][1];
	atomdihedpairs[2][numkai*4-2]=kai[3][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[4][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"HZ1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[4][0];
	atomdihedpairs[2][numkai*4-3]=kai[4][1];
	atomdihedpairs[2][numkai*4-2]=kai[4][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"PHE",3)==0) {
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CD1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"TYR",3)==0) {
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CD1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CE1",3)==0){
    	kai[2][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CZ",2)==0){
	kai[2][1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"OH",2)==0){
	kai[2][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"HH",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[2][0];
	atomdihedpairs[2][numkai*4-3]=kai[2][1];
	atomdihedpairs[2][numkai*4-2]=kai[2][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"TRP",3)==0) {
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
	kai[1][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
	kai[1][1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CG",2)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
	kai[1][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CD1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[1][0];
	atomdihedpairs[2][numkai*4-3]=kai[1][1];
	atomdihedpairs[2][numkai*4-2]=kai[1][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }

    if (numres==1 && ACEFLAG==1 /*2014-03-23*/&& ACEFLAG2==1/*2014-03-23*/) {
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
	iomega[2]=i;
	atomdihedpairs[3][0]=ipsi[0];
	atomdihedpairs[3][1]=ipsi[1];
	atomdihedpairs[3][2]=ipsi[2];
	atomdihedpairs[3][3]=i;
      }
      else if (strncmp(AP.IGRAPH[i],"C",1)==0){
	iomega[3]=i;
	atomdihedpairs[3][4]=iomega[0];
	atomdihedpairs[3][5]=iomega[1];
	atomdihedpairs[3][6]=iomega[2];
	atomdihedpairs[3][7]=i;
	//	ACEFLAG=0; // 2014-03-23
	ACEFLAG2=0; // 2014-03-23
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

int PTL_res_ca(int numres) {
  int i;

  for (i=AP.IPRES[numres]-1;i<AP.IPRES[numres+1]-1;++i) {
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
      return i;
      break;
    }
  }

  if (i==AP.IPRES[numres+1]-1)
    return -1;
  else           // 2014-07-09
    return 0;    // 2014-07-09
}

/**********************************************************************/
/* int *PTL_ca_res_set(void) {					      */
/*   int i,ii;							      */
/*   int numatom;						      */
/*   int *index_ca_res;						      */
/* 								      */
/*   numatom=AP.NATOM;						      */
/*   index_ca_res=(int *)gcemalloc(sizeof(int)*1);		      */
/* 								      */
/*   ii=0;							      */
/*   for (i=0;i<numatom;++i) {					      */
/*     if (strncmp(AP.IGRAPH[i],"CA",2)==0) {			      */
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
/*   for (i=AP.IPRES[numres]-1;i<AP.IPRES[numres+1]-1;++i) {	      */
/*     if (strncmp(AP.IGRAPH[i],"CA",1)==0) {			      */
/*       return i;						      */
/*       break;							      */
/*     }							      */
/*   }								      */
/* 								      */
/*   if (i==AP.IPRES[numres+1]-1)				      */
/*     return -1;						      */
/* }								      */
/**********************************************************************/

int PTL_joinatomtores(int numatom, char LABERES[4]) {
  int i,j;

  for (i=0;i<AP.NRES;++i) {
    if (numatom >= AP.IPRES[i]-1 ) {
      if (i==(AP.NRES-1)){
	for (j=0;j<3;++j) {
	  LABERES[j]=AP.LABERES[i][j];
	}
	LABERES[3]='\0';
      }
      else {
	if (numatom < AP.IPRES[i+1]-1) {
	  for (j=0;j<3;++j) {
	    LABERES[j]=AP.LABERES[i][j];
	  }
	  LABERES[3]='\0';
	  break;
	}
      }
    }
  }
  
  return i;
}

int PTL_resnum(int numatom) {
  int i;

  if (numatom>AP.NATOM) return -1;

  for (i=0;i<AP.NRES;++i) {
    if (AP.IPRES[i]-1 > numatom)
      return i;
  }
  return AP.NRES;
}

int PTL_resnum2(int numatom) {
  int i;

  if (numatom>AP.NATOM) return -1;

  for (i=0;i<AP.NRES;++i) {
    //    if (AP.IPRES[i]-1 > numatom)
    if (AP.IPRES[i]-1 > numatom) //12-02-20
      return i;
  }
  return AP.NRES;
}

int PTL_canum_fromresnum(int numres) {
  int i;
  int fin;

  if (numres==AP.NRES-1) 
    fin=AP.NATOM;
  else 
    fin=AP.IPRES[numres+1];

  for (i=AP.IPRES[numres]-1;i<fin;++i) {
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
      return i;
      break;
    }
  }

  return -1;
}

int PTL_which_include(int numres,int *listres,int numlistres){
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
