#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>

#include "PTL.h"
#include "EF.h"

int readParmtopL(FILE *parmfile){
  int i,j;
  char *line;
  size_t len=0;
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
    fscanf(parmfile,"%4s",&AP.IGRAPH[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  AP.CHRG=(double *)gcemalloc(sizeof(double)*AP.NATOM);
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
    fscanf(parmfile,"%4s",&AP.LABERES[i]);
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
    fscanf(parmfile,"%4s",&AP.ISYMBL[i]);
  }
  
  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NATOM;++i){
    fscanf(parmfile,"%4s",&AP.ITREE[i]);
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
      fscanf(parmfile,"%12d",&AP.BOX[i]);
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
    for (i=0;i<AP.NDPER;++i)  AP.PPER[i]=(int *)gcemalloc(sizeof(int)*4); //AP.PPER=(int *)gcemalloc(sizeof(int)*4);
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
      fscanf(parmfile,"%4s",&AP.LABERES[i]);
    }

    for (i=0;i<AP.NATOM;++i){
      fscanf(parmfile,"%4s",&AP.IGRPER[i]);
    }

    for (i=0;i<AP.NATOM;++i){
      fscanf(parmfile,"%4s",&AP.ISMPER[i]);
    }

    for (i=0;i<AP.NATOM;++i){
      fscanf(parmfile,"%lf",&AP.ALMPER[i]);
    }

    AP.IAPER=(/*int*/double *)gcemalloc(sizeof(/*int*/double)*AP.NATOM);
    for (i=0;i<AP.NATOM;++i){
      fscanf(parmfile,"%lf",&AP.IAPER[i]);
    }

    AP.IACPER=(int *)gcemalloc(sizeof(int)*AP.NATOM);
    for (i=0;i<AP.NATOM;++i){
      fscanf(parmfile,"%d",&AP.IACPER[i]);
    }

    AP.CGPER=(/*int*/double *)gcemalloc(sizeof(/*int*/double)*AP.NATOM);
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
      fscanf(parmfile,"%lf",&AP.LES_TYPE[i]);
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

}

