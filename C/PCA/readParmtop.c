//#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include "ParmTop.h"

int readParmtop(FILE *parmfile/*, struct AmberParm AP*/){
  int i,j;
  char *line;
  size_t len=0;

  getline(&line,&len,parmfile);
  getline(&line,&len,parmfile);
  getline(&line,&len,parmfile);

  fscanf(parmfile,"%4s",&AP.ITITLE);

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  fscanf(parmfile,"%6d",&AP.NATOM);
  fscanf(parmfile,"%6d",&AP.NTYPES);
  fscanf(parmfile,"%6d",&AP.NBONH);
  fscanf(parmfile,"%6d",&AP.MBONA);
  fscanf(parmfile,"%6d",&AP.NTHETH);
  fscanf(parmfile,"%6d",&AP.MTHETA);
  fscanf(parmfile,"%6d",&AP.NPHIH);
  fscanf(parmfile,"%6d",&AP.MPHIA);
  fscanf(parmfile,"%6d",&AP.NHPARM);
  fscanf(parmfile,"%6d",&AP.NPARM);
  fscanf(parmfile,"%6d",&AP.NNB);
  AP.NEXT = AP.NNB;
  fscanf(parmfile,"%6d",&AP.NRES);
  fscanf(parmfile,"%6d",&AP.NBONA);
  fscanf(parmfile,"%6d",&AP.NTHETA);
  fscanf(parmfile,"%6d",&AP.NPHIA);
  fscanf(parmfile,"%6d",&AP.NUMBND);
  fscanf(parmfile,"%6d",&AP.NUMANG);
  fscanf(parmfile,"%6d",&AP.NPTRA);
  fscanf(parmfile,"%6d",&AP.NATYP);
  fscanf(parmfile,"%6d",&AP.NPHB);
  fscanf(parmfile,"%6d",&AP.IFPERT);
  fscanf(parmfile,"%6d",&AP.NBPER);
  fscanf(parmfile,"%6d",&AP.NGPER);
  fscanf(parmfile,"%6d",&AP.NDPER);
  fscanf(parmfile,"%6d",&AP.MBPER);
  fscanf(parmfile,"%6d",&AP.MGPER);
  fscanf(parmfile,"%6d",&AP.MDPER);
  fscanf(parmfile,"%6d",&AP.IFBOX);
  fscanf(parmfile,"%6d",&AP.NMXPS);
  fscanf(parmfile,"%6d",&AP.IFCAP);
  fscanf(parmfile,"%6d",&AP.NEXTRA);

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NATOM;++i){
    fscanf(parmfile,"%4s",&AP.IGRAPH[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NATOM;++i){
    fscanf(parmfile,"%lf",&AP.CHRG[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NATOM;++i){
    fscanf(parmfile,"%lf",&AP.AMASS[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NATOM;++i){
    fscanf(parmfile,"%d",&AP.IAC[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NATOM;++i){
    fscanf(parmfile,"%lf",&AP.NUMEX[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<(AP.NTYPES)*(AP.NTYPES);++i){
    fscanf(parmfile,"%d",&AP.ICO[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NRES;++i){
    fscanf(parmfile,"%4s",&AP.ICO[i]);
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
  for (i=0;i<AP.NRES;++i){
    fscanf(parmfile,"%4s",&AP.IPRES[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NUMBND;++i){
    fscanf(parmfile,"%lf",AP.RK[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NUMBND;++i){
    fscanf(parmfile,"%lf",&AP.REQ[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NUMANG;++i){
    fscanf(parmfile,"%lf",&AP.TK[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NUMANG;++i){
    fscanf(parmfile,"%lf",&AP.TEQ[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NPTRA;++i){
    fscanf(parmfile,"%lf",&AP.PK[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NPTRA;++i){
    fscanf(parmfile,"%lf",&AP.PN[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NPTRA;++i){
    fscanf(parmfile,"%lf",&AP.PHASE[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NATYP;++i){
    fscanf(parmfile,"%lf",&AP.SOLTY[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<(AP.NTYPES)*(AP.NTYPES+1)/2;++i){
    fscanf(parmfile,"%lf",&AP.CN1[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<(AP.NTYPES)*(AP.NTYPES+1)/2;++i){
    fscanf(parmfile,"%lf",&AP.CN2[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NBONH;++i){
    for (j=0;j<3;++j){
      fscanf(parmfile,"%d",&AP.BH[i][j]);
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NBONA;++i){
    for (j=0;j<3;++j){
      fscanf(parmfile,"%d",&AP.BA[i][j]);
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NTHETH;++i){
    for (j=0;j<4;++j){
      fscanf(parmfile,"%d",&AP.TH[i][j]);
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NTHETA;++i){
    for (j=0;j<4;++j){
      fscanf(parmfile,"%d",&AP.TA[i][j]);
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NPHIH;++i){
    for (j=0;j<5;++j){
      fscanf(parmfile,"%d",&AP.TA[i][j]);
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NPHIA;++i){
    for (j=0;j<5;++j){
      fscanf(parmfile,"%d",&AP.TA[i][j]);
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NEXT;++i){
    fscanf(parmfile,"%lf",&AP.NATEX[i]);
    }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NEXT;++i){
    fscanf(parmfile,"%lf",&AP.NATEX[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NPHB;++i){
    fscanf(parmfile,"%lf",&AP.ASOL[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NPHB;++i){
    fscanf(parmfile,"%lf",&AP.BSOL[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
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
  for (i=0;i<AP.NATOM;++i){
    fscanf(parmfile,"%d",&AP.JOIN[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NATOM;++i){
    fscanf(parmfile,"%d",&AP.IROTAT[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  if (AP.IFBOX > 0) {
    fscanf(parmfile,"%d",&AP.IPTRES);

    fscanf(parmfile,"%d",&AP.NSPM);

    fscanf(parmfile,"%d",&AP.NSPSOL);

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
    for (i=0;i<AP.NBPER;++i){
      for (j=0;j<2;++j){
	fscanf(parmfile,"%d",&AP.BPER[i][j]);
      }
    }

    for (i=0;i<AP.NBPER*2;++i){
      fscanf(parmfile,"%d",&AP.ICBPER[i]);
    }

    for (i=0;i<AP.NGPER;++i){
      for (j=0;j<3;++j){
	fscanf(parmfile,"%d",&AP.TPER[i][j]);
      }
    }

    for (i=0;i<AP.NGPER*2;++i){
      fscanf(parmfile,"%d",&AP.ICTPER[i]);
    }

    for (i=0;i<AP.NDPER;++i){
      for (j=0;j<4;++j){
	fscanf(parmfile,"%d",&AP.PPER[i][j]);
      }
    }

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

    for (i=0;i<AP.NATOM;++i){
      fscanf(parmfile,"%lf",&AP.IAPER[i]);
    }

    for (i=0;i<AP.NATOM;++i){
      fscanf(parmfile,"%d",&AP.IACPER[i]);
    }

    for (i=0;i<AP.NATOM;++i){
      fscanf(parmfile,"%lf",&AP.CGPER[i]);
    }
  } 

  /* if (AP.IPOL == 1) {
    fscanf(parmfile,"%18.8lf",&AP.ATPOL[i]);

    fscanf(parmfile,"%18.8lf",&AP.ATPOL1[i]);
    }*/

  if (AP.NPARM == 1) {
    fscanf(parmfile,"%d",&AP.NLES_NTYP);

    for (i=0;i<AP.NATOM;++i){
      fscanf(parmfile,"%lf",&AP.LES_TYPE[i]);
    }

    for (i=0;i<(AP.NLES_NTYP)*(AP.NLES_NTYP);++i){
      fscanf(parmfile,"%lf",&AP.LES_FAC[i]);
    }

    for (i=0;i<AP.NATOM;++i){
      fscanf(parmfile,"%lf",&AP.LES_CNUM[i]);
    }

    for (i=0;i<AP.NATOM;++i){
      fscanf(parmfile,"%lf",&AP.LES_ID[i]);
    }

  }

}

