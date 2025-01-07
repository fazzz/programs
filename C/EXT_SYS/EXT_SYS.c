#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>

#include "PTL.h"
#include "EF.h"

#include "EXT_SYS.h"

int extend_test_system(int num){
  int i,j,k;

  for (i=1;i<num;++i) for (j=0;j<AP.NATOM;++j) for (k=0;k<4;++k) AP.IGRAPH[i*AP.NATOM+j][k]=AP.IGRAPH[j][k];

  AP.CHRG=(double *)gcerealloc(AP.CHRG,sizeof(double)*AP.NATOM*num);
  for (i=1;i<num;++i) for (j=0;j<AP.NATOM;++j) AP.CHRG[i*AP.NATOM+j]=AP.CHRG[j];

  AP.AMASS=(double *)gcerealloc(AP.AMASS,sizeof(double)*AP.NATOM*num);
  for (i=1;i<num;++i) for (j=0;j<AP.NATOM;++j) AP.AMASS[i*AP.NATOM+j]=AP.AMASS[j];

  AP.IAC=(int *)gcerealloc(AP.IAC,sizeof(int)*AP.NATOM*num);
  for (i=1;i<num;++i) for (j=0;j<AP.NATOM;++j) AP.IAC[i*AP.NATOM+j]=AP.IAC[j];

  AP.NUMEX=(int *)gcerealloc(AP.NUMEX,sizeof(int)*AP.NATOM*num);
  for (i=1;i<num;++i) for (j=0;j<AP.NATOM;++j) AP.NUMEX[i*AP.NATOM+j]=AP.NUMEX[j];

  for (i=1;i<num;++i) for (j=0;j<AP.NRES;++j) for (k=0;k<4;++k) AP.LABERES[i*AP.NRES+j][k]=AP.LABERES[j][k];

  AP.IPRES=(int *)gcerealloc(AP.IPRES,sizeof(int)*AP.NRES*num);
  for (i=1;i<num;++i) for (j=0;j<AP.NRES;++j) AP.IPRES[i*AP.NRES+j]=AP.IPRES[j];

  if (AP.NBONH>0) {
    AP.BH=(int **)gcerealloc(AP.BH,sizeof(int *)*AP.NBONH*num);
    for (i=AP.NBONH;i<AP.NBONH*num;++i) AP.BH[i]=(int *)gcemalloc(sizeof(int)*3);
    for (i=1;i<num;++i){
      for (j=0;j<AP.NBONH;++j){
	AP.BH[i*AP.NBONH+j][0]=((i*AP.NATOM+(AP.BH[j][0]/3+1))-1)*3;
	AP.BH[i*AP.NBONH+j][1]=((i*AP.NATOM+(AP.BH[j][1]/3+1))-1)*3;
	AP.BH[i*AP.NBONH+j][2]=AP.BH[j][2];
      }
    }
  }

  if (AP.MBONA>0) {
    AP.BA=(int **)gcerealloc(AP.BA,sizeof(int *)*AP.MBONA*num);
    for (i=AP.MBONA;i<AP.MBONA*num;++i) AP.BA[i]=(int *)gcemalloc(sizeof(int)*3);
    for (i=1;i<num;++i){
      for (j=0;j<AP.MBONA;++j){
	AP.BA[i*AP.MBONA+j][0]=((i*AP.NATOM+(AP.BA[j][0]/3+1))-1)*3;
	AP.BA[i*AP.MBONA+j][1]=((i*AP.NATOM+(AP.BA[j][1]/3+1))-1)*3;
	AP.BA[i*AP.MBONA+j][2]=AP.BA[j][2];
      }
    }
  }

  if (AP.NTHETH>0) {
    AP.TH=(int **)gcerealloc(AP.TH,sizeof(int *)*AP.NTHETH*num);
    for (i=AP.NTHETH;i<AP.NTHETH*num;++i) AP.TH[i]=(int *)gcemalloc(sizeof(int)*4);
    for (i=1;i<num;++i){
      for (j=0;j<AP.NTHETH;++j){
	AP.TH[i*AP.NTHETH+j][0]=((i*AP.NATOM+(AP.TH[j][0]/3+1))-1)*3;
	AP.TH[i*AP.NTHETH+j][1]=((i*AP.NATOM+(AP.TH[j][1]/3+1))-1)*3;
	AP.TH[i*AP.NTHETH+j][2]=((i*AP.NATOM+(AP.TH[j][2]/3+1))-1)*3;
	AP.TH[i*AP.NTHETH+j][3]=AP.TH[j][3];
      }
    }
  }

  if (AP.MTHETA>0) {
    AP.TA=(int **)gcerealloc(AP.TA,sizeof(int *)*AP.MTHETA*num);
    for (i=AP.MTHETA;i<AP.MTHETA*num;++i) AP.TA[i]=(int *)gcemalloc(sizeof(int)*4);
    for (i=1;i<num;++i){
      for (j=0;j<AP.MTHETA;++j){
	AP.TA[i*AP.MTHETA+j][0]=((i*AP.NATOM+(AP.TA[j][0]/3+1))-1)*3;
	AP.TA[i*AP.MTHETA+j][1]=((i*AP.NATOM+(AP.TA[j][1]/3+1))-1)*3;
	AP.TA[i*AP.MTHETA+j][2]=((i*AP.NATOM+(AP.TA[j][2]/3+1))-1)*3;
	AP.TA[i*AP.MTHETA+j][3]=AP.TA[j][3];
      }
    }
  }

  if (AP.NPHIH>0) {
    AP.PH=(int **)gcerealloc(AP.PH,sizeof(int *)*AP.NPHIH*num);
    for (i=AP.NPHIH;i<AP.NPHIH*num;++i) AP.PH[i]=(int *)gcemalloc(sizeof(int)*5);
    for (i=1;i<num;++i){
      for (j=0;j<AP.NPHIH;++j){
	AP.PH[i*AP.NPHIH+j][0]=((i*AP.NATOM+(AP.PH[j][0]/3+1))-1)*3;
	AP.PH[i*AP.NPHIH+j][1]=((i*AP.NATOM+(AP.PH[j][1]/3+1))-1)*3;
	AP.PH[i*AP.NPHIH+j][2]=((i*AP.NATOM+(AP.PH[j][2]/3+1))-1)*3;
	AP.PH[i*AP.NPHIH+j][3]=((i*AP.NATOM+(AP.PH[j][3]/3+1))-1)*3;
	AP.PH[i*AP.NPHIH+j][4]=AP.PH[j][4];
      }
    }
  }

  if (AP.MPHIA>0) {
    AP.PA=(int **)gcerealloc(AP.PA,sizeof(int *)*AP.MPHIA*num);
    for (i=AP.MPHIA;i<AP.MPHIA*num;++i) AP.PA[i]=(int *)gcemalloc(sizeof(int)*5);
    for (i=1;i<num;++i){
      for (j=0;j<AP.MPHIA;++j){
	AP.PA[i*AP.MPHIA+j][0]=((i*AP.NATOM+(AP.PA[j][0]/3+1))-1)*3;
	AP.PA[i*AP.MPHIA+j][1]=((i*AP.NATOM+(AP.PA[j][1]/3+1))-1)*3;
	AP.PA[i*AP.MPHIA+j][2]=((i*AP.NATOM+(AP.PA[j][2]/3+1))-1)*3;
	AP.PA[i*AP.MPHIA+j][3]=((i*AP.NATOM+(AP.PA[j][3]/3+1))-1)*3;
	AP.PA[i*AP.MPHIA+j][4]=AP.PA[j][4];
      }
    }
  }

  if (AP.NEXT>0) {
    AP.NATEX=(int *)gcerealloc(AP.NATEX,sizeof(int)*AP.NEXT*num);
    for (i=1;i<num;++i){
      for (j=0;j<AP.NEXT;++j){
	AP.NATEX[i*AP.NEXT+j]=i*AP.NATOM+AP.NATEX[j];
      }
    }
  }

  AP.NATOM=AP.NATOM*num;

  AP.NBONH=AP.NBONH*num;
  AP.MBONA=AP.MBONA*num;
  AP.NTHETH=AP.NTHETH*num;
  AP.MTHETA=AP.MTHETA*num;
  AP.NPHIH=AP.NPHIH*num;
  AP.MPHIA=AP.MPHIA*num;

  AP.NRES=AP.NRES*num;
  AP.NBONA=AP.NBONA*num;
  AP.NTHETA=AP.NTHETA*num;
  AP.NPHIA=AP.NPHIA*num;

}

