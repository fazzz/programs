
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GOLMAA_PROTEINS2008_Amber_hybrid.h"
#include "FFL.h"

#include "EF.h"
#include "MB.h"
#include "TOPO.h"
#include "LA.h"
#include "mymath.h"

#include "PTL.h"

#define UNITIN 4.184070*100.0

double GOLMAA_PROTEINS2008_Amber_hybrid_ff_calcff_b(double *crd, int numatom,
						    struct potential_GOLMAA_PROTEINS2008 *ene) {
  int i,j;

  (*ene).p_d_t=GOLMAA_PROTEINS2008_Amber_hybrid_calcDIHE_force_Cartesian((*ene).f_d,crd);
  (*ene).p_a_t=GOLMAA_PROTEINS2008_Amber_hybrid_calcANGLE_force_Cartesian((*ene).f_a,crd);
  (*ene).p_b_t=GOLMAA_PROTEINS2008_Amber_hybrid_calcBOND_force_Cartesian((*ene).f_b,crd);

    //  (*ene).p_b_t=GOLMAA_PROTEINS2008_ff_calcBOND(crd,numatom,(*ene).p_b,(*ene).f_b,(*ene).Kb,(*ene).bon_equ,(*ene).pairs_bond,(*ene).num_bond);
  //  (*ene).p_a_t=GOLMAA_PROTEINS2008_ff_calcANGLE(crd,numatom,(*ene).p_a,(*ene).f_a,(*ene).Ka,(*ene).ang_equ,(*ene).pairs_angl,(*ene).num_angl);
  //  (*ene).p_d_t=GOLMAA_PROTEINS2008_ff_calcDIHE(crd,numatom,(*ene).p_d,(*ene).f_d,(*ene).Kd1,(*ene).Kd2,(*ene).Ki,(*ene).dih_equ,(*ene).pairs_dihe,(*ene).num_dihe,(*ene).impindex);


  GOLMAA_PROTEINS2008_ff_calcff_nonlocal_b(crd,numatom,
  					   (*ene).NC_index,(*ene).numNC,(*ene).NotNC_index,(*ene).numNotNC,
  					   (*ene).ALJ_natatt,(*ene).BLJ_natatt,(*ene).ep_natatt,(*ene).ALJ_repul,
  					   &((*ene).p_natatt_t),&((*ene).p_repul_t),(*ene).f_natatt,(*ene).f_repul);

  (*ene).p_t=(*ene).p_natatt_t+(*ene).p_repul_t+(*ene).p_b_t+(*ene).p_a_t+(*ene).p_d_t;

  for (i=0;i<numatom;++i) 
    for (j=0;j<3;++j) 
      (*ene).f_t[i][j]=(*ene).f_b[i][j]+(*ene).f_a[i][j]+(*ene).f_d[i][j]+(*ene).f_natatt[i][j]+(*ene).f_repul[i][j];

  return (*ene).p_t;
}

double GOLMAA_PROTEINS2008_Amber_hybrid_calcDIHE_force_Cartesian(double **f_d,double *cord) {
  int i,j,k,l,flag;
  int dtype;

  double fa,fb[3],fc[3];
  double *n1,*n2,ln1,ln2;
  double vij[3],vkj[3],vki[3],vjl[3],vji[3],vik[3],vkl[3],vlk[3];
  double op1[3],op2[3],op3[3],op4[3],op5[3],op6[3];
  
  double atom[4][3];
  double cosdih,sindih;
  double dihedang;

  double p_d_t=0.0;

  n1=(double *)gcemalloc(sizeof(double)*3);
  n2=(double *)gcemalloc(sizeof(double)*3);

  for (i=0;i<AP.NATOM;++i) for (j=0;j<3;++j) f_d[i][j] = 0.0;
  
  for (i=0;i<AP.MPHIA;++i) {
    dtype = AP.PA[i][4]-1;
    for (j=0;j<4;++j) for (k=0;k<3;++k)	atom[j][k]=cord[abs(AP.PA[i][j])+k];

    for (j=0;j<3;++j) {
      vij[j] = atom[1][j]-atom[0][j];
      vkj[j] = atom[1][j]-atom[2][j];
      vki[j] = atom[0][j]-atom[2][j];
      vjl[j] = atom[3][j]-atom[1][j];
      vji[j] = atom[0][j]-atom[1][j];
      vik[j] = atom[2][j]-atom[0][j];
      vkl[j] = atom[3][j]-atom[2][j];
      vlk[j] = atom[2][j]-atom[3][j];
    }

    outprod(vij,vkj,n1);
    outprod(vkj,vkl,n2);
    ln1=sqrt(inprod(n1,n1,3));
    ln2=sqrt(inprod(n2,n2,3));
  
    csdih(atom[0],atom[1],atom[2],atom[3],&cosdih,&sindih);

    if (AP.PN[dtype]==1) fa=-AP.PK[dtype]*AP.PN[dtype]*cos(AP.PHASE[dtype]);
    else if (AP.PN[dtype]==2) fa=-AP.PK[dtype]*AP.PN[dtype]*2.0*cosdih*cos(AP.PHASE[dtype]);
    else if (AP.PN[dtype]==3) fa=-AP.PK[dtype]*AP.PN[dtype]*(-4.0*sindih*sindih+3.0)*cos(AP.PHASE[dtype]);
    else if (AP.PN[dtype]==4) fa=-AP.PK[dtype]*AP.PN[dtype]*4.0*(cosdih*(2.0*cosdih*cosdih-1.0))*cos(AP.PHASE[dtype]);
    else {
      printf("error:periodicity must be 1~4\n");
      exit(1);
    }

    for (j=0;j<3;++j) fb[j]=(n2[j]/ln2-cosdih*n1[j]/ln1)/ln1;
    for (j=0;j<3;++j) fc[j]=(n1[j]/ln1-cosdih*n2[j]/ln2)/ln2;

    outprod(fb,vkj,op1);
    //    outprod(fc,vki,op2); // 1111
    //    outprod(fb,vik,op3); // 1111
    outprod(fb,vki,op2);       // 1111
    outprod(fc,vlk,op3);       // 1111
    outprod(fb,vij,op4);
    outprod(fc,vjl,op5);
    outprod(fc,vkj,op6);

    for (j=0;j<3;++j) {
      f_d[abs(AP.PA[i][0])/3][j] += fa*op1[j]*UNITIN;
      f_d[abs(AP.PA[i][1])/3][j] += fa*(-op2[j]+op3[j])*UNITIN;
      f_d[abs(AP.PA[i][2])/3][j] += fa*(-op4[j]+op5[j])*UNITIN;
      f_d[abs(AP.PA[i][3])/3][j] += fa*op6[j]*UNITIN;
    }

    dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
    p_d_t += AP.PK[dtype]*(1.0+cos(AP.PN[dtype]*dihedang-AP.PHASE[dtype]));
  }

  return p_d_t;
}

double GOLMAA_PROTEINS2008_Amber_hybrid_calcANGLE_force_Cartesian(double **f_a,double *cord){
  int i,j,k,l;
  int numatom;
  int type;
  double kang,ang_eq;

  double atom[3][3];
  double *f_temp;

  double p_a=0.0;
  double angle;

  numatom=AP.NATOM;
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) f_a[i][j] = 0.0;
  f_temp=(double *)gcemalloc(sizeof(double)*3*3);

  for (i=0;i<AP.MTHETA;++i) {
    for (j=0;j<3;++j)  for (k=0;k<3;++k) atom[j][k]=cord[abs(AP.TA[i][j])+k];
    type = AP.TA[i][3]-1;
    kang = AP.TK[type];
    ang_eq = AP.TEQ[type];
    calcANGKE_force2(atom[0],atom[1],atom[2],kang,ang_eq,f_temp);
    for (j=0;j<3;++j) {
      f_a[abs(AP.TA[i][0])/3][j] += f_temp[j];
      f_a[abs(AP.TA[i][1])/3][j] += f_temp[3+j];
      f_a[abs(AP.TA[i][2])/3][j] += f_temp[6+j];
    }
    angle = pick_angle(atom[0],atom[1],atom[2],0,0.0);
    p_a += AP.TK[type]*(angle-AP.TEQ[type])*(angle-AP.TEQ[type]);

  }

  return p_a;
}

double calcANGKE_force2(double atomi[3],double atomj[3],double atomk[3],double kang,double ang_eq,double *f) {
  int i,j,k;
  double lenij,lenkj;
  double vij[3],vkj[3];
  double cosijk,angijk;
  double f1,f2;

  lenij = len(atomi,atomj);
  lenkj = len(atomk,atomj);
  for (j=0;j<3;++j) {
    vij[j]=atomj[j]-atomi[j];
    vkj[j]=atomj[j]-atomk[j];
  }
  cosijk=inprod(vij,vkj,3);
  cosijk=cosijk/lenij/lenkj;
  angijk = acos(cosijk);
  //  angijk=ang(atomi,atomj,atomk);
  //  cosijk=cos(angijk);

  for (j=0;j<3;++j) {
    /*******************************************************************************************/
    /* f1 = -kang*(angijk-ang_eq)/(lenij*sin(angijk))*(vkj[j]/lenkj-cosijk*vij[j]/lenij)*UNITIN; */
    /* f2 = -kang*(angijk-ang_eq)/(lenkj*sin(angijk))*(vij[j]/lenij-cosijk*vkj[j]/lenkj)*UNITIN; */
    /*******************************************************************************************/
    f1 = -2.0*kang*(angijk-ang_eq)/(lenij*sin(angijk))*(vkj[j]/lenkj-cosijk*vij[j]/lenij)*UNITIN;
    f2 = -2.0*kang*(angijk-ang_eq)/(lenkj*sin(angijk))*(vij[j]/lenij-cosijk*vkj[j]/lenkj)*UNITIN;

    f[j] = f1;
    f[6+j] = f2;
    f[3+j] = -f1-f2;
  }
}

double GOLMAA_PROTEINS2008_Amber_hybrid_calcBOND_force_Cartesian(double **f_b,double *cord){
  int i,j,k;
  int numatom,type;
  double f,lenij;
  double atom[2][3];

  double p_b_t=0.0;

  numatom=AP.NATOM;
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) f_b[i][j] = 0.0;
  
  for (i=0;i<AP.MBONA;++i) {
    type = AP.BA[i][2]-1;
    for (j=0;j<2;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(AP.BA[i][j])+k];

    lenij=0.0; for (j=0;j<3;++j) lenij+=(atom[0][j]-atom[1][j])*(atom[0][j]-atom[1][j]);  lenij=sqrt(lenij);

    for (j=0;j<3;++j) {
      f = /*-*/2.0*AP.RK[type]*(lenij-AP.REQ[type])*(atom[1][j]-atom[0][j])/lenij*UNITIN;
      f_b[abs(AP.BA[i][0])/3][j] += f;
      f_b[abs(AP.BA[i][1])/3][j] += -f;
    }

    p_b_t += AP.RK[type]*(lenij-AP.REQ[type])*(lenij-AP.REQ[type]);
  }

  return p_b_t;
}
