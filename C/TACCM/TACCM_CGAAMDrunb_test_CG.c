
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"
#include "FFLc.h"
#include "PTLb.h"

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"

#include "TACCM_CGAAMDrun_test_CG.h"
#include "TACCM_MDrun.h"
#include "TACCM_MD.h"

#include "MD_NHC_MP1996.h"
#include "MDrun.h"
#include "MD.h"

#include "MB.h"
#include "LA.h"
#include "TOPO.h"
#include "mymath.h"

double ffL_calcffandforce_AA(double *crd, int numatom,struct potential *ene,struct force *f) {
  int i;
  int numnb,num14;
  double *n_d;

  numnb=(*ene).parm.numnb;
  num14=(*ene).parm.num14;

  ffL_calcDIHE((*ene).p_d,n_d,crd,1,0,0);
  ffL_calcDIHE_force_Cartesian((*f).f_d,crd); // 1111

  ffL_calcANGLE((*ene).p_a,crd);
  ffL_calcANGLE_force_Cartesian((*f).f_a,crd);

  ffL_calcBOND((*ene).p_b,crd);
  ffL_calcBOND_force_Cartesian((*f).f_b,crd);
  
  (*ene).p_t=0.0;
  (*ene).p_e_t=0.0;
  (*ene).p_LJ_t=0.0;
  (*ene).p_e_14_t=0.0;
  (*ene).p_LJ_14_t=0.0;
  (*ene).p_d_t=0.0;
  (*ene).p_a_t=0.0;
  (*ene).p_b_t=0.0;
  for (i=0;i<numatom*3;++i) (*f).f_t[i]=0.0;

  for (i=0;i<AP.NPHIH+AP.MPHIA;++i) {
    (*ene).p_t+=(*ene).p_d[i];
    (*ene).p_d_t+=(*ene).p_d[i];
  }
  for (i=0;i<AP.NTHETH+AP.MTHETA;++i) {
    (*ene).p_t+=(*ene).p_a[i];
    (*ene).p_a_t+=(*ene).p_a[i];             // 0911
  }                                          // 0911
  for (i=0;i<AP.NBONH+AP.MBONA;++i) {        // 0911
    (*ene).p_t+=(*ene).p_b[i];               // 0911
    (*ene).p_b_t+=(*ene).p_b[i];
  }
  for (i=0;i<numatom*3;++i) {
    (*f).f_e_14[i]=1.0/1.2*(*f).f_e_14[i];
    (*f).f_LJ_14[i]=0.5*(*f).f_LJ_14[i];
  }

  for (i=0;i<numatom*3;++i)
    (*f).f_t[i]+=(*f).f_d[i]+(*f).f_a[i]+(*f).f_b[i];
  
  return (*ene).p_t;
}

double TACCM_MD_Propagetor_NH_MP1998_CG_test(double *crd,double *vel,double *mass,
					     double *zeta,double *V_zeta,double Q,
					     double NfKT,int numatom,double *KE,double *KEv,double *PEv,
					     double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
					     double parameterCG, struct potential *e, struct force *f,
					     double *Z,  int numZ,double *theta,double Kapa,
					     int **pairs, double *PEZ,double pi) {
  int i,j,k;
  double *frc;
  double **frcZ;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  frcZ=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) frcZ[i]=(double *)gcemalloc(sizeof(double)*3);

  TACCM_CTheta(crd,numatom,theta,numZ,pairs,pi);
  *PEZ=TACCM_calc_eff_FF_MD(crd,numatom,theta,Z,numZ,Kapa,frcZ,pairs,pi);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+frcZ[i][j];
  
  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];

  ffL_calcffandforce_CG(crd,numatom,e,f,parameterCG);
  TACCM_CTheta(crd,numatom,theta,numZ,pairs,pi);
  *PEZ=TACCM_calc_eff_FF_MD(crd,numatom,theta,Z,numZ,Kapa,frcZ,pairs,pi);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+frcZ[i][j];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  *KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) *KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return 0.0;
}

double ffL_calcffandforce_CG(double *crd, int numatom,struct potential *ene,struct force *f,double parameterCG) {
  int i;
  int numnb,num14;
  double *n_d;

  numnb=(*ene).parm.numnb;
  num14=(*ene).parm.num14;

  ffL_calcDIHE_CG((*ene).p_d,n_d,crd,1,0,0,parameterCG);
  ffL_calcDIHE_force_Cartesian_CG((*f).f_d,crd,parameterCG); // 1111

  ffL_calcANGLE((*ene).p_a,crd);
  ffL_calcANGLE_force_Cartesian((*f).f_a,crd);

  ffL_calcBOND((*ene).p_b,crd);
  ffL_calcBOND_force_Cartesian((*f).f_b,crd);
  
  (*ene).p_t=0.0;
  (*ene).p_e_t=0.0;
  (*ene).p_LJ_t=0.0;
  (*ene).p_e_14_t=0.0;
  (*ene).p_LJ_14_t=0.0;
  (*ene).p_d_t=0.0;
  (*ene).p_a_t=0.0;
  (*ene).p_b_t=0.0;
  for (i=0;i<numatom*3;++i) (*f).f_t[i]=0.0;

  for (i=0;i<AP.NPHIH+AP.MPHIA;++i) {
    (*ene).p_t+=(*ene).p_d[i];
    (*ene).p_d_t+=(*ene).p_d[i];
  }
  for (i=0;i<AP.NTHETH+AP.MTHETA;++i) {
    (*ene).p_t+=(*ene).p_a[i];
    (*ene).p_a_t+=(*ene).p_a[i];             // 0911
  }                                          // 0911
  for (i=0;i<AP.NBONH+AP.MBONA;++i) {        // 0911
    (*ene).p_t+=(*ene).p_b[i];               // 0911
    (*ene).p_b_t+=(*ene).p_b[i];
  }
  for (i=0;i<numatom*3;++i) {
    (*f).f_e_14[i]=1.0/1.2*(*f).f_e_14[i];
    (*f).f_LJ_14[i]=0.5*(*f).f_LJ_14[i];
  }

  for (i=0;i<numatom*3;++i)
    (*f).f_t[i]+=(*f).f_d[i]+(*f).f_a[i]+(*f).f_b[i];
  
  return (*ene).p_t;
}

int ffL_calcDIHE_CG(double *p_d,
		    double *n_d,
		    double *cord,
		    int flagp, int flagf, int flaginp,double parameterCG) {
  int i,j,k,l,flag;
  int dtype;
  
  double atom[4][3];
  double dihedang;
  
  for (i=0;i<AP.NPHIH;++i){
    p_d[i] = 0.0;
    if (flagf==1)
      n_d[i] = 0.0;
  }
  
  for (i=0;i<AP.MPHIA;++i){
    p_d[i+AP.NPHIH] = 0.0;
    if (flagf==1)
      n_d[i+AP.NPHIH] = 0.0;
  }
  
  for (i=0;i<AP.NPHIH;++i) {
    dtype = AP.PH[i][4]-1;
    flag=1;
    /************************************/
    /* for (j=0;j<inpnum;++j) {	        */
    /*   if (i == inpindex[j]) {        */
    /* 	flag=0;break;		        */
    /*   }			        */
    /* }			        */
    /************************************/
    if (flag==1 || flaginp==0) {
      for (j=0;j<4;++j) {
  	for (k=0;k<3;++k) {
  	  atom[j][k]=cord[abs(AP.PH[i][j])+k];
  	}
      }
  
      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      p_d[i] = AP.PK[dtype]*(1.0+cos(AP.PN[dtype]*dihedang-AP.PHASE[dtype]))*parameterCG;
      if (flagf==1)
  	n_d[i] = -AP.PK[dtype]*4.18407*100.0*(sin(AP.PN[dtype]*dihedang-AP.PHASE[dtype])*AP.PN[dtype])*parameterCG;
    }
  }
  
  for (i=0;i<AP.MPHIA;++i) {
    dtype = AP.PA[i][4]-1;
    flag=1;
    /************************************/
    /* for (j=0;j<inpnum;++j) {	        */
    /*   if (i == inpindex[j]) {        */
    /* 	flag=0;break;		        */
    /*   }			        */
    /* }			        */
    /************************************/
    if (flag==1 || flaginp==0) {
      for (j=0;j<4;++j) {
  	for (k=0;k<3;++k) {
  	  atom[j][k]=cord[abs(AP.PA[i][j])+k];
  	}
      }
  
      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      p_d[i+AP.NPHIH] = AP.PK[dtype]*(1.0+cos(AP.PN[dtype]*dihedang-AP.PHASE[dtype]))*parameterCG;
      if (flagf==1)
  	n_d[i+AP.NPHIH] = -AP.PK[dtype]*4.18407*100.0*(sin(AP.PN[dtype]*dihedang-AP.PHASE[dtype])*AP.PN[dtype])*parameterCG;
    }
  }
}

int ffL_calcDIHE_force_Cartesian_CG(double *f_d,double *cord,double parameterCG) {
  int i,j,k,l,flag;
  int dtype;

  double fa,fb[3],fc[3];
  double *n1,*n2,ln1,ln2;
  double vij[3],vkj[3],vki[3],vjl[3],vji[3],vik[3],vkl[3],vlk[3];
  double op1[3],op2[3],op3[3],op4[3],op5[3],op6[3];
  
  double atom[4][3];
  double cosdih,sindih;
  double dihedang;

  n1=(double *)gcemalloc(sizeof(double)*3);
  n2=(double *)gcemalloc(sizeof(double)*3);

  //  for (i=0;i<(AP.NPHIH+AP.MPHIA)*3;++i) f_d[i] = 0.0;
  for (i=0;i<AP.NATOM*3;++i) f_d[i] = 0.0;
  
  for (i=0;i<AP.NPHIH;++i) {
    dtype = AP.PH[i][4]-1;
    for (j=0;j<4;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(AP.PH[i][j])+k];

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

    if (AP.PN[dtype]==1) fa=-AP.PK[dtype]*AP.PN[dtype]*cos(AP.PHASE[dtype])*parameterCG;
    else if (AP.PN[dtype]==2) fa=-AP.PK[dtype]*AP.PN[dtype]*2.0*cosdih*cos(AP.PHASE[dtype])*parameterCG;
    else if (AP.PN[dtype]==3) fa=-AP.PK[dtype]*AP.PN[dtype]*(-4.0*sindih*sindih+3.0)*cos(AP.PHASE[dtype])*parameterCG;
    else if (AP.PN[dtype]==4) fa=-AP.PK[dtype]*AP.PN[dtype]*4.0*(cosdih*(2.0*cosdih*cosdih-1.0))*cos(AP.PHASE[dtype])*parameterCG;
    else {
      printf("error:phase must be 1~4\n");
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
      f_d[abs(AP.PH[i][0])+j] += fa*op1[j]*UNIT;
      f_d[abs(AP.PH[i][1])+j] += fa*(-op2[j]+op3[j])*UNIT;
      f_d[abs(AP.PH[i][2])+j] += fa*(-op4[j]+op5[j])*UNIT;
      f_d[abs(AP.PH[i][3])+j] += fa*op6[j]*UNIT;
    }
  }
  
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

    if (AP.PN[dtype]==1) fa=-AP.PK[dtype]*AP.PN[dtype]*cos(AP.PHASE[dtype])*parameterCG;
    else if (AP.PN[dtype]==2) fa=-AP.PK[dtype]*AP.PN[dtype]*2.0*cosdih*cos(AP.PHASE[dtype])*parameterCG;
    else if (AP.PN[dtype]==3) fa=-AP.PK[dtype]*AP.PN[dtype]*(-4.0*sindih*sindih+3.0)*cos(AP.PHASE[dtype])*parameterCG;
    else if (AP.PN[dtype]==4) fa=-AP.PK[dtype]*AP.PN[dtype]*4.0*(cosdih*(2.0*cosdih*cosdih-1.0))*cos(AP.PHASE[dtype])*parameterCG;
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
      f_d[abs(AP.PA[i][0])+j] += fa*op1[j]*UNIT;
      f_d[abs(AP.PA[i][1])+j] += fa*(-op2[j]+op3[j])*UNIT;
      f_d[abs(AP.PA[i][2])+j] += fa*(-op4[j]+op5[j])*UNIT;
      f_d[abs(AP.PA[i][3])+j] += fa*op6[j]*UNIT;
    }
    // for debug
    //    FILE *db;
    //    db=efopen("db_z_fa","a");
    //    fprintf(db,"%e\n",fa);
    //    fclose(db);
    //    FILE *db2;
    //    db2=efopen("db_z_sd_cd","a");
    //    fprintf(db2,"%e %e\n",sindih,cosdih);
    //    fclose(db2);
    // for debug
  }
}
