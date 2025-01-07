
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "FFL.h"
#include "MB.h"
#include "PTL.h"
#include "LA.h"
#include "TOPO.h"
#include "EF.h"

#define YES 1
#define NO  0

int ffL_set_calcffandforce(struct potential *ene, struct force *f){
  int i,numatom,fd;
  int numnb,num14;
  char dummy;

  ffL_set_non_bonding_index_1(&numnb,&num14);
  (*ene).parm.numnb=numnb;
  (*ene).parm.num14=num14;
  (*ene).parm.indexnb=(int *)gcemalloc(sizeof(int)*numnb*2);
  (*ene).parm.index14=(int *)gcemalloc(sizeof(int)*num14*2);
  ffL_set_non_bonding_index_2((*ene).parm.indexnb,(*ene).parm.index14);

  numatom=AP.NATOM;
  (*ene).parm.ele=(double *)gcemalloc(sizeof(double)*numatom);
  (*ene).parm.ALJ=(double *)gcemalloc(sizeof(double)*numatom*numatom);
  (*ene).parm.BLJ=(double *)gcemalloc(sizeof(double)*numatom*numatom);
  (*ene).p_e=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_LJ=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_e_14=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_LJ_14=(double *)gcemalloc(sizeof(double)*numnb);
  (*ene).p_d=(double *)gcemalloc(sizeof(double)*(AP.NPHIH+AP.MPHIA));
  (*ene).p_a=(double *)gcemalloc(sizeof(double)*(AP.NTHETH+AP.MTHETA));
  (*ene).p_b=(double *)gcemalloc(sizeof(double)*(AP.NBONH+AP.MBONA));

  ffL_set_NB_PARM((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,numatom);

  (*f).f_t=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_e=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_LJ=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_e_14=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_LJ_14=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_d=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_a=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_b=(double *)gcemalloc(sizeof(double)*numatom*3);

}

void ffL_set_NB_PARM(double *ele, double *ALJ, double *BLJ, int numatom){
  int i,j,index;
  
  for (i=0;i<numatom;++i) {
    ele[i]= AP.CHRG[i];
    for (j=0;j<numatom;++j) {
      index=AP.ICO[(AP.IAC[i]-1)*AP.NTYPES+(AP.IAC[j]-1)]-1;
      ALJ[i*numatom+j] = AP.CN1[index];
      BLJ[i*numatom+j] = AP.CN2[index];
    }
  }
}

void ffL_set_non_bonding_index_1(int *numindexnb, int *numindex14) {
  int i,j,k,l;
  int numatom;

  int gnumnb=0;
  int gnum14=0;

  numatom=AP.NATOM;

  gnumnb=0;
  gnum14=0;

  for(i=0;i<numatom;++i) {
    for(j=i+1;j<numatom;++j) {
      if (which_calc_nb(i,j) == YES)
	++gnumnb;
      else if(which_calc_14_nb(i,j)==YES)
      	++gnum14;
    }
  }

  (*numindexnb)=gnumnb;
  (*numindex14)=gnum14;

}

void ffL_set_non_bonding_index_2(int *gindexnb,int *gindex14) {
  int i,j,k,l;
  int numatom;

  int gnumnb=0;
  int gnum14=0;

  numatom=AP.NATOM;

  k=0;
  l=0;
  for(i=0;i<numatom;++i)
    for(j=i+1;j<numatom;++j)
      if (which_calc_nb(i,j) == YES) {
  	gindexnb[k*2]=i;
  	gindexnb[k*2+1]=j;
  	++k;
      }
      else if(which_calc_14_nb(i,j)==YES){
  	gindex14[l*2]=i;
  	gindex14[l*2+1]=j;
  	++l;
      }

}

int which_calc_nb(int atomi,int atomj) {
  int k;
  int sum=0;

  for(k=0;k<atomi;++k) sum+=AP.NUMEX[k];

  for(k=sum;k<sum+AP.NUMEX[atomi];++k){
    if (AP.NATEX[k]-1 == atomj)
      return NO;
  }

  return YES;
}

int which_calc_14_nb(int atomi,int atomj) {
  int k,l;
  int flag;

  if ((flag=which_calc_nb(atomi,atomj))==YES)
    return NO;
  else {
    for(l=0;l<AP.NBONH;++l){
      if ((AP.BH[l][0])/3 == atomi ){
  	if ((AP.BH[l][1])/3 == atomj)
  	  return NO;
      }
      else if ((AP.BH[l][1])/3 == atomi){
  	if ((AP.BH[l][0])/3 == atomj)
  	  return NO;
      }
      else if ((AP.BH[l][0])/3 == atomj){
  	if ((AP.BH[l][1])/3 == atomi)
  	  return NO;
      }
      else if ((AP.BH[l][1])/3 == atomj){
  	if ((AP.BH[l][0])/3 == atomi)
  	  return NO;
      }
    }
  
    for(l=0;l<AP.MBONA;++l){
      if ((AP.BA[l][0])/3 == atomi ){
  	if ((AP.BA[l][1])/3 == atomj)
  	  return NO;
      }
      else if ((AP.BA[l][1])/3 == atomi){
  	if ((AP.BA[l][0])/3 == atomj)
  	  return NO;
      }
      else if ((AP.BA[l][0])/3 == atomj){
  	if ((AP.BA[l][1])/3 == atomi)
  	  return NO;
      }
      else if ((AP.BA[l][1])/3 == atomj){
  	if ((AP.BA[l][0])/3 == atomi)
  	  return NO;
      }
    }
  
    for(l=0;l<AP.NTHETH;++l){
      if ((AP.TH[l][0])/3 == atomi ){
  	if ((AP.TH[l][2])/3 == atomj)
  	  return NO;
      }
      else if ((AP.TH[l][2])/3 == atomi){
  	if ((AP.TH[l][0])/3 == atomj)
  	  return NO;
      }
      else if ((AP.TH[l][0])/3 == atomj){
  	if ((AP.TH[l][2])/3 == atomi)
  	  return NO;
      }
      else if ((AP.TH[l][2])/3 == atomj){
  	if ((AP.TH[l][0])/3 == atomi)
  	  return NO;
      }
    }
  
    for(l=0;l<AP.MTHETA;++l){
      if ((AP.TA[l][0])/3 == atomi ) {
  	if ((AP.TA[l][2])/3 == atomj)
  	  return NO;
      }
      //      else if ((AP.TH[l][2])/3 == atomi){
      else if ((AP.TA[l][2])/3 == atomi){ //1011
  	if ((AP.TA[l][0])/3 == atomj)
  	  return NO;
      }
      else if ((AP.TA[l][0])/3 == atomj){
  	if ((AP.TA[l][2])/3 == atomi)
  	  return NO;
      }
      else if ((AP.TA[l][2])/3 == atomj){
  	if ((AP.TA[l][0])/3 == atomi)
  	  return NO;
      }
    }
  }
  
  return YES;
}

double ffL_calcffandforce(double *crd, int numatom,struct potential *ene,struct force *f) {
  int i;
  int numnb,num14;
  double *n_d;

  numnb=(*ene).parm.numnb;
  num14=(*ene).parm.num14;

  ffL_calcFFNB((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e,(*ene).p_LJ,(*f).f_e,(*f).f_LJ,numnb,(*ene).parm.indexnb,numatom,crd,1,1);
  ffL_calcFFNB/*_14*/((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,(*ene).p_e_14,(*ene).p_LJ_14,(*f).f_e_14,(*f).f_LJ_14,num14,(*ene).parm.index14,numatom,crd,1,1);

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

  for (i=0;i<numnb;++i) {
    (*ene).p_t+=(*ene).p_e[i]+(*ene).p_LJ[i];
    (*ene).p_e_t+=(*ene).p_e[i];
    (*ene).p_LJ_t+=(*ene).p_LJ[i];
  }
  for (i=0;i<num14;++i) {
    //    (*ene).p_t+=1.0/1.2*(*ene).p_e_14[i]+0.5*(*ene).p_LJ_14[i]; // 1011
    //    (*ene).p_e_14_t+=/*1.0/1.2**/(*ene).p_e_14[i];                  // 1011
    //    (*ene).p_LJ_14_t+=/*0.5**/(*ene).p_LJ_14[i];                    // 1011
    (*ene).p_e_14_t+=(*ene).p_e_14[i];             
    (*ene).p_LJ_14_t+=(*ene).p_LJ_14[i];               
  }
  (*ene).p_e_14_t=1.0/1.2*(*ene).p_e_14_t;
  (*ene).p_LJ_14_t=0.5*(*ene).p_LJ_14_t;
  (*ene).p_t=(*ene).p_e_14_t+(*ene).p_LJ_14_t;

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
    (*f).f_t[i]+=(*f).f_e[i]+(*f).f_LJ[i]+(*f).f_e_14[i]+(*f).f_LJ_14[i]+(*f).f_d[i]+(*f).f_a[i]+(*f).f_b[i];
  
  return (*ene).p_t;
}

int ffL_calcFFNB(double *ele, double *ALJ, double *BLJ,
		double *p_e,double *p_LJ,
		double *f_e,double *f_LJ,
		int numnb, int *indexnb,
		int num_atom,
		double *cord,
		int flagp, int flagf) {
  int i,j;
  int num_a_prot,NUM_A_PROT;
  double potedummy;
  double f[3];
  double len,len2,len6;
  double vec[3],fdummy[3];
  double p12,p6;

  if (flagf != 0 && flagf != 1 && flagf !=2 ) {
    printf("error !\n");exit(1);
  }
  if (flagp != 0 && flagp != 1 && flagp !=2 ) {
    printf("error !\n");exit(1);
  }

  if (flagp==1 ) {
    for(i=0;i<num_atom;++i) {
      p_e[i]=0.0;p_LJ[i]=0.0;
      if (flagf==1) {
	for(j=0;j<3;++j) {
	  f_LJ[i*3+j]=0.0;
	  f_e[i*3+j]=0.0;
	}
      }
    }
  }

  if (flagp==2 ) {
    for(i=0;i<numnb;++i) {
      p_e[i]=0.0;
      p_LJ[i]=0.0; 
      if (flagf==2) {
	for(j=0;j<3;++j) {
	  f_LJ[i*3+j]=0.0;
	  f_e[i*3+j]=0.0;
	}
      }
    }
  }

  for(i=0;i<numnb;++i){
    num_a_prot=indexnb[i*2];
    NUM_A_PROT=indexnb[i*2+1];
    len2 = 0.0;
    for(j=0;j<3;++j){
      vec[j] = cord[NUM_A_PROT*3+j]-cord[num_a_prot*3+j];
      len2 += vec[j]*vec[j];
    }
    len = sqrt(len2);
    potedummy=ele[num_a_prot]*ele[NUM_A_PROT]/(len);
    if (flagp==1) {
      p_e[num_a_prot] += potedummy;
      p_e[NUM_A_PROT] += potedummy;
    }
    else if (flagp==2) p_e[i] = potedummy;

    if (flagf==1) {
      for(j=0; j<3; ++j) {
	fdummy[j] = -potedummy*UNIT*vec[j]/len2;
	f_e[num_a_prot*3+j] += fdummy[j];
	f_e[NUM_A_PROT*3+j] -= fdummy[j];
      }
    }
    else if (flagf==2) {
      for(j=0; j<3; ++j) {
	fdummy[j] = -potedummy*UNIT*vec[j]/len2;
	f_e[i*3+j] += fdummy[j];
      }
    }
    len6=len2;
    for (j=0;j<2;++j)  len6 = len6*len2;
    p12 = ALJ[num_a_prot*num_atom+NUM_A_PROT]/(len6*len6);
    p6  = BLJ[num_a_prot*num_atom+NUM_A_PROT]/len6;
    if (flagp==1) {
      p_LJ[num_a_prot] += p12-p6;
      p_LJ[NUM_A_PROT] += p12-p6;
    }
    else if (flagp==2) p_LJ[i] = p12-p6;

    if (flagf==1) {
      for(j=0;j<3;++j) {
	fdummy[j]=-6.0*(2.0*p12-p6)/(len2)*vec[j]*UNIT;
	f_LJ[num_a_prot*3+j] +=fdummy[j];
	f_LJ[NUM_A_PROT*3+j] -= fdummy[j];
      }
    }
    else if (flagf==2) {
      for(j=0;j<3;++j) {
	fdummy[j]=-6.0*(2.0*p12-p6)/(len2)*vec[j]*UNIT;
	f_LJ[i*3+j] +=fdummy[j];
      }
    }
  }

  return 0;
}

int ffL_calcDIHE(double *p_d,
		double *n_d,
		double *cord,
		int flagp, int flagf, int flaginp) {
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
      p_d[i] = AP.PK[dtype]*(1.0+cos(AP.PN[dtype]*dihedang-AP.PHASE[dtype]));
      if (flagf==1)
  	n_d[i] = -AP.PK[dtype]*4.18407*100.0*(sin(AP.PN[dtype]*dihedang-AP.PHASE[dtype])*AP.PN[dtype]);
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
      p_d[i+AP.NPHIH] = AP.PK[dtype]*(1.0+cos(AP.PN[dtype]*dihedang-AP.PHASE[dtype]));
      if (flagf==1)
  	n_d[i+AP.NPHIH] = -AP.PK[dtype]*4.18407*100.0*(sin(AP.PN[dtype]*dihedang-AP.PHASE[dtype])*AP.PN[dtype]);
    }
  }
}

int ffL_calcDIHE_force_Cartesian(double *f_d,double *cord) {
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

    if (AP.PN[dtype]==1) fa=-AP.PK[dtype]*AP.PN[dtype]*cos(AP.PHASE[dtype]);
    else if (AP.PN[dtype]==2) fa=-AP.PK[dtype]*AP.PN[dtype]*2.0*cosdih*cos(AP.PHASE[dtype]);
    else if (AP.PN[dtype]==3) fa=-AP.PK[dtype]*AP.PN[dtype]*(-4.0*sindih*sindih+3.0)*cos(AP.PHASE[dtype]);
    else if (AP.PN[dtype]==4) fa=-AP.PK[dtype]*AP.PN[dtype]*4.0*(cosdih*(2.0*cosdih*cosdih-1.0))*cos(AP.PHASE[dtype]);
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

int ffL_calcANGLE(double *p_a,double *cord){
  int i,j,k,l;
  int type;
  
  double atom[3][3];
  double ang;
  
  for (i=0;i<AP.NTHETH+AP.MTHETA;++i)
    p_a[i] = 0.0;
  
  for (i=0;i<AP.NTHETH;++i) {
    type = AP.TH[i][3]-1;
     for (j=0;j<3;++j)
      for (k=0;k<3;++k)
	atom[j][k]=cord[abs(AP.TH[i][j])+k];
  
    ang = pick_angle(atom[0],atom[1],atom[2],0,0.0);
    p_a[i] = AP.TK[type]*(ang-AP.TEQ[type])*(ang-AP.TEQ[type]);
  }
  
  for (i=0;i<AP.MTHETA;++i) {
    type = AP.TA[i][3]-1;
     for (j=0;j<3;++j)
      for (k=0;k<3;++k)
	atom[j][k]=cord[abs(AP.TA[i][j])+k];
  
    ang = pick_angle(atom[0],atom[1],atom[2],0,0.0);
    p_a[i+AP.NTHETH] = AP.TK[type]*(ang-AP.TEQ[type])*(ang-AP.TEQ[type]);
  }

  return 0;
}

int ffL_calcANGLE_force_Cartesian(double *f_a,double *cord){
  int i,j,k,l;
  int numatom;
  int type;
  double kang,ang_eq;

  double atom[3][3];
  double *f_temp;

  numatom=AP.NATOM;
  for (i=0;i<numatom*3;++i) f_a[i] = 0.0;
  f_temp=(double *)gcemalloc(sizeof(double)*3*3);

  for (i=0;i<AP.NTHETH;++i) {
    for (j=0;j<3;++j)  for (k=0;k<3;++k) atom[j][k]=cord[abs(AP.TH[i][j])+k];
    type = AP.TH[i][3]-1;
    kang = AP.TK[type];
    ang_eq = AP.TEQ[type];
    calcANGKE_force(atom[0],atom[1],atom[2],kang,ang_eq,f_temp);
    for (j=0;j<3;++j) {
      f_a[abs(AP.TH[i][0])+j] += f_temp[j];
      f_a[abs(AP.TH[i][1])+j] += f_temp[3+j];
      f_a[abs(AP.TH[i][2])+j] += f_temp[6+j];
    }
  }

  for (i=0;i<AP.MTHETA;++i) {
    for (j=0;j<3;++j)  for (k=0;k<3;++k) atom[j][k]=cord[abs(AP.TA[i][j])+k];
    type = AP.TA[i][3]-1;
    kang = AP.TK[type];
    ang_eq = AP.TEQ[type];
    calcANGKE_force(atom[0],atom[1],atom[2],kang,ang_eq,f_temp);
    for (j=0;j<3;++j) {
      f_a[abs(AP.TA[i][0])+j] += f_temp[j];
      f_a[abs(AP.TA[i][1])+j] += f_temp[3+j];
      f_a[abs(AP.TA[i][2])+j] += f_temp[6+j];
    }
  }

  return 0;
}

double calcANGKE_force(double atomi[3],double atomj[3],double atomk[3],double kang,double ang_eq,double *f) {
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
    /* f1 = -kang*(angijk-ang_eq)/(lenij*sin(angijk))*(vkj[j]/lenkj-cosijk*vij[j]/lenij)*UNIT; */
    /* f2 = -kang*(angijk-ang_eq)/(lenkj*sin(angijk))*(vij[j]/lenij-cosijk*vkj[j]/lenkj)*UNIT; */
    /*******************************************************************************************/
    f1 = -2.0*kang*(angijk-ang_eq)/(lenij*sin(angijk))*(vkj[j]/lenkj-cosijk*vij[j]/lenij)*UNIT;
    f2 = -2.0*kang*(angijk-ang_eq)/(lenkj*sin(angijk))*(vij[j]/lenij-cosijk*vkj[j]/lenkj)*UNIT;

    f[j] = f1;
    f[6+j] = f2;
    f[3+j] = -f1-f2;
  }
}

int ffL_calcBOND(double *p_b,double *cord){
  int i,j,k;
  int type;
  double len;
  double atom[2][3];

  for (i=0;i<AP.NBONH+AP.MBONA;++i)
    p_b[i] = 0.0;
  
  for (i=0;i<AP.NBONH;++i) {
    type = AP.BH[i][2]-1;
    for (j=0;j<2;++j)
      for (k=0;k<3;++k)
	atom[j][k]=cord[abs(AP.BH[i][j])+k];
  
    len = pick_bond_leng(atom[0],atom[1]);
    p_b[i] = AP.RK[type]*(len-AP.REQ[type])*(len-AP.REQ[type]);
  }

  for (i=0;i<AP.MBONA;++i) {
    type = AP.BA[i][2]-1;
    for (j=0;j<2;++j)
      for (k=0;k<3;++k)
	atom[j][k]=cord[abs(AP.BA[i][j])+k];  

    len = pick_bond_leng(atom[0],atom[1]);
    p_b[i+AP.NBONH] = AP.RK[type]*(len-AP.REQ[type])*(len-AP.REQ[type]);
  }

  return 0;
}

int ffL_calcBOND_force_Cartesian(double *f_b,double *cord){
  int i,j,k;
  int numatom;
  int type;
  double f;
  double lenij;
  double atom[2][3];

  numatom=AP.NATOM;
  for (i=0;i<numatom*3;++i) f_b[i] = 0.0;
  
  for (i=0;i<AP.NBONH;++i) {
    type = AP.BH[i][2]-1;
    for (j=0;j<2;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(AP.BH[i][j])+k];
  
    lenij = len(atom[0],atom[1]);
    for (j=0;j<3;++j) {
      f = -2.0*AP.RK[type]*(lenij-AP.REQ[type])*(atom[1][j]-atom[0][j])/lenij*UNIT;
      f_b[abs(AP.BH[i][0])+j] += f;
      f_b[abs(AP.BH[i][1])+j] += -f;
    }
  }

  for (i=0;i<AP.MBONA;++i) {
    type = AP.BA[i][2]-1;
    for (j=0;j<2;++j) for (k=0;k<3;++k) atom[j][k]=cord[abs(AP.BA[i][j])+k];
  
    lenij = len(atom[0],atom[1]);
    for (j=0;j<3;++j) {
      f = -2.0*AP.RK[type]*(lenij-AP.REQ[type])*(atom[1][j]-atom[0][j])/lenij*UNIT;
      f_b[abs(AP.BA[i][0])+j] += f;
      f_b[abs(AP.BA[i][1])+j] += -f;
    }
  }

  return 0;
}

///////////////////////////////////////////////////////////////////
