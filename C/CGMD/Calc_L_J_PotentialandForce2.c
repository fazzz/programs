#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "MD.h"
#include "force.h"
#include "physics.h"
#include "math.h"

//  non bonding 
double Calc_L_J_PotentialandForce2(double *ele, double *ALJ, double *BLJ,
				   double *p_e, double *p_1_4_e, 
				   double *p_LJ,double *p_1_4_LJ,
				   double *f_e, double *f_1_4_e, 
				   double *f_LJ,double *f_1_4_LJ,
				   int numnb, int *indexnb,
				   int num14, int *index14,
				   double *cord,double *eig,double *eig_14)
{
  int i,j,ii; 
  int alpha;
  
  int nNumClut, nNumClut2;
  
  double potedummy;
  double f[3];
  double d_i_j,l,len,len2,len6;
  double vec[3],fdummy[3];
  double p12,p6;

  int tDOF=prot.DOF;
  int tnum_atom=prot.num_atom;
  int tnum_atom_clust;
  int tnum_atom_clust2;

  for(i = 0;i < tnum_atom; ++i){
      p_e[i] = 0.0;
      p_1_4_e[i] = 0.0;
      p_LJ[i] = 0.0;
      p_1_4_LJ[i] = 0.0;
      for(alpha=0; alpha<3; ++alpha) {
	f_LJ[i*3+alpha] = 0.0;
	f_1_4_LJ[i*3+alpha] = 0.0;
	f_e[i*3+alpha] = 0.0;
	f_1_4_e[i*3+alpha] = 0.0;
      }
  }

  if (nbstopflag!=3) {
    for(i=0;i<numnb;++i){
      num_a_prot=indexnb[i*2];
      NUM_A_PROT=indexnb[i*2+1];
      len2 = 0.0;
      for(alpha=0;alpha<3;++alpha){
	vec[alpha] = cord[NUM_A_PROT*3+alpha]-cord[num_a_prot*3+alpha];
	len2 += vec[alpha]*vec[alpha];
      }
      len = sqrt(len2);
      potedummy=ele[num_a_prot]*ele[NUM_A_PROT]/(len);
      if ( PEPCAACCMODE == ON ) {
	potedummy+=potedummy*eig[i*2]*fact;
      }
      else if ( PEPCAACCMODE == 2) {
	potedummy=potedummy*eig[i*2]*fact;
      }
      else if ( PEPCAACCMODE == 3) {
	potedummy=potedummy*fact+potedummy*eig[i*2]*fact;
      }
      if (MASSIVEOUT==ON && (nNumStep%out_put_steps_thomo)==0){
	//	fprintf(data,"%e  ",potedummy);
      }
      p_e[num_a_prot] += potedummy;
      p_e[NUM_A_PROT] += potedummy;
      for(alpha=0; alpha<3; ++alpha) {
	fdummy[alpha] = -potedummy*4.184070*100.0*vec[alpha]/len2;
	f_e[num_a_prot*3+alpha] += fdummy[alpha];
	f_e[NUM_A_PROT*3+alpha] -= fdummy[alpha];
      }
      len6=len2;
      for (j=0;j<2;++j) {
	len6 = len6*len2;
      }

      p12 = ALJ[num_a_prot*tnum_atom+NUM_A_PROT]/(len6*len6);
      p6  = BLJ[num_a_prot*tnum_atom+NUM_A_PROT]/len6;
      if ( PEPCAACCMODE == ON ) {
	p12+=p12*eig[i*2+1]*fact;
	p6 +=p6*eig[i*2+1]*fact;
      }
      else if ( PEPCAACCMODE == 2 ) {
	p12=p12*eig[i*2+1]*fact;
	p6 =p6*eig[i*2+1]*fact;
      }
      else if ( PEPCAACCMODE == 3) {
	p12=p12*eig[i*2+1]*fact;
	p6 =p6*eig[i*2+1]*fact;
	potedummy=potedummy*fact+potedummy*eig[i*2]*fact;
      }
      p_LJ[num_a_prot] += p12-p6;
      p_LJ[NUM_A_PROT] += p12-p6;
      if (MASSIVEOUT==ON && (nNumStep % out_put_steps_thomo)==0){
	//	fprintf(data2,"%e  ",p12-p6);
      }
      for(alpha=0;alpha<3;++alpha) {
	fdummy[alpha]=-6.0*(2.0*p12-p6)/(len2)*vec[alpha]*4.184070*100.0;
	f_LJ[num_a_prot*3+alpha] +=fdummy[alpha];
	f_LJ[NUM_A_PROT*3+alpha] -= fdummy[alpha];
      }
    }
  }
 
if (MASSIVEOUT==ON && (nNumStep%out_put_steps_thomo)==0){
  //  fprintf(data,"\n  ");
  //    fprintf(data2,"\n  ");
  //    fclose(data);
  //    fclose(data2);
  }
if (nbstopflag!=2) {
    for(i = 0;i < num14; ++i){
      num_a_prot=index14[i*2];
      NUM_A_PROT=index14[i*2+1];
      len2 = 0.0;
      for(alpha=0;alpha<3;++alpha){
	vec[alpha] = cord[NUM_A_PROT*3+alpha]-cord[num_a_prot*3+alpha];
	len2 += vec[alpha]*vec[alpha];
      }
      len = sqrt(len2);
      potedummy=ele[num_a_prot]*ele[NUM_A_PROT]/(len);
      if ( PEPCAACCMODE == ON ) {
	potedummy+=potedummy*eig_14[i*2]*fact;
      }
      else if ( PEPCAACCMODE == 2 ) {
	potedummy=potedummy*eig_14[i*2]*fact;
      }
      if (MASSIVEOUT==ON && (nNumStep%out_put_steps_thomo)==0){
	//	fprintf(data3,"%e  ",1.0/1.2*potedummy);
      }
      p_1_4_e[num_a_prot] += 1.0/1.2*potedummy;
      p_1_4_e[NUM_A_PROT] += 1.0/1.2*potedummy;
      for(alpha=0; alpha<3; ++alpha) {
	fdummy[alpha] = -1.0/1.2*potedummy*4.184070*100.0*vec[alpha]/len2;
	f_1_4_e[num_a_prot*3+alpha] += fdummy[alpha];
	f_1_4_e[NUM_A_PROT*3+alpha] -= fdummy[alpha];
      }
      len6=len2;
      for (j=0;j<2;++j) {
	len6 = len6*len2;
      }
      p12 = ALJ[num_a_prot*tnum_atom+NUM_A_PROT]/(len6*len6);
      p6  = BLJ[num_a_prot*tnum_atom+NUM_A_PROT]/len6;
      if ( PEPCAACCMODE == ON ) {
	p12+=p12*eig_14[i*2+1]*fact;
	p6 +=p6*eig_14[i*2+1]*fact;
      }
      if ( PEPCAACCMODE == 2 ) {
	p12=p12*eig_14[i*2+1]*fact;
	p6 =p6*eig_14[i*2+1]*fact;
      }
      if (MASSIVEOUT==ON && (nNumStep % out_put_steps_thomo)==0){
	//	fprintf(data4,"%e  ",0.5*(p12-p6));
      }
    
      p_1_4_LJ[num_a_prot] += 0.5*(p12-p6);
      p_1_4_LJ[NUM_A_PROT] += 0.5*(p12-p6);
    
      for(alpha=0;alpha<3;++alpha) {
	fdummy[alpha] = -3.0*(2.0*p12-p6)/(len*len)*vec[alpha]*4.184070*100.0;
	f_1_4_LJ[num_a_prot*3+alpha] += fdummy[alpha];
	f_1_4_LJ[NUM_A_PROT*3+alpha] -= fdummy[alpha];
      }
    }
  }
if (MASSIVEOUT==ON && (nNumStep%out_put_steps_thomo)==0){
  //  fprintf(data3,"\n  ");
  //  fprintf(data4,"\n  ");
  //  fclose(data3);
  //  fclose(data4);
 }

  num_a_prot=0;
  for(nNumClut = 0;nNumClut < tDOF; ++nNumClut) {
    tnum_atom_clust = clust[nNumClut].num_atom_clust;
    for(i=0;i < tnum_atom_clust;++i) {
      potential_pro.p_elesta[num_a_prot] = p_e[num_a_prot];
      if (nbstopflag!=2)
	potential_pro.p_1_4_elesta[num_a_prot] = p_1_4_e[num_a_prot];
      if (nbstopflag!=3)
	potential_pro.p_L_J[num_a_prot] = p_LJ[num_a_prot];
      if (nbstopflag!=2)
	potential_pro.p_1_4_L_J[num_a_prot] = p_1_4_LJ[num_a_prot];
      for(alpha=0; alpha<3; ++alpha) {
	if (nbstopflag!=3)
  	clust[nNumClut].f_c.f_L_J[i][alpha] = f_LJ[num_a_prot*3+alpha];
	if (nbstopflag!=2)
	  clust[nNumClut].f_c.f_1_4_L_J[i][alpha] = f_1_4_LJ[num_a_prot*3+alpha];
	if (nbstopflag!=3)
	  clust[nNumClut].f_c.f_elesta[i][alpha] = f_e[num_a_prot*3+alpha];
	if (nbstopflag!=2)
	  clust[nNumClut].f_c.f_1_4_elesta[i][alpha] = f_1_4_e[num_a_prot*3+alpha];
      }
      ++num_a_prot;
    }
  }

}

void set_non_bonding_parameters(double *ele, double *ALJ, double *BLJ){
  int i,j,index;
  int nNumClut,num_a_prot;

  num_a_prot=0;
  for (nNumClut=0;nNumClut<prot.DOF;++nNumClut) {
    for (i=0;i<clust[nNumClut].num_atom_clust;++i) {
      ele[num_a_prot]= clust[nNumClut].f_p_clust.e_f[i];
      ++num_a_prot;
    }
  }

  for(i = 0;i < prot.num_atom; ++i) {
    for(j = 0;j < prot.num_atom; ++j) {
      index=prot.L_J_parm.atomtypeIndex[prot.L_J_parm.atomtype[i]-1][prot.L_J_parm.atomtype[j]-1]-1;
      ALJ[i*prot.num_atom+j] = prot.L_J_parm.A[index];
      BLJ[i*prot.num_atom+j] = prot.L_J_parm.B[index];
    }
  }
}

void set_non_bonding_index(void) {
  int i,j;
  int tDOF,tnum_atom_clust,tnum_atom_clust2,tnum_atom;
  int num_a_prot,NUM_A_PROT;
  int nNumClut,nNumClut2;
  int numnbl=0;
  int num14l=0;
  FILE *log;

  gnumnb=0;
  gnum14=0;

  tnum_atom=prot.num_atom;
  tDOF=prot.DOF;
  
  num_a_prot = 0;
  for(nNumClut = 0;nNumClut < tDOF; ++nNumClut){
    tnum_atom_clust = clust[nNumClut].num_atom_clust;
    for(i=0;i < tnum_atom_clust;++i) {
      NUM_A_PROT = num_a_prot+1;
      if (i==tnum_atom_clust-1) {
	nNumClut2=nNumClut+1;
      }
      else {
	nNumClut2=nNumClut;
      }
    
      for(;nNumClut2 < tDOF; ++nNumClut2){
	if (nNumClut2>nNumClut) {
	  j=0;
	}
	else {
	  j = i + 1;
	}
	tnum_atom_clust2 = clust[nNumClut2].num_atom_clust;
	for(;j < tnum_atom_clust2;++j){
	  if (which_calc_non_bonding_pote(nNumClut, i, NUM_A_PROT) == 1) {
	    ++gnumnb;
	  }
	  else if(which_calc_1_4_non_bonding_pote(nNumClut,i, NUM_A_PROT)==1){
	    ++gnum14;
	  }
	  ++NUM_A_PROT;
	}
      }
      ++num_a_prot;
    }
  }

  /**************************/
  /* printf("%d\n",gnumnb); */
  /* printf("%d\n",gnum14); */
  /**************************/


  gindexnb=(int *)calloc(gnumnb*2,sizeof(int));
  gindex14=(int *)calloc(gnum14*2,sizeof(int));

  num_a_prot = 0;
  for(nNumClut = 0;nNumClut < tDOF; ++nNumClut){
    tnum_atom_clust = clust[nNumClut].num_atom_clust;
    for(i=0;i < tnum_atom_clust;++i) {
      NUM_A_PROT = num_a_prot+1;
      if (i==tnum_atom_clust-1) {
	nNumClut2=nNumClut+1;
      }
      else {
	nNumClut2=nNumClut;
      }
    
      for(;nNumClut2 < tDOF; ++nNumClut2){
	if (nNumClut2>nNumClut) {
	  j=0;
	}
	else {
	  j = i + 1;
	}
	tnum_atom_clust2 = clust[nNumClut2].num_atom_clust;
	for(;j < tnum_atom_clust2;++j){
	  if (which_calc_non_bonding_pote(nNumClut, i, NUM_A_PROT) == 1) {
	    gindexnb[numnbl*2]=num_a_prot;
	    gindexnb[numnbl*2+1]=NUM_A_PROT;
	    ++numnbl;
	  }
	  else if(which_calc_1_4_non_bonding_pote(nNumClut,i, NUM_A_PROT)==1){
	    gindex14[num14l*2]=num_a_prot;
	    gindex14[num14l*2+1]=NUM_A_PROT;
	    ++num14l;
	  }
	  ++NUM_A_PROT;
	}
      }
      ++num_a_prot;
    }
  }

  //  if (MASSIVEOUT==ON){
    if((log=fopen("log.txt","w"))==NULL)
      printf("error\n");
    for(i=0;i<numnbl;++i){
      fprintf(log,"%d - %d \n",gindexnb[i*2]+1,gindexnb[i*2+1]+1);
    }
    fprintf(log,"\n");
    for(i=0;i<num14l;++i){
      fprintf(log,"%d - %d \n",gindex14[i*2]+1,gindex14[i*2+1]+1);
    }

    fprintf(log,"\n  ");
    fclose(log);
    //  }


}
