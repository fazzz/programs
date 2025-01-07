
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GOLMAA.h"
#include "GOLMAA_set.h"
#include "NC.h"
#include "PTL.h"
#include "TOPO.h"
#include "MB.h"

double GOLMAAff_calcff(double *crd, int numatom,struct potential_GOLMAA *ene,
		       int flagd, int flagnc, int flagnn,int **nb_matrix) {
  int i,j;

  if (flagd==ON) GOLMAApote_calcDIHE((*ene).p_d,crd,(*ene).DEQ,(*ene).FC_dihed);

  GOLMAApoteforc_calcNatAttandRepul((*ene).p_natatt,(*ene).p_repul,(*ene).f_natatt,(*ene).f_repul,crd,(*ene).ALJ_natatt,(*ene).BLJ_natatt,(*ene).ep_natatt,(*ene).ALJ_repul,(*ene).num_natatt,numatom,(*ene).ncmap,nb_matrix);

  (*ene).p_t=0.0;
  (*ene).p_b_t=0.0;
  (*ene).p_a_t=0.0;
  (*ene).p_d_t=0.0;
  (*ene).p_natatt_t=0.0;
  (*ene).p_repul_t=0.0;

  if (flagd==ON) {
    for (i=0;i<AP.MPHIA;++i) {
      (*ene).p_t+=(*ene).p_d[i];
      (*ene).p_d_t+=(*ene).p_d[i];
    }
  }

  if (flagnc==ON) {
    for (i=0;i<numatom;++i) {
      (*ene).p_t+=0.5*(*ene).p_natatt[i];
      (*ene).p_natatt_t+=0.5*(*ene).p_natatt[i];
    }
  }

  if (flagnn==ON) {
    for (i=0;i<numatom;++i) {
      (*ene).p_t+=0.5*(*ene).p_repul[i];
      (*ene).p_repul_t+=0.5*(*ene).p_repul[i];
    }
  }

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) (*ene).f_t[i][j]=0.0;
  //  if (flagd==ON) for (i=0;i<numatom;++i) for (j=0;j<3;++j) (*ene).f_t[i][j]+=(*ene).f_d[i][j];
  if (flagnc==ON) for (i=0;i<numatom;++i) for (j=0;j<3;++j) (*ene).f_t[i][j]+=(*ene).f_natatt[i][j];
  if (flagnn==ON) for (i=0;i<numatom;++i) for (j=0;j<3;++j) (*ene).f_t[i][j]+=(*ene).f_repul[i][j];
 
  return (*ene).p_t;
}

void GOLMAApote_calcDIHE(double *p_d,double *crd,double *DEQ,double FC_dihed){
  int i,j,k,l;
  double atom[4][3];
  double cosdih,sindih;
  double dihedang,dang;
  double pi;

  pi=acos(-1.0);
  for (i=0;i<AP.MPHIA;++i) {
    for (j=0;j<4;++j) for (k=0;k<3;++k)	atom[j][k]=crd[abs(AP.PA[i][j])+k];
  
    dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
    if (dihedang>pi) dihedang-=2.0*pi;
    else if (dihedang<-1.0*pi) dihedang+=2.0*pi;

    if ((dang=dihedang-DEQ[i])>pi) dang-=2.0*pi;
    else if ((dang=dihedang-DEQ[i])<-1.0*pi) dang+=2.0*pi;
    p_d[i] = FC_dihed*((1.0-cos(dang))+0.5*(1.0-cos(3.0*(dang))));
  }

}

int GOLMAApote_calcNatAtt(double *p_natatt,double *crd,double *ALJ_natatt,double *BLJ_natatt,int num_natatt,int numatom,int *index_natatt) {
  int i,j;
  int num_a_prot,NUM_A_PROT;
  double vec[3];
  double len,len2,len6,len12;
  double p12,p6;

  for(i=0;i<numatom;++i) p_natatt[i]=0.0;

  for(i=0;i<num_natatt;++i){
    num_a_prot=index_natatt[i*2];
    NUM_A_PROT=index_natatt[i*2+1];
    len2 = 0.0;
    for(j=0;j<3;++j){
      vec[j] = crd[NUM_A_PROT*3+j]-crd[num_a_prot*3+j];
      len2 += vec[j]*vec[j];
    }
    len = sqrt(len2);
    len6=len2;
    len12=len2;
    for (j=0;j<2;++j)  len6 = len6*len2;
    for (j=0;j<5;++j)  len12 = len12*len2;
    p12 = ALJ_natatt[i]/len12;
    p6 = 2.0*BLJ_natatt[i]/len6;
    p_natatt[num_a_prot] += p12-p6;
    p_natatt[NUM_A_PROT] += p12-p6;
  }

  return 0;
}

int GOLMAApote_calcRepul(double *p_repul,double *crd,double ALJ_repul,int numatom) {
  int i,j,k;
  double vec[3];
  double len,len2,len12;
  double p12;

  for(i=0;i<numatom;++i) p_repul[i]=0.0;

  for (i=0;i<numatom;++i) {
    for (j=i;j<numatom;++j) {
      if (j>i+3) {
	if (strncmp(AP.IGRAPH[i],"H",1)!=0 && strncmp(AP.IGRAPH[j],"H",1)!=0) {
	  len2 = 0.0;
	  for(k=0;k<3;++k){
	    vec[k] = crd[j*3+k]-crd[i*3+k];
	    len2 += vec[k]*vec[k];
	  }
	  len = sqrt(len2);
	  len12=len2;
	  for (k=0;k<5;++k)  len12 = len12*len2;
	  p12 = ALJ_repul/len12;
	  p_repul[i] += p12;
	  p_repul[j] += p12;
  	}
      }
    }
  }

  return 0;
}

int GOLMAAforc_calcNatAtt(double **f_natatt,double *crd,double *ALJ_natatt,double *BLJ_natatt,int num_natatt,int numatom,int *index_natatt) {
  int i,j;
  int num_a_prot,NUM_A_PROT;
  double vec[3];
  double len,len2,len6,len12;
  double p12,p6;

  for(i=0;i<numatom;++i) for(j=0;j<3;++j) f_natatt[i][j]=0.0;

  for(i=0;i<num_natatt;++i){
    num_a_prot=index_natatt[i*2];
    NUM_A_PROT=index_natatt[i*2+1];
    len2 = 0.0;
    for(j=0;j<3;++j){
      vec[j] = crd[NUM_A_PROT*3+j]-crd[num_a_prot*3+j];
      len2 += vec[j]*vec[j];
    }
    len = sqrt(len2);
    len6=len2;
    len12=len2;
    for (j=0;j<2;++j)  len6 = len6*len2;
    for (j=0;j<5;++j)  len12 = len12*len2;
    p12 = ALJ_natatt[i]/len12;
    p6 = 2.0*BLJ_natatt[i]/len6;
    if (1.0<p12-p6+1.0) {
      for (j=0;j<3;++j) {
	f_natatt[num_a_prot][j] += -(12.0*p12-6.*p6)/(len2)*vec[j]*UNIT;
	f_natatt[NUM_A_PROT][j] += (12.0*p12-6.*p6)/(len2)*vec[j]*UNIT;
      }
    }
  }

  return 0;
}

int GOLMAAforc_calcRepul(double **f_repul,double *crd,double ALJ_repul,int numatom) {
  int i,j,k;
  int ii,jj;
  double vec[3];
  double len,len2,len12;
  double p12;

  for(i=0;i<numatom;++i) for(j=0;j<3;++j) f_repul[i][j]=0.0;

  ii=0;
  for (i=0;i<numatom;++i) {
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
      ++ii;
    }
    jj=ii;
    for (j=i;j<numatom;++j) {
      if ( strncmp(AP.IGRAPH[j],"CA",2)==0 ) {
  	++jj;
      }
      if (strncmp(AP.IGRAPH[i],"CA",2)==0 && strncmp(AP.IGRAPH[j],"CA",2)==0) {
  	if (jj>ii+3) {
	  len2 = 0.0;
	  for(k=0;k<3;++k){
	    vec[k] = crd[j*3+k]-crd[i*3+k];
	    len2 += vec[k]*vec[k];
	  }
	  len = sqrt(len2);
	  len12=len2;
	  for (k=0;k<5;++k)  len12 = len12*len2;
	  p12 = ALJ_repul/len12;
	  for (k=0;k<3;++k) {
	    f_repul[i][k] += -(12.0*p12)/(len2)*vec[k]*UNIT;
	    f_repul[j][k] += (12.0*p12)/(len2)*vec[k]*UNIT;
	  }
  	}
      }
    }
  }

  return 0;
}

int GOLMAApoteforc_calcNatAttandRepul(double *p_natatt,double *p_repul,double **f_natatt,double **f_repul,double *crd,double *ALJ_natatt,double *BLJ_natatt,double e_natatt,double ALJ_repul,int num_natatt,int numatom,int **ncmap,int **nb_matrix) {
  int i,j,k,na;
  int num_a_prot,NUM_A_PROT;
  double vec[3];
  double len,len2,len6,len12;
  double p12,p6;

  for(i=0;i<numatom;++i) {
    for(j=0;j<3;++j){
      f_natatt[i][j]=0.0;
      f_repul[i][j]=0.0;
    }
    p_natatt[i]=0.0;
    p_repul[i]=0.0;
  }

  na=0;
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (ncmap[i][j]==0) {
	len2 = 0.0;
	for(k=0;k<3;++k){
	  vec[k] = crd[i*3+k]-crd[j*3+k];
	  len2 += vec[k]*vec[k];
	}
	len = sqrt(len2);
	len6=len2;
	len12=len2;
	for (k=0;k<2;++k)  len6 = len6*len2;
	for (k=0;k<5;++k)  len12 = len12*len2;
	p12 = ALJ_natatt[na]/len12;
	p6 = 2.0*BLJ_natatt[na]/len6;
	p_natatt[j] += e_natatt*(p12-p6);
	p_natatt[i] += e_natatt*(p12-p6);
	for (k=0;k<3;++k) {
	  f_natatt[i][k] +=  e_natatt*(12.0*p12-6.0*p6)/(len2)*vec[k]*UNIT;
	  f_natatt[j][k] += -e_natatt*(12.0*p12-6.0*p6)/(len2)*vec[k]*UNIT;
	}
	++na;
      }
      else {
	if (nb_matrix[i][j]==-1 && strncmp(AP.IGRAPH[i],"H",1)!=0 && strncmp(AP.IGRAPH[j],"H",1)!=0) {

          /**********************************************************/
          /* for (k=0;k<2;++k) printf("%c",AP.IGRAPH[i][k]);	    */
	  /* printf("(%d)",i);					    */
	  /* printf("-");					    */
	  /* for (k=0;k<2;++k) printf("%c",AP.IGRAPH[j][k]);	    */
	  /* printf("(%d)",j);					    */
	  /* printf("\n ");					    */
          /**********************************************************/
	  /***********************************************************************/
          /* for (k=0;k<3;++k) printf("%c",AP.LABERES[PTL_resnum(i)][k]);	 */
	  /* printf("%d",PTL_resnum(i));					 */
	  /* printf("-");							 */
	  /* for (k=0;k<3;++k) printf("%c",AP.IGRAPH[PTL_resnum(j)][k]);	 */
	  /* printf("%d",PTL_resnum(j));					 */
	  /* printf("\n ");							 */
          /***********************************************************************/

	  len2 = 0.0;
	  for(k=0;k<3;++k){
	    vec[k] = crd[i*3+k]-crd[j*3+k];
	    len2 += vec[k]*vec[k];
	  }
	  len = sqrt(len2);
	  len12=len2;
	  for (k=0;k<5;++k)  len12 = len12*len2;
	  p12 = ALJ_repul/len12;
	  p_repul[i] += p12;
	  p_repul[j] += p12;
	  for (k=0;k<3;++k) {
	    f_repul[i][k] +=  (12.0*p12)/(len2)*vec[k]*UNIT;
	    f_repul[j][k] += -(12.0*p12)/(len2)*vec[k]*UNIT;
	  }
	}
      }
    }
  }

  return 0;
}

double GOLMAA_calcTorque(double *Q,double *crd,double *DEQ,double FC_dihed,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a) {
  int i,j,k,l,n;
  int dtype,flag;  
  double atom[4][3];
  double dihedang,dang;
  double pot_d=0.0;
  int *inpindex,inpnumA,*indexclut;
  double pi;

  pi=acos(-1.0);

  indexclut=(int *)gcemalloc(sizeof(int)*(AP.MPHIA));
  inpindex=GOLMAAL_make_inpindex(&inpnumA,indexclut,numclut,nNumClutOfParent,terminal_atom_a,origin_atom_a);

  for (i=0;i<AP.MPHIA;++i) {
    flag=ON;
    for (j=0;j<inpnumA;++j) {
      if (i == inpindex[j]) {
	flag=OFF;
	break;
      }
    }

    if (flag==ON) {
      dtype = AP.PA[i][4]-1;
      n = indexclut[i]-1;
      for (j=0;j<4;++j) for (k=0;k<3;++k) atom[j][k]=crd[abs(AP.PA[i][j])+k];
  
      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      if (dihedang<-1.0*pi) dihedang+=2.0*pi;
      if (dihedang> pi) dihedang-=2.0*pi;

      /***************************************************************************************************************************/
      /* printf("%d(%d) %d-%d-%d-%d  %lf %lf\n",										 */
      /* 	     i,n,abs(AP.PA[i][0])/3-1,abs(AP.PA[i][1])/3-1,abs(AP.PA[i][2])/3-1,abs(AP.PA[i][3])/3-1,dihedang,DEQ[i]);	 */
      /***************************************************************************************************************************/

      dang=dihedang-DEQ[i];
      if (dang<-1.0*pi) dang+=2.0*pi;
      if (dang> pi) dang-=2.0*pi;
      pot_d+= FC_dihed*((1.0-cos(dang))+0.5*(1.0-cos(3.0*dang)));
      Q[n] += /*-*/FC_dihed*UNIT*(sin(dang)+0.5*(sin(3.0*dang))*3.0);
    }
  }

  return pot_d;
}

int *GOLMAAL_make_inpindex(int *inpnumA,int *indexclut,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a) {
  int i,j,k,l,ll,p;
  int flag;
  int *atom_dihed_pair;
  int *inpindexA;
  int *inpindex;

  (*inpnumA)=0;
  
  atom_dihed_pair=(int *)gcemalloc(sizeof(int)*(AP.MPHIA)*6);
  inpindexA=(int *)gcemalloc(sizeof(int)*1);
  
  for (k=0;k<AP.MPHIA;++k) {
    // check for improper dihed
    flag=OFF;
    if (flag==OFF) {
      for (i=0;i<AP.NBONA;++i) {
  	if (   ((AP.BA[i][0] == abs(AP.PA[k][1]) && AP.BA[i][1] == abs(AP.PA[k][0])) || (AP.BA[i][0] == abs(AP.PA[k][0]) && AP.BA[i][1] == abs(AP.PA[k][1])))) {
  	  flag = ON;
  	  break;
  	}
      }
    }
    if (flag==ON) {
      flag=OFF;
      if (flag==OFF ) {
  	for (i=0;i<AP.NBONA;++i) {
  	  if (((AP.BA[i][0] == abs(AP.PA[k][2]) && AP.BA[i][1] == abs(AP.PA[k][3])) || (AP.BA[i][0] == abs(AP.PA[k][3]) && AP.BA[i][1] == abs(AP.PA[k][2]))) ) {
  	    flag = ON;
  	    break;
  	  }
  	}
      }
    }
    if (flag==OFF) {
      inpindexA[(*inpnumA)]=k;
      ++(*inpnumA);
      inpindexA=(int *)gcerealloc(inpindexA,sizeof(int)*(*inpnumA));
    }
  
    for (i=0;i<4;++i) {
      atom_dihed_pair[k/*][*/*6+i] = abs(AP.PA[k][i])/3+1;
    }
    atom_dihed_pair[k/*][*/*6+4] = AP.PA[k][4];
    l = 1;
    ll = 2;
    if (atom_dihed_pair[k/*][*/*6+1] > atom_dihed_pair[k/*][*/*6+2]) {
      l = 2;
      ll = 1;
    }
    for (i=0;i<numclut;++i){
      p = nNumClutOfParent[i]-1;
      if (atom_dihed_pair[k/*][*/*6+l]==terminal_atom_a[p] && atom_dihed_pair[k/*][*/*6+ll]==origin_atom_a[i]){  
  	atom_dihed_pair[k/*][*/*6+5] = i+1;
      }
    }
  }
  
  for (i=0;i<AP.MPHIA;++i) indexclut[i]=atom_dihed_pair[i*6+5];


  //  inpindex=(int *)gcerealloc(inpindex,sizeof(int)*((*inpnumH)+(*inpnumA)));
  inpindex=(int *)gcemalloc(sizeof(int)*(*inpnumA));
  for (i=0;i<(*inpnumA);++i) inpindex[i]=inpindexA[i];

  return inpindex;  
}

void GOLMAA_out_formated(FILE *outputfile,struct potential_GOLMAA e,double KE,double KEv,double PEv,int i,double dt) {
  fprintf(outputfile,"/***********************************************/\n");
  fprintf(outputfile,"steps            = %d  \n",i);
  fprintf(outputfile,"total time       = %10.3lf ps  \n",dt*(double)i);
  fprintf(outputfile,"toal_energy               = %10.8e kcal/mol  \n",e.p_t+KE+KEv+PEv);

  fprintf(outputfile,"kinetic_energy            = %10.8e kcal/mol  \n",KE);
  fprintf(outputfile,"kinetic_v_energy          = %10.8e kcal/mol  \n",KEv);
  fprintf(outputfile,"potential_v_energy        = %10.8e kcal/mol  \n",PEv);
  fprintf(outputfile,"potential_energy          = %10.8e kcal/mol  \n",e.p_t);

  fprintf(outputfile,"dihedral_angle_energy     = %10.8e kcal/mol  \n",e.p_d_t);
  fprintf(outputfile,"native_contact_energy     = %10.8e kcal/mol  \n",e.p_natatt_t);
  fprintf(outputfile,"non_native_repul_energy   = %10.8e kcal/mol  \n",e.p_repul_t);

}


