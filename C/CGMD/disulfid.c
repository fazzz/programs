
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "disulfid.h"
//#include "ParmTop.h"
#include "/home/yamamori/work/programs/readParmtop/PTL.h"
#include "EF.h"

#define UNIT 4.184070*100.0

void disulfid_csdih(double atom_i[3],double atom_j[3],double  atom_k[3],double atom_l[3], double *cs, double *sn);
int disulfid_outprod(double v1[3],double v2[3],double *v1x2);
double disulfid_inprod(double *v1, double *v2, int n);
double disulfid_len(double a[3],double b[3]);

#define ON 0
#define OFF 1

int disulfid_count_numcys(void) {
  int i,j,k,l;
  int nres=0;
  int ns=0,nc=0;
  int *index_S;
  int *index_CT;
  int numcys;
  int findS[2];
  int nbond,nangl,ndihed;

  index_S=(int *)gcemalloc(sizeof(int)*1);
  index_CT=(int *)gcemalloc(sizeof(int)*1);
  for (i=0;i<AP.NATOM;++i) {
    if (AP.IPRES[nres+1]-1==i && i>AP.IPRES[0]-1 ) {
      ++nres;
    }
    if (strncmp(AP.LABERES[nres],"CYX",3)==0 /*&& strncmp(AP.IGRAPH[i],"S",1)==0*/ ) {
      if (strncmp(AP.IGRAPH[i],"S",1)==0 ) {
	index_S[ns]=i;
	++ns;
	index_S=(int *)gcerealloc(index_S,sizeof(int)*(ns+1));
      }
    }
    if (strncmp(AP.LABERES[nres],"CYX",3)==0 /*&& strncmp(AP.IGRAPH[i],"CB",2)==0*/ ) {
      if (strncmp(AP.IGRAPH[i],"CB",2)==0 ) {
	index_CT[nc]=i;
	++nc;
	index_CT=(int *)gcerealloc(index_CT,sizeof(int)*(nc+1));
      }
    }
  }
  if (ns==nc)
    numcys=ns;
  else {
    printf("error\n");
    exit(1);
  }

  if (numcys%2!=0) {
    printf("error\n");
    exit(1);
  }
  else 
    numcys=(int)numcys/2;

  return numcys;
}

int disulfid_read_parm(int **atoms_b,double *K_bond,double *eq_bond,int *nbond,
		       int **atoms_a,double *K_angl,double *eq_angl,int *nangl,
		       int **atoms_d,double *V_dihe,double *n_dihe,double *the_dihe,int *ndihed
		       ) {
  int i,j,k,l;
  int nres=0;
  int ns=0,nc=0;
  int numcys;
  int *index_S;
  int *index_CT;
  int findS[2];
  //  int nbond,nangl,ndihed;

  index_S=(int *)gcemalloc(sizeof(int)*1);
  index_CT=(int *)gcemalloc(sizeof(int)*1);
  for (i=0;i<AP.NATOM;++i) {
    if (AP.IPRES[nres+1]-1==i && i>AP.IPRES[0]-1 ) {
      ++nres;
    }
    if (strncmp(AP.LABERES[nres],"CYX",3)==0 /*&& strncmp(AP.IGRAPH[i],"S",1)==0*/ ) {
      if (strncmp(AP.IGRAPH[i],"S",1)==0 ) {
	index_S[ns]=i;
	++ns;
	index_S=(int *)gcerealloc(index_S,sizeof(int)*(ns+1));
      }
    }
    if (strncmp(AP.LABERES[nres],"CYX",3)==0 /*&& strncmp(AP.IGRAPH[i],"CB",2)==0*/ ) {
      if (strncmp(AP.IGRAPH[i],"CB",2)==0 ) {
	index_CT[nc]=i;
	++nc;
	index_CT=(int *)gcerealloc(index_CT,sizeof(int)*(nc+1));
      }
    }
  }
  if (ns==nc)
    numcys=ns;
  else {
    printf("error\n");
    exit(1);
  }

  if (numcys%2!=0) {
    printf("error\n");
    exit(1);
  }
  else 
    numcys=(int)numcys/2;
  
  *nbond=0;
  for (i=0;i<AP.MBONA;++i) {
    findS[0]=OFF;
    findS[1]=OFF;
    for (j=0;j<numcys*2;++j) {
      if (AP.BA[i][0]==index_S[j]*3) {
	findS[0]=ON;
	break;
      }
    }
    if (findS[0]==ON) {
      findS[1]=OFF;
      for (k=0;k<numcys*2;++k) {
	if (AP.BA[i][1]==index_S[k]*3) {
	  findS[1]=ON;
	  break;
	}
      }
      if (findS[1]==ON) {
	for (l=0;l<3;++l) {
	  atoms_b[*nbond][l]=AP.BA[i][l];
	}
	++(*nbond);
      }
    }
  }

  *nangl=0;
  for (i=0;i<AP.MTHETA;++i) {
    for (j=0;j<3;++j) findS[j]=OFF;
    // CT-CT-S
    for (j=0;j<numcys*2;++j) {
      if (AP.TA[i][0]==index_CT[j]*3) {
	findS[0]=ON;
	break;
      }
    }
    if (findS[0]==ON) {
      for (j=0;j<numcys*2;++j) {
	if (AP.TA[i][1]==index_CT[j]*3) {
	  findS[1]=ON;
	  break;
	}
      }
    }
    if (findS[1]==ON) {
      for (j=0;j<numcys*2;++j) {
	if (AP.TA[i][2]==index_S[j]*3) {
	  findS[2]=ON;
	  break;
	}
      }
    }
    if (findS[2]==OFF) {
      for (j=0;j<3;++j) findS[j]=OFF;
      // CT-S-CT
      for (j=0;j<numcys*2;++j) {
	if (AP.TA[i][0]==index_CT[j]*3) {
	  findS[0]=ON;
	  break;
	}
      }
      if (findS[0]==ON) {
	for (j=0;j<numcys*2;++j) {
	  if (AP.TA[i][1]==index_S[j]*3) {
	    findS[1]=ON;
	    break;
	  }
	}
      }
      if (findS[1]==ON) {
	for (j=0;j<numcys*2;++j) {
	  if (AP.TA[i][2]==index_CT[j]*3) {
	    findS[2]=ON;
	    break;
	  }
	}
      }
    }
    if (findS[2]==OFF) {
      for (j=0;j<3;++j) findS[j]=OFF;
      // CT-S-S
      for (j=0;j<numcys*2;++j) {
	if (AP.TA[i][0]==index_CT[j]*3) {
	  findS[0]=ON;
	  break;
	}
      }
      if (findS[0]==ON) {
	for (j=0;j<numcys*2;++j) {
	  if (AP.TA[i][1]==index_S[j]*3) {
	    findS[1]=ON;
	    break;
	  }
	}
      }
      if (findS[1]==ON) {
	for (j=0;j<numcys*2;++j) {
	  if (AP.TA[i][2]==index_S[j]*3) {
	    findS[2]=ON;
	    break;
	  }
	}
      }
    }
    if (findS[2]==OFF) {
      for (j=0;j<3;++j) findS[j]=OFF;
      // S-S-CT
      for (j=0;j<numcys*2;++j) {
	if (AP.TA[i][0]==index_S[j]*3) {
	  findS[0]=ON;
	  break;
	}
      }
      if (findS[0]==ON) {
	for (j=0;j<numcys*2;++j) {
	  if (AP.TA[i][1]==index_S[j]*3) {
	    findS[1]=ON;
	    break;
	  }
	}
      }
      if (findS[1]==ON) {
	for (j=0;j<numcys*2;++j) {
	  if (AP.TA[i][2]==index_CT[j]*3) {
	    findS[2]=ON;
	    break;
	  }
	}
      }
    }
    if (findS[2]==ON) {
      for (l=0;l<4;++l) {
	atoms_a[*nangl][l]=AP.TA[i][l];
      }
      ++(*nangl);
    }
  }

  *ndihed=0;
  for (i=0;i<AP.MTHETA;++i) {
    for (j=0;j<4;++j) findS[j]=OFF;
    // CT-S-S-CT
    for (j=0;j<numcys*2;++j) {
      if (abs(AP.PA[i][0])==index_CT[j]*3) {
	findS[0]=ON;
	break;
      }
    }
    if (findS[0]==ON) {
      for (j=0;j<numcys*2;++j) {
	if (abs(AP.PA[i][1])==index_S[j]*3) {
	  findS[1]=ON;
	  break;
	}
      }
    }
    if (findS[1]==ON) {
      for (j=0;j<numcys*2;++j) {
	if (abs(AP.PA[i][2])==index_S[j]*3) {
	  findS[2]=ON;
	  break;
	}
      }
    }
    if (findS[2]==ON) {
      for (j=0;j<numcys*2;++j) {
	if (abs(AP.PA[i][3])==index_CT[j]*3) {
	  findS[3]=ON;
	  break;
	}
      }
    }
    if (findS[3]==ON) {
      for (j=0;j<5;++j) {
	atoms_d[*ndihed][j]=abs(AP.PA[i][j]);
      }
      ++(*ndihed);
    }
  }

  for (i=0;i<*nbond;++i) {
    K_bond[i]=AP.RK[atoms_b[i][2]-1];
    eq_bond[i]=AP.REQ[atoms_b[i][2]-1];
  }
  for (i=0;i<*nangl;++i) {
    K_angl[i]=AP.TK[atoms_a[i][3]-1];
    eq_angl[i]=AP.TEQ[atoms_a[i][3]-1];
  }
  for (i=0;i<*ndihed;++i) {
    V_dihe[i]=AP.PK[atoms_d[i][4]-1];
    n_dihe[i]=(double)(AP.PN[atoms_d[i][4]-1]);
    the_dihe[i]=AP.PHASE[atoms_d[i][4]-1];
  }

  return numcys;
}

void disulfid_calc_pf(double *crd,double *p,double *f,
		      int **atoms_b,double *K_bond,double *eq_bond,int nbond,
		      int flagb,
		      int **atoms_a,double *K_angl,double *eq_angl,int nangl,
		      int flaga,
		      int **atoms_d,double *V_dihe,
		      double *n_dihe,double *the_dihe,int ndihed,
		      int flagd
		      ) {
  int i,j,k;
  double p_b,*f_b,p_a,*f_a,p_d,*f_d;
  
  f_b=(double *)gcemalloc(sizeof(double)*AP.NATOM*3);
  f_a=(double *)gcemalloc(sizeof(double)*AP.NATOM*3);
  f_d=(double *)gcemalloc(sizeof(double)*AP.NATOM*3);

  if (flagb==1)
    disulfid_calcBOND(crd,&p_b,f_b,atoms_b,K_bond,eq_bond,nbond);
  if (flaga==1)
    disulfid_calcANGL(crd,&p_a,f_a,atoms_a,K_angl,eq_angl,nangl);
  if (flagd==1)
    disulfid_calcDIHE(crd,&p_d,f_d,atoms_d,V_dihe,n_dihe,the_dihe,ndihed);

  *p=0.0;
  *p+=p_b+p_a+p_d;
  for (i=0;i<AP.NATOM;++i) {
    for (j=0;j<3;++j)
      f[i*3+j]=f_b[i*3+j]+f_a[i*3+j]+f_d[i*3+j];
  }
}

void disulfid_calcBOND(double *crd,double *p_b,double *f_b,
		       int **atoms_b,double *K_bond,double *eq_bond,
		       int numbond) {
  int i,j,k;
  double f;
  double lenij;
  double atom[2][3];
  
  *p_b = 0.0;
  for (i=0;i<AP.NATOM*3;++i) f_b[i] = 0.0;
  for (i=0;i<numbond;++i) {
    for (j=0;j<2;++j)
      for (k=0;k<3;++k)
	atom[j][k]=crd[atoms_b[i][j]+k];
    
    lenij = disulfid_len(atom[0],atom[1]);
    *p_b+=K_bond[i]*(lenij-eq_bond[i])*(lenij-eq_bond[i]);
    for (j=0;j<3;++j) {
      f = -2.0*K_bond[i]*(lenij-eq_bond[i])*(atom[1][j]-atom[0][j])/lenij*UNIT;
      f_b[atoms_b[i][0]+j] += f;
      f_b[atoms_b[i][1]+j] += -f;
    }
  }
}

void disulfid_calcANGL(double *crd,double *p_a,double *f_a,
		       int **atoms_a,double *K_angl,double *eq_angl,
		       int numangl){
  int i,j,k,l;
  
  double atom[3][3];

  double lenij,lenkj;
  double vij[3],vkj[3];
  double cosijk,angijk;
  double f1,f2;

  *p_a = 0.0;
  for (i=0;i<AP.NATOM;++i) f_a[i] = 0.0;  
  for (i=0;i<numangl;++i) {
    for (j=0;j<3;++j)
      for (k=0;k<3;++k)
	atom[j][k]=crd[atoms_a[i][j]+k];


    lenij = disulfid_len(atom[0],atom[1]);
    lenkj = disulfid_len(atom[2],atom[1]);
    for (j=0;j<3;++j) {
      vij[j]=atom[1][j]-atom[0][j];
      vkj[j]=atom[1][j]-atom[2][j];
    }
    cosijk=disulfid_inprod(vij,vkj,3);
    cosijk=cosijk/lenij/lenkj;
    angijk = acos(cosijk);

    *p_a += K_angl[i]*(angijk-eq_angl[i])*(angijk-eq_angl[i]);

    for (j=0;j<3;++j) {
      f1 = -2.0*K_angl[i]*(angijk-eq_angl[i])/(lenij*sin(angijk))*(vkj[j]/lenkj-cosijk*vij[j]/lenij)*UNIT;
      f2 = -2.0*K_angl[i]*(angijk-eq_angl[i])/(lenkj*sin(angijk))*(vij[j]/lenij-cosijk*vkj[j]/lenkj)*UNIT;

      f_a[atoms_a[i][0]+j] += f1;
      f_a[atoms_a[i][1]+j] += -f1-f2;
      f_a[atoms_a[i][1]+j] += f2;

    }
  }
}

void disulfid_calcDIHE(double *crd,double *p_d,double *f_d,
		       int **atoms_d,double *V_dihe,double *n_dihe,double *the_dihe,
		       int numdihed) {
  int i,j,k,l;

  double fa,fb[3],fc[3];
  double *n1,*n2,ln1,ln2;
  double vij[3],vkj[3],vki[3],vjl[3],vji[3],vik[3],vkl[3];
  double op1[3],op2[3],op3[3],op4[3],op5[3],op6[3];
  
  double atom[4][3];
  double cosdih,sindih;
  double dihedang;

  n1=(double *)gcemalloc(sizeof(double)*3);
  n2=(double *)gcemalloc(sizeof(double)*3);

  *p_d=0.0;
  for (i=0;i<AP.NATOM;++i) f_d[i] = 0.0;
  for (i=0;i<numdihed;++i) {
    for (j=0;j<4;++j)
      for (k=0;k<3;++k)
	atom[j][k]=crd[atoms_d[i][j]+k];
  
    dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
    *p_d += V_dihe[i]*(1.0+cos(n_dihe[i]*dihedang-the_dihe[i]));

    for (j=0;j<3;++j) {
      vij[j] = atom[1][j]-atom[0][j];
      vkj[j] = atom[1][j]-atom[2][j];
      vki[j] = atom[0][j]-atom[2][j];
      vjl[j] = atom[3][j]-atom[1][j];
      vji[j] = atom[0][j]-atom[1][j];
      vik[j] = atom[2][j]-atom[0][j];
      vkl[j] = atom[3][j]-atom[2][j];
    }

    disulfid_outprod(vij,vkj,n1);
    disulfid_outprod(vkj,vkl,n2);
    ln1=sqrt(disulfid_inprod(n1,n1,3));
    ln2=sqrt(disulfid_inprod(n2,n2,3));
  
    disulfid_csdih(atom[0],atom[1],atom[2],atom[3],&cosdih,&sindih);

    if (n_dihe[i]==1) fa=-V_dihe[i]*n_dihe[i]*cos(the_dihe[i]);
    else if (n_dihe[i]==2) fa=-V_dihe[i]*n_dihe[i]*2.0*cosdih*cos(the_dihe[i]);
    else if (n_dihe[i]==3) fa=-V_dihe[i]*n_dihe[i]*(-4.0*sindih*sindih+3.0)*cos(the_dihe[i]);
    else if (n_dihe[i]==4) fa=-V_dihe[i]*n_dihe[i]*4.0*(cosdih*(2.0*cosdih*cosdih-1.0))*cos(the_dihe[i]);
    else {
      printf("error:periodicity must be 1~4\n");
      exit(1);
    }

    for (j=0;j<3;++j) fb[j]=(n2[j]/ln2-cosdih*n1[j]/ln1)/ln1;
    for (j=0;j<3;++j) fc[j]=(n1[j]/ln1-cosdih*n2[j]/ln2)/ln2;

    disulfid_outprod(fb,vkj,op1);
    disulfid_outprod(fc,vki,op2);
    disulfid_outprod(fb,vik,op3);
    disulfid_outprod(fb,vij,op4);
    disulfid_outprod(fc,vjl,op5);
    disulfid_outprod(fc,vkj,op6);

    for (j=0;j<3;++j) {
      f_d[atoms_d[i][0]+j] += fa*op1[j]*UNIT;
      f_d[atoms_d[i][1]+j] += fa*(-op2[j]+op3[j])*UNIT;
      f_d[atoms_d[i][2]+j] += fa*(-op4[j]+op5[j])*UNIT;
      f_d[atoms_d[i][3]+j] += fa*op6[j]*UNIT;
    }
  }
}

void disulfid_csdih(double atom_i[3],double atom_j[3],double  atom_k[3],double atom_l[3], double *cs, double *sn){
  int alpha;
  double vec_ij[3],vec_jk[3],vec_kl[3];
  double out_ij_jk[3],out_jk_kl[3],out_ijkl_jkkl[3];

  double d_ij_jk=0.0,d_jk_kl=0.0,d_ijkl_jkkl=0.0,d_jk=0.0;

  double det=0.0;

  double pi;

  double out_ijjk_jk[3];
  double in_ijjkjk_jkkl=0.0,in_ijjk_jkkl=0.0;

  *cs=0.0;  

  for (alpha=0;alpha<3;++alpha) {
    vec_ij[alpha] = atom_j[alpha]-atom_i[alpha];
    vec_jk[alpha] = atom_k[alpha]-atom_j[alpha];
    vec_kl[alpha] = atom_l[alpha]-atom_k[alpha];
  }

  out_ij_jk[0]=vec_ij[1]*vec_jk[2]-vec_ij[2]*vec_jk[1];
  out_ij_jk[1]=vec_ij[2]*vec_jk[0]-vec_ij[0]*vec_jk[2];
  out_ij_jk[2]=vec_ij[0]*vec_jk[1]-vec_ij[1]*vec_jk[0];
			
  out_jk_kl[0]=vec_jk[1]*vec_kl[2]-vec_jk[2]*vec_kl[1];
  out_jk_kl[1]=vec_jk[2]*vec_kl[0]-vec_jk[0]*vec_kl[2];
  out_jk_kl[2]=vec_jk[0]*vec_kl[1]-vec_jk[1]*vec_kl[0];

  d_ij_jk += out_ij_jk[0]*out_ij_jk[0]+out_ij_jk[1]*out_ij_jk[1]+out_ij_jk[2]*out_ij_jk[2];
  d_jk_kl += out_jk_kl[0]*out_jk_kl[0]+out_jk_kl[1]*out_jk_kl[1]+out_jk_kl[2]*out_jk_kl[2];

  d_ij_jk = sqrt(d_ij_jk);
  d_jk_kl = sqrt(d_jk_kl);

  for (alpha=0;alpha<3;++alpha) {
    out_ij_jk[alpha] = out_ij_jk[alpha]/d_ij_jk;
    out_jk_kl[alpha] = out_jk_kl[alpha]/d_jk_kl;
  }

  for(alpha=0;alpha<3;++alpha) *cs += out_ij_jk[alpha]*out_jk_kl[alpha];

  if (*cs < -1.0 ) *cs = -1.0;
  else if (*cs > 1.0 ) *cs = 1.0;

  out_ijkl_jkkl[0] = out_ij_jk[1]*out_jk_kl[2]-out_ij_jk[2]*out_jk_kl[1];
  out_ijkl_jkkl[1] = out_ij_jk[2]*out_jk_kl[0]-out_ij_jk[0]*out_jk_kl[2];
  out_ijkl_jkkl[2] = out_ij_jk[0]*out_jk_kl[1]-out_ij_jk[1]*out_jk_kl[0];

  det = out_ijkl_jkkl[0]*vec_jk[0]+out_ijkl_jkkl[1]*vec_jk[1]+out_ijkl_jkkl[2]*vec_jk[2];
  
  d_ijkl_jkkl += out_ijkl_jkkl[0]*out_ijkl_jkkl[0]+out_ijkl_jkkl[1]*out_ijkl_jkkl[1]+out_ijkl_jkkl[2]*out_ijkl_jkkl[2];
  d_ijkl_jkkl = sqrt(d_ijkl_jkkl);
  d_jk += vec_jk[0]*vec_jk[0]+vec_jk[1]*vec_jk[1]+vec_jk[2]*vec_jk[2];
  d_jk = sqrt(d_jk);

  det = det/(d_ijkl_jkkl*d_jk);
  if (det <0) {
    *sn = 1.0 - (*cs)*(*cs);
    *sn = -1.0*sqrt(*sn);
  }
  else {
    *sn = 1.0 - (*cs)*(*cs);
    *sn = sqrt(*sn);
  }
}

int disulfid_outprod(double v1[3],double v2[3],double *v1x2) {

  v1x2[0]=v1[1]*v2[2]-v1[2]*v2[1];
  v1x2[1]=v1[2]*v2[0]-v1[0]*v2[2];
  v1x2[2]=v1[0]*v2[1]-v1[1]*v2[0];

}

double disulfid_inprod(double *v1, double *v2, int n) {
  int i;
  double in=0.0;

  for (i=0;i<n;++i)
    in += v1[i]*v2[i];

  return in;
}

double disulfid_len(double a[3],double b[3]){
  int i;
  double l=0.0;

  for (i=0;i<3;++i)
    l+=(a[i]-b[i])*(a[i]-b[i]);
  return sqrt(l);
}

double disulfid_check_force_calc(double *crd,double *f,double dx,
				 int **atoms_b,double *K_bond,double *eq_bond,int nbond,
				 int flagb,
				 int **atoms_a,double *K_angl,double *eq_angl,int nangl,
				 int flaga,
				 int **atoms_d,double *V_dihe,
				 double *n_dihe,double *the_dihe,int ndihed,
				 int flagd) {
  int i,j,k;
  double *f_b,*f_a,*f_d;

  f_b=(double *)gcemalloc(sizeof(double)*AP.NATOM*3);
  f_a=(double *)gcemalloc(sizeof(double)*AP.NATOM*3);
  f_d=(double *)gcemalloc(sizeof(double)*AP.NATOM*3);

  for (i=0;i<AP.NATOM*3;++i) {
    f_b[i]=0.0;
    f_a[i]=0.0;
    f_d[i]=0.0;
  }

  if (flagb==1) {
    for (i=0;i<nbond;++i) {
      for (j=0;j<2;++j) {
	disulfid_calcBOND_check(crd,(int)(atoms_b[i][j]/3),dx,
				f_b,atoms_b,K_bond,eq_bond,nbond);
      }
    }
  }
  if (flaga==1) {
    for (i=0;i<nangl;++i) {
      for (j=0;j<3;++j) {
	disulfid_calcANGL_check(crd,(int)(atoms_a[i][j]/3),dx,
				f_a,atoms_a,K_angl,eq_angl,nangl);
      }
    }
  }
  if (flagd==1) {
    for (i=0;i<ndihed;++i) {
      for (j=0;j<4;++j) {
	disulfid_calcDIHE_check(crd,(int)(atoms_d[i][j]/3),dx,
				f_d,atoms_d,V_dihe,n_dihe,the_dihe,ndihed);
      }
    }
  }

  for (i=0;i<AP.NATOM;++i) {
    for (j=0;j<3;++j)
      f[i*3+j]=f_b[i*3+j]+f_a[i*3+j]+f_d[i*3+j];
  }
}

void disulfid_calcBOND_check(double *crd,int numatom,double dx,
			     double *f_b,int **atoms_b,
			     double *K_bond,double *eq_bond,
			     int numbond) {
  int i,j,k,a;
  double *crdd;
  double p,pd;
  double f[3];
  double lenij,lenijd;
  double atom[2][3],atomd[2][3];

  crdd=(double *)gcemalloc(sizeof(double)*AP.NATOM*3);

  for (a=0;a<3;++a) {
    for (i=0;i<AP.NATOM*3;++i) crdd[i]=crd[i];
    crdd[numatom*3+a]+=dx;
  
    p = 0.0;
    pd = 0.0;
    for (i=0;i<numbond;++i) {
      for (j=0;j<2;++j) {
	for (k=0;k<3;++k) {
	  atom[j][k]=crd[atoms_b[i][j]+k];
	  atomd[j][k]=crdd[atoms_b[i][j]+k];
	}
      }
      
      lenij = disulfid_len(atom[0],atom[1]);
      lenijd = disulfid_len(atomd[0],atomd[1]);
      p+=K_bond[i]*(lenij-eq_bond[i])*(lenij-eq_bond[i]);
      pd+=K_bond[i]*(lenijd-eq_bond[i])*(lenijd-eq_bond[i]);
    }
    f[a] = (pd-p)/dx*UNIT;
    f_b[numatom*3+a] += f[a];
  }
}

void disulfid_calcANGL_check(double *crd,int numatom,double dx,
			     double *f_a,
			     int **atoms_a,double *K_angl,double *eq_angl,
			     int numangl){
  int i,j,k,l,a;
  double *crdd;
  double p,pd;
  double f[3];  
  double atom[3][3],atomd[3][3];
  double lenij,lenkj,lenijd,lenkjd;
  double angijk,angijkd;
  double vij[3],vkj[3],vijd[3],vkjd[3];
  double cosijk,cosijkd;

  crdd=(double *)gcemalloc(sizeof(double)*AP.NATOM*3);

  for (a=0;a<3;++a) {
    for (i=0;i<AP.NATOM*3;++i) crdd[i]=crd[i];
    crdd[numatom*3+a]+=dx;

    p = 0.0;
    pd= 0.0;
    for (i=0;i<numangl;++i) {
      for (j=0;j<3;++j) {
	for (k=0;k<3;++k) {
	  atom[j][k]=crd[atoms_a[i][j]+k];
	  atomd[j][k]=crdd[atoms_a[i][j]+k];
	}
      }

      lenij = disulfid_len(atom[0],atom[1]);
      lenkj = disulfid_len(atom[2],atom[1]);
      for (j=0;j<3;++j) {
	vij[j]=atom[1][j]-atom[0][j];
	vkj[j]=atom[1][j]-atom[2][j];
      }
      cosijk=disulfid_inprod(vij,vkj,3);
      cosijk=cosijk/lenij/lenkj;
      angijk = acos(cosijk);

      lenijd = disulfid_len(atomd[0],atomd[1]);
      lenkjd = disulfid_len(atomd[2],atomd[1]);
      for (j=0;j<3;++j) {
	vijd[j]=atomd[1][j]-atomd[0][j];
	vkjd[j]=atomd[1][j]-atomd[2][j];
      }
      cosijkd=disulfid_inprod(vijd,vkjd,3);
      cosijkd=cosijkd/lenijd/lenkjd;
      angijkd = acos(cosijkd);

      p += K_angl[i]*(angijk-eq_angl[i])*(angijk-eq_angl[i]);
      pd += K_angl[i]*(angijkd-eq_angl[i])*(angijkd-eq_angl[i]);

    }
    f[a] = (pd-p)/dx*UNIT;
    f_a[numatom*3+a] += f[a];
  }
}

void disulfid_calcDIHE_check(double *crd,int numatom,double dx,
			     double *f_d,
			     int **atoms_d,double *V_dihe,double *n_dihe,double *the_dihe,
			     int numdihed) {
  int i,j,k,l,a;
  double *crdd;
  double p,pd;
  double f[3];
  double dihedang,dihedangd;  
  double atom[4][3],atomd[4][3];

  crdd=(double *)gcemalloc(sizeof(double)*AP.NATOM*3);

  for (a=0;a<3;++a) {
    for (i=0;i<AP.NATOM*3;++i) crdd[i]=crd[i];
    crdd[numatom*3+a]+=dx;

    p=0.0;
    pd=0.0;
    for (i=0;i<AP.NATOM;++i) f_d[i] = 0.0;
    for (i=0;i<numdihed;++i) {
      for (j=0;j<4;++j) {
	for (k=0;k<3;++k) {
	  atom[j][k]=crd[atoms_d[i][j]+k];
	  atomd[j][k]=crdd[atoms_d[i][j]+k];
	}
      }

      dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
      dihedangd = pick_dihed(atomd[0],atomd[1],atomd[2],atomd[3],0,0.0);
      p += V_dihe[i]*(1.0+cos(n_dihe[i]*dihedang-the_dihe[i]));
      pd += V_dihe[i]*(1.0+cos(n_dihe[i]*dihedangd-the_dihe[i]));
    }
    f[a] = (pd-p)/dx*UNIT;
    f_d[numatom*3+a] += f[a];
  }
}
