#include <stdio.h>
#include <math.h>

#include "MB.h"

double pick_dihed(  double atom_i[3],double atom_j[3],double  atom_k[3],double atom_l[3], int flag, double old_value) {
  int alpha;
   
  double vec_ij[3],vec_jk[3],vec_kl[3];
  double out_ij_jk[3],out_jk_kl[3],out_ijkl_jkkl[3];

  double d_ij_jk=0.0,d_jk_kl=0.0,d_ijkl_jkkl=0.0,d_jk=0.0;

  double det=0.0;

  double cs=0.0;
  double theta;
  double pi;

  double out_ijjk_jk[3];
  double in_ijjkjk_jkkl=0.0,in_ijjk_jkkl=0.0;

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

  for(alpha=0;alpha<3;++alpha) {
    cs += out_ij_jk[alpha]*out_jk_kl[alpha];
  }

  if (cs < -1.0 ){
    cs = -1.0;
  }
  else if (cs > 1.0 ) {
    cs = 1.0;
  }

  out_ijkl_jkkl[0] = out_ij_jk[1]*out_jk_kl[2]-out_ij_jk[2]*out_jk_kl[1];
  out_ijkl_jkkl[1] = out_ij_jk[2]*out_jk_kl[0]-out_ij_jk[0]*out_jk_kl[2];
  out_ijkl_jkkl[2] = out_ij_jk[0]*out_jk_kl[1]-out_ij_jk[1]*out_jk_kl[0];

  det = out_ijkl_jkkl[0]*vec_jk[0]+out_ijkl_jkkl[1]*vec_jk[1]+out_ijkl_jkkl[2]*vec_jk[2];
  
  d_ijkl_jkkl += out_ijkl_jkkl[0]*out_ijkl_jkkl[0]+out_ijkl_jkkl[1]*out_ijkl_jkkl[1]+out_ijkl_jkkl[2]*out_ijkl_jkkl[2];
  d_ijkl_jkkl = sqrt(d_ijkl_jkkl);
  d_jk += vec_jk[0]*vec_jk[0]+vec_jk[1]*vec_jk[1]+vec_jk[2]*vec_jk[2];
  d_jk = sqrt(d_jk);

  det = det/(d_ijkl_jkkl*d_jk);
  if (det <0)
    theta = -1.0*acos(cs);
  else
    theta = acos(cs);

  if (flag==0) {
    pi = acos(-1.0);
    if (det<0) {
      theta = 2.0*pi+theta;
    }
  }
  else if (flag==1){
    pi = acos(-1.0);
    if (det<0) {
      theta = 2.0*pi+theta;
    }

    if (theta-old_value>pi) {
      theta-=2.0*pi;
    }
    else if (old_value-theta>pi) {
      theta+=2.0*pi;
    }

  }    
	
  return theta/*rad*/;

}

double pick_angle(double atom_i[3],double atom_j[3],double atom_k[3],int flag,double old_value){
  int i;
  double abv1=0.0,abv2=0.0;
  double v1[3],v2[3],theta=0.0;
  double pi;

  pi=acos(-1.0);
  
  for (i=0;i<3;++i) {
    v1[i]=atom_j[i]-atom_i[i];
    v2[i]=atom_j[i]-atom_k[i];
  }
  
  for (i=0;i<3;++i) {
    abv1+=v1[i]*v1[i];
    abv2+=v2[i]*v2[i];
  }

  for (i=0;i<3;++i) {
    v1[i]=v1[i]/sqrt(abv1);
    v2[i]=v2[i]/sqrt(abv2);
  }

  for (i=0;i<3;++i) 
    theta+=v1[i]*v2[i];
  theta=acos(theta);

  if (flag==1){
    if (theta-old_value>pi) {
      theta-=2.0*pi;
    }
    else if (old_value-theta>pi) {
      theta+=2.0*pi;
    }
  }    

  return theta;
}

double pick_bond_leng(  double atom_i[3],double atom_j[3]) {
  int alpha;
  double len=0.0;

  for (alpha=0;alpha<3;++alpha) {
    len += (atom_j[alpha]-atom_i[alpha])*(atom_j[alpha]-atom_i[alpha]);
  }
  len=sqrt(len);
  return len;
}

