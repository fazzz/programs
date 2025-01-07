
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"

#define ON 1
#define OFF 0

double dih(double atom_i[3],double atom_j[3],double  atom_k[3],double atom_l[3]);

int main(int argc, char *argv[]) {
  int i,j,k;
  int numstep;
  int flag1=OFF,flag2=OFF;

  double pi;
  double p1,h1,p2,h2,p,h;
  double r1[5][3],r2[5][3],crd[5][3];
  char *trjfilename,*inifilename;
  FILE *trjfile,*inifile;

  pi=acos(-1.0);

  if (argc < 8) {
    printf("USAGE: %s numstep tjfilename inifilename p1 h1 p2 h2 \n",argv[0]);
    exit(1);
  }
  numstep = atoi(*++argv);
  trjfilename  = *++argv;
  inifilename  = *++argv;
  p1=atof(*++argv);
  h1=atof(*++argv);
  p2=atof(*++argv);
  h2=atof(*++argv);
  trjfile = efopen(trjfilename,"r");
  inifile = efopen(inifilename,"w");

  for (i=0;i<numstep;++i) {
    for (j=0;j<5;++j)
      for (k=0;k<3;++k)
	fscanf(trjfile,"%lf",&crd[j][k]);

    p=dih(crd[0],crd[1],crd[2],crd[3])*180.0/pi;
    h=dih(crd[1],crd[2],crd[3],crd[4])*180.0/pi;

    if (  p > p1-5.0  && p < p1+5.0 && h > h1-5.0  && h < h1+5.0   ) {
      flag1=ON;
      for (j=0;j<5;++j)
	for (k=0;k<3;++k)
	  r1[j][k]=crd[j][k];

    }
    else if ( p > p2-5.0  && p < p2+5.0 && h > h2-5.0  && h < h2+5.0   ) {
      flag2=ON;
      for (j=0;j<5;++j)
	for (k=0;k<3;++k)
	  r2[j][k]=crd[j][k];
    }
    if ( flag1==ON && flag2==ON )
      break;
  }

  if ( flag1==ON && flag2==ON ) {
    for (j=0;j<5;++j)
      fprintf(inifile,"%12.8lf %12.8lf %12.8lf\n",r1[j][0],r1[j][1],r1[j][2]);
    for (j=0;j<5;++j)
      fprintf(inifile,"%12.8lf %12.8lf %12.8lf\n",r2[j][0],r2[j][1],r2[j][2]);
    fprintf(inifile,"%12.8lf %12.8lf \n",p1,h1);
    fprintf(inifile,"%12.8lf %12.8lf \n",p2,h2);
  }
  else {
    printf("can not detect !");
  }

  fclose(inifile);
  fclose(trjfile);
  
  return 0;
}

double dih(double atom_i[3],double atom_j[3],double  atom_k[3],double atom_l[3]){
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

  pi = acos(-1.0);
  if (det<0) {
    theta = 2.0*pi+theta;
  }
	
  return theta/*rad*/;
}


