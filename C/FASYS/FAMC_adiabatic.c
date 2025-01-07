
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"
#include "RAND.h"
#include "BOXMULL.h"

struct potential {
  double v;
  double w;

  double bond;
  double angl;
  double dihd;

  double eata;
  double guid_potential;
};

void cv(double r1[5][3],double e1[2],double *p1, double *h1, double k1, int sflag, struct potential *energy);
void cv_e(double p1,double h1,double e1[2], double k1,int sflag, double p_t1,double h_t1,double p_t2,double h_t2, struct potential *energy);
double len(double a[3],double b[3]);
double ang(double atom_i[3],double atom_j[3],double atom_k[3]);
double dih(double atom_i[3],double atom_j[3],double  atom_k[3],double atom_l[3]);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int c,c1,c2,sflag='O',gflag='x';
  int numstep,interval;

  double f,beta1,beta2,delta;
  double dx=0.01,dex=0.01,k1=0.0,k2=0.0;
  double r1[5][3],e1[2],v;
  double p1,h1,pi;
  double p_t1,h_t1,p_t2,h_t2;
  double p1_ave,h1_ave;
  double r1_trial[5][3],e1_trial[2],v_trial;
  double p1_trial,h1_trial;
  struct potential ene,ene_trial,ene_eata,ene_eata_trial;

  char *inifilename,*trjfilename,*dihfilename,*eatafilename,*enefilename,*enefile_detailname,*logfilename;
  FILE *inifile,*trjfile,*dihfile,*eatafile,*enefile,*ene_detailfile,*logfile;

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  extern char *optarg;
  extern int optind,opterr,optopt;

  pi=acos(-1.0);

  while((c=getopt(argc,argv,"GOT"))!=-1) {
    switch(c) {
    case 'O':
      sflag='O';
      break;
    case 'T':
      sflag='T';
      break; 
    case 'G':
      gflag='g';
      break; 
   default:
      printf("USAGE: %s numstep interval beta k1 k2 inifilename trjfilename dihfilename eatafilename enefilename enefile_detailname logfilename \n", argv[0]);
      exit(1);
    }
  }

  argc-=optind;
  argv+=optind;

  if (argc < 8) {
    printf("USAGE: %s numstep interval beta1 beta2  k1 inifilename trjfilename dihfilename eatafilename enefilename enefile_detailname logfilename \n", argv[0]);
    exit(1);
  }
  numstep = atoi(*argv);
  interval = atoi(*++argv);
  beta1= atof(*++argv);
  beta2= atof(*++argv);
  k1=atof(*++argv);
  inifilename  = *++argv;
  trjfilename  = *++argv;
  dihfilename  = *++argv;
  eatafilename = *++argv;
  enefilename  = *++argv;
  enefile_detailname = *++argv;
  logfilename  = *++argv;
  inifile=efopen(inifilename,"r");
  for (i=0;i<5;++i)
    for (j=0;j<3;++j)
      fscanf(inifile,"%lf",&r1[i][j]);
  for (i=0;i<2;++i) {
    fscanf(inifile,"%lf",&f);
    e1[i]=f*pi/180.0;
  }
  if (gflag=='g') {
    fscanf(inifile,"%lf",&f);
    p_t1=f/180*pi;
    fscanf(inifile,"%lf",&f);
    h_t1=f/180*pi;
    fscanf(inifile,"%lf",&f);
    p_t2=f/180*pi;
    fscanf(inifile,"%lf",&f);
    h_t2=f/180*pi;
  }
  fclose(inifile);
  trjfile=efopen(trjfilename,"w");
  dihfile=efopen(dihfilename,"w");
  eatafile=efopen(eatafilename,"w");
  enefile=efopen(enefilename,"w");
  ene_detailfile=efopen(enefile_detailname,"w");
  logfile=efopen(logfilename,"w");
  fprintf(logfile,"%s %s %s %s %s\n",inifilename,trjfilename,enefilename,dihfilename,eatafilename);
  fprintf(logfile,"numstep=%d beta1=%4.2e beta2=%4.2e k1=%4.2e k2=%4.2e\n",numstep,beta1,beta2,k1,k2);

  cv(r1,e1,&p1,&h1,k1,sflag,&ene);
  cv_e(p1,h1,e1,k1,sflag,p_t1,h_t1,p_t2,h_t2,&ene_eata);
  for (i=0;i<5;++i)
    for (j=0;j<3;++j)
      r1_trial[i][j]=r1[i][j];
  for (j=0;j<2;++j)
    e1_trial[j]=e1[j];
  p1_ave=0.0;
  h1_ave=0.0;
  for (i=0;i<numstep;++i) {
    for (j=0;j<5;++j) {
      for (k=0;k<3;++k)
	r1_trial[j][k]=r1[j][k]+dx*Box_Muller(i,0.0,1.0);
      cv(r1_trial,e1,&p1_trial,&h1_trial,k1,sflag,&ene_trial);
      delta=(ene_trial.v+ene_trial.eata)-(ene.v+ene.eata);
      if((c1=Metropolis(beta1*delta))==1) {
	for (k=0;k<3;++k)
	  r1[j][k]=r1_trial[j][k];
	p1=p1_trial;
	h1=h1_trial;
	ene.v=ene_trial.v;
	ene.eata=ene_trial.eata;
      }
      p1_ave=((i*5+j)*p1_ave+p1)/(i*5+j+1);
      h1_ave=((i*5+j)*h1_ave+h1)/(i*5+j+1);
    }

    if (i%10==0) {
      for (j=0;j<2;++j) {
    	e1_trial[j]=e1[j]+dex*Box_Muller(i,0.0,1.0);
    	cv(r1_trial,e1,&p1_trial,&h1_trial,k1,sflag,&ene_trial);
	cv_e(p1_ave,h1_ave,e1_trial,k1,sflag,p_t1,h_t1,p_t2,h_t2,&ene_eata_trial);
    	delta=ene_eata_trial.eata-ene_eata.eata;
	if (gflag=='g')
	  delta+=ene_eata_trial.guid_potential-ene_eata.guid_potential;
    	if((c2=Metropolis(beta2*delta))==1) {
    	  e1[j]=e1_trial[j];
    	  ene.eata=ene_trial.eata;
    	  ene_eata.eata=ene_eata_trial.eata;
	  ene_eata.guid_potential=ene_eata_trial.guid_potential;
	  ene.eata=ene_trial.eata;
    	}
      }
      p1_ave=0.0;
      h1_ave=0.0;
    }

    if (i%interval==0) {
      for (j=0;j<5;++j) {
	for (k=0;k<3;++k) {
	  fprintf(trjfile,"%12.8lf ",r1[j][k]);
	if ((j*3+k+1)%10==0)
	  fprintf(trjfile,"\n");
	}
      }
      if (p1>pi)
	p1-=2.0*pi;
      if (h1>pi)
	h1-=2.0*pi;
      if (e1[0]>pi)
	e1[0]-=2.0*pi;
      else if (e1[0]<-1.0*pi)
	e1[0]+=2.0*pi;
      if (e1[1]>pi)
	e1[1]-=2.0*pi;
      else if (e1[1]<-1.0*pi)
	e1[1]+=2.0*pi;
      fprintf(dihfile,"%d %12.8lf %12.8lf \n",i+1,p1*180.0/pi,h1*180.0/pi);
      fprintf(eatafile,"%d %12.8lf %12.8lf \n",i+1,e1[0]*180.0/pi,e1[1]*180.0/pi);
      if (gflag=='g')
	fprintf(enefile,"%d %12.8lf %12.8lf  %12.8lf %d %d\n",i+1,ene.v,ene.eata,ene_eata.guid_potential,c1,c2);
      else
	fprintf(enefile,"%d %12.8lf %12.8lf %d %d\n",i+1,ene.v,ene.eata,c1,c2);
      fprintf(ene_detailfile,"%d %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf \n",i+1,ene.bond,ene.angl,ene.dihd,ene.eata,ene.v);
    }
  }

  fclose(trjfile);
  fclose(dihfile);
  fclose(eatafile);
  fclose(enefile);
  fclose(ene_detailfile);
  fclose(logfile);
  
  return 0;
}

void cv(double r1[5][3],double e1[2],double *p1, double *h1, double k1, int sflag, struct potential *energy) {
  int i;
  int step;
  double pi;
  double kb=100.0,ka=50.0;
  double kd1=0.0,n1=0.0,kd2=0.0,n2=0.0;
  double l_eq=1.53,a_eq;
  double dp,dh;
  double A,B,len2,len6,len12;

  pi=acos(-1.0);
  a_eq=111.0/180*pi;

  (*energy).bond=0.0;
  (*energy).angl=0.0;
    
  for (i=0;i<4;++i)
    (*energy).bond+=0.5*kb*(len(r1[i],r1[i+1])-l_eq)*(len(r1[i],r1[i+1])-l_eq);

  for (i=0;i<3;++i)
    (*energy).angl+=0.5*ka*(ang(r1[i],r1[i+1],r1[i+2])-a_eq)*(ang(r1[i],r1[i+1],r1[i+2])-a_eq);

  if (sflag=='O') {
    kd1=2.5;
    n1=3.0;
    kd2=2.5;
    n2=3.0;
  }
  else {
    kd1=10.0;
    n1=3.0;
    kd2=10.0;
    n2=3.0;
  }

  *p1=dih(r1[0],r1[1],r1[2],r1[3]);
  *h1=dih(r1[1],r1[2],r1[3],r1[4]);

  (*energy).dihd=0.5*kd1*(1.0+cos(n1*(*p1)))+0.5*kd2*(1.0+cos(n2*(*h1)));
  if ((dp=(*p1)-e1[0]) > pi)
    dp-=2.0*pi;
  else if ((dp=(*p1)-e1[0])<-pi)
    dp+=2.0*pi;

  if ((dh=(*h1)-e1[1])>pi)
    dh-=2.0*pi;
  else if ((dh=*(h1)-e1[1])<-pi)
    dh+=2.0*pi;

  (*energy).eata=0.5*k1*dp*dp+0.5*k1*dh*dh;
  (*energy).v=(*energy).bond+(*energy).angl+(*energy).dihd+(*energy).eata;

  (*energy).w=(*energy).eata;

}

void cv_e(double p1,double h1,double e1[2], double k1,int sflag,double p_t1,double h_t1,double p_t2,double h_t2, struct potential *energy) {
  int i;

  double kd1=0.0,n1=0.0,kd2=0.0,n2=0.0;
  double dp,dh;
  double k3=1.0,pi;
  double dp1,dp2,dh1,dh2;

  pi=acos(-1.0);

  (*energy).bond=0.0;
  (*energy).angl=0.0;
  (*energy).dihd=0.0;
  (*energy).v=0.0;
  (*energy).w=0.0;

  if (sflag=='O') {
    kd1=2.5;
    n1=3.0;
    kd2=2.5;
    n2=3.0;
  }
  else {
    kd1=10.0;
    n1=3.0;
    kd2=10.0;
    n2=3.0;
  }

  if ((dp=p1-e1[0]) > pi)
    dp-=2.0*pi;
  else if ((dp=p1-e1[0])<-pi)
    dp+=2.0*pi;

  if ((dh=h1-e1[1])>pi)
    dh-=2.0*pi;
  else if ((dh=h1-e1[1])<-pi)
    dh+=2.0*pi;

  (*energy).eata=0.5*k1*dp*dp+0.5*k1*dh*dh;
  
  if ((dp1=p_t1-e1[0]) > pi)
    dp1-=2.0*pi;
  else if (dp1<-pi)
    dp1+=2.0*pi;

  if ((dh1=h_t1-e1[1])>pi)
    dh1-=2.0*pi;
  else if (dh1<-pi)
    dh1+=2.0*pi;

  if ((dp2=p_t2-e1[0]) > pi)
    dp2-=2.0*pi;
  else if (dp2<-pi)
    dp2+=2.0*pi;

  if ((dh2=h_t2-e1[1])>pi)
    dh2-=2.0*pi;
  else if (dh2<-pi)
    dh2+=2.0*pi;

  (*energy).guid_potential=k3*(dp1*dp1+dh1*dh1)*(dp2*dp2+dh2*dh2);

}

double len(double a[3],double b[3]){
  int i;
  double l=0.0;

  for (i=0;i<3;++i)
    l+=(a[i]-b[i])*(a[i]-b[i]);
  return sqrt(l);
}

double ang(double atom_i[3],double atom_j[3],double atom_k[3]) {
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

  return theta;
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

