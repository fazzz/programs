
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"
#include "RAND.h"
#include "BOXMULL.h"

struct potential {
  double v[2];

  double bond[2];
  double angl[2];
  double dihd[2];
  double eata[3];
};

double cv(double r1[5][3],double r2[5][3],double e1[2], double e2[2], double *p1, double *h1,double *p2, double *h2, double k1, double k2, int sflag, int flag, struct potential *energy);
double len(double a[3],double b[3]);
double ang(double atom_i[3],double atom_j[3],double atom_k[3]);
double dih(double atom_i[3],double atom_j[3],double  atom_k[3],double atom_l[3]);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int c,sflag='O',flag='1';
  int numstep,interval;

  double beta,delta;
  double dx=0.01,dex=0.01,k1=0.0,k2=0.0;
  double r1[5][3],r2[5][3],e1[2],e2[2],v;
  double p1,h1,p2,h2,pi;
  double r1_trial[5][3],r2_trial[5][3],e1_trial[2],e2_trial[2],v_trial;
  double p1_trial,h1_trial,p2_trial,h2_trial;
  struct potential ene;

  char *inifilename,*trjfilename,*dihfilename,*eatafilename,*enefilename,*enefile_detailname,*logfilename;
  FILE *inifile,*trjfile,*dihfile,*eatafile,*enefile,*ene_detailfile,*logfile;

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  extern char *optarg;
  extern int optind,opterr,optopt;

  pi=acos(-1.0);

  while((c=getopt(argc,argv,"OT12"))!=-1) {
    switch(c) {
    case 'O':
      sflag='O';
      break;
    case 'T':
      sflag='T';
      break; 
    case '1':
      flag='1';
      break;
    case '2':
      flag='2';
      break;
   default:
      printf("USAGE: %s numstep interval beta k1 k2 inifilename trjfilename dihfilename eatafilename enefilename enefile_detailname logfilename \n", argv[0]);
      exit(1);
    }
  }

  argc-=optind;
  argv+=optind;

  if (argc < 8) {
    printf("USAGE: %s numstep interval beta k1 k2 inifilename trjfilename dihfilename eatafilename enefilename enefile_detailname logfilename \n", argv[0]);
    exit(1);
  }
  numstep = atoi(*argv);
  interval = atoi(*++argv);
  beta = atof(*++argv);
  k1=atof(*++argv);
  k2=atof(*++argv);
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
  for (i=0;i<5;++i)
    for (j=0;j<3;++j)
      fscanf(inifile,"%lf",&r2[i][j]);
  for (i=0;i<2;++i)
    fscanf(inifile,"%lf",&e1[i]);
  for (i=0;i<2;++i)
    fscanf(inifile,"%lf",&e2[i]);
  fclose(inifile);
  trjfile=efopen(trjfilename,"w");
  dihfile=efopen(dihfilename,"w");
  eatafile=efopen(eatafilename,"w");
  enefile=efopen(enefilename,"w");
  ene_detailfile=efopen(enefile_detailname,"w");
  logfile=efopen(logfilename,"w");
  fprintf(logfile,"%s %s %s\n",inifilename,trjfilename,enefilename);
  fprintf(logfile,"numstep=%d beta=%4.2lf k1=%4.2lf k2=%4.2lf\n",numstep,beta,k1,k2);

  v=cv(r1,r2,e1,e2,&p1,&h1,&p2,&h2,k1,k2,sflag,flag,&ene);
  for (i=0;i<5;++i) {
    for (j=0;j<3;++j) {
      r1_trial[i][j]=r1[i][j];
      r2_trial[i][j]=r2[i][j];
    }
  }
  for (i=0;i<2;++i) {
    e1_trial[i]=e1[i];
    e2_trial[i]=e2[i];
  }

  for (i=0;i<numstep;++i) {
    for (j=0;j<5;++j) {
      for (k=0;k<3;++k)
	r1_trial[j][k]=r1[j][k]+dx*Box_Muller(i,0.0,1.0);
      v_trial=cv(r1_trial,r2_trial,e1_trial,e2_trial,&p1_trial,&h1_trial,&p2_trial,&h2_trial,k1,k2,sflag,flag,&ene);
      delta=v_trial-v;
      if((c=Metropolis(beta*delta))==1) {
	for (k=0;k<3;++k) {
	  r1[j][k]=r1_trial[j][k];
	}
	p1=p1_trial;
	h1=h1_trial;
	v=v_trial;
      }
    }
    for (j=0;j<2;++j) {
      e1_trial[j]=e1[j]+dex*Box_Muller(i,0.0,1.0);
      v_trial=cv(r1_trial,r2_trial,e1_trial,e2_trial,&p1_trial,&h1_trial,&p2_trial,&h2_trial,k1,k2,sflag,flag,&ene);
      delta=v_trial-v;
      if((c=Metropolis(beta*delta))==1) {
	for (k=0;k<2;++k) {
	  e1[k]=e1_trial[k];
	}
	p1=p1_trial;
	h1=h1_trial;
	v=v_trial;
      }
    }
    for (j=0;j<5;++j) {
      for (k=0;k<3;++k)
	r2_trial[j][k]=r2[j][k]+dx*Box_Muller(i,0.0,1.0);
      v_trial=cv(r1_trial,r2_trial,e1_trial,e2_trial,&p1_trial,&h1_trial,&p2_trial,&h2_trial,k1,k2,sflag,flag,&ene);
      delta=v_trial-v;
      if((c=Metropolis(beta*delta))==1) {
	for (k=0;k<3;++k) {
	  r2[j][k]=r2_trial[j][k];
	}
	p1=p1_trial;
	h1=h1_trial;
	v=v_trial;
      }
    }
    for (j=0;j<2;++j) {
      e2_trial[j]=e2[j]+dex*Box_Muller(i,0.0,1.0);
      v_trial=cv(r1_trial,r2_trial,e1_trial,e2_trial,&p1_trial,&h1_trial,&p2_trial,&h2_trial,k1,k2,sflag,flag,&ene);
      delta=v_trial-v;
      if((c=Metropolis(beta*delta))==1) {
	for (k=0;k<2;++k) {
	  e2[k]=e2_trial[k];
	}
	p1=p1_trial;
	h1=h1_trial;
	v=v_trial;
      }
    }

    if (i%interval==0) {
      for (k=0;k<5;++k) {
	for (l=0;l<3;++l) {
	  fprintf(trjfile,"%12.8lf ",r1[k][l]);
	  if ((k*3+l)%10==0)
	    fprintf(trjfile,"\n");
	}
      }
      if (flag=='2') {
	for (k=0;k<5;++k) {
	  for (l=0;l<3;++l) {
	    fprintf(trjfile,"%12.8lf ",r2[k][l]);
	    if ((k*3+l)%10==0)
	      fprintf(trjfile,"\n");
	    }
	}
      }
      if (p1>pi)
	p1-=2.0*pi;
      if (h1>pi)
	h1-=2.0*pi;
      if (p2>pi)
	p2-=2.0*pi;
      if (h2>pi)
	h2-=2.0*pi;
      fprintf(dihfile,"%d %12.8lf %12.8lf\n",i+1,p1*180.0/pi,h1*180.0/pi);
      if (flag=='2') {
	fprintf(dihfile,"%d %12.8lf %12.8lf\n",i+1,p2*180.0/pi,h2*180.0/pi);
      }
      fprintf(eatafile,"%12.8lf %12.8lf\n",e1[0],e1[1]);
      if (flag=='2') {
	fprintf(eatafile,"%12.8lf %12.8lf\n",e2[0],e2[1]);
      }
      fprintf(enefile,"%d %12.8lf %12.8lf %12.8lf %d\n",i+1,v,ene.v[0],ene.v[1],c);
      fprintf(ene_detailfile,"%d %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf \n",i*5+j+1,v,ene.bond[0],ene.bond[1],ene.angl[0],ene.angl[1],ene.dihd[0],ene.angl[1],ene.eata[0],ene.eata[1],ene.eata[2]);
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

double cv(double r1[5][3],double r2[5][3],double e1[2], double e2[2], double *p1, double *h1,double *p2, double *h2, double k1, double k2, int sflag, int flag, struct potential *energy) {
  int i;
  double pi;
  double v=0.0;
  double kb=100.0,ka=50.0;
  double kd1=0.0,n1=0.0,kd2=0.0,n2=0.0;
  double l_eq=1.53,a_eq;
  pi=acos(-1.0);
  a_eq=111.0/180*pi;

  for (i=0;i<2;++i) {
    (*energy).bond[i]=0.0;
    (*energy).angl[i]=0.0;
  }
    
  for (i=0;i<4;++i)
    (*energy).bond[0]+=0.5*kb*(len(r1[i],r1[i+1])-l_eq)*(len(r1[i],r1[i+1])-l_eq);

  for (i=0;i<3;++i)
    (*energy).angl[0]+=0.5*ka*(ang(r1[i],r1[i+1],r1[i+2])-a_eq)*(ang(r1[i],r1[i+1],r1[i+2])-a_eq);

  if (flag=='2') {
    for (i=0;i<4;++i)
      (*energy).bond[1]+=0.5*kb*(len(r2[i],r2[i+1])-l_eq)*(len(r2[i],r2[i+1])-l_eq);

    for (i=0;i<3;++i)
      (*energy).angl[1]+=0.5*ka*(ang(r1[i],r1[i+1],r1[i+2])-a_eq)*(ang(r1[i],r1[i+1],r1[i+2])-a_eq);
  }

  if (sflag='O') {
    kd1=2.5;
    n1=3.0;
    kd2=2.5;
    n2=3.0;
  }
  else {
    kd1=2.0;
    n1=3.0;
    kd2=2.0;
    n2=1.0;
  }

  *p1=dih(r1[0],r1[1],r1[2],r1[3]);
  *h1=dih(r1[1],r1[2],r1[3],r1[4]);
  *p2=dih(r2[0],r2[1],r2[2],r2[3]);
  *h2=dih(r2[1],r2[2],r2[3],r2[4]);

  (*energy).dihd[0]=0.5*kd1*(1.0+cos(n1*(*p1)))+0.5*kd2*(1.0+cos(n2*(*h1)));
  if (flag=='2')
    (*energy).dihd[1]=0.5*kd1*(1.0+cos(n1*(*p2)))+0.5*kd2*(1.0+cos(n2*(*h2)));

  (*energy).eata[0]=0.5*k1*((*p1)-e1[0])*((*p1)-e1[0])+0.5*k1*((*h1)-e1[1])*((*h1)-e1[1]);
  if (flag=='2') {
    (*energy).eata[1]=0.5*k1*((*p2)-e2[0])*((*p2)-e2[0])+0.5*k2*((*h2)-e2[1])*((*h2)-e2[1]);
    (*energy).eata[2]=0.5*k2*(e1[0]-e2[0])*(e1[0]-e2[0])+0.5*k2*(e1[1]-e2[1])*(e1[1]-e2[1]);
  }

  (*energy).v[0]=(*energy).bond[0]+(*energy).angl[0]+(*energy).dihd[0];
  if (flag=='2')
    (*energy).v[1]=(*energy).bond[1]+(*energy).angl[1]+(*energy).dihd[1];

  v=(*energy).bond[0]+(*energy).angl[0]+(*energy).dihd[0]+(*energy).eata[0];
  if (flag=='2')
  v+=(*energy).bond[1]+(*energy).angl[1]+(*energy).dihd[1]+(*energy).eata[1]+(*energy).eata[2];

  return v;
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

