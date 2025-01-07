#include <stdio.h>
#include <math.h>

#include "MB.h"
#include "PT.h"
#include "EF.h"

int countdihed(int flagKOP) {
  int i;
  int numatom,numdihedtype,*numdihed,numdihedtotal;
  int **adpairs;

  numdihedtotal=0;

  if (flagKOP=='P')  numdihedtype=1;
  else if(flagKOP=='O') numdihedtype=2;
  else numdihedtype=5;
  
  numatom=AP.NATOM;

  adpairs=(int **)gcemalloc(sizeof(int *)*5);
  adpairs[0]=(int *)gcemalloc(sizeof(int)*4); // PHI,PSI
  adpairs[1]=(int *)gcemalloc(sizeof(int)*4); // OMEGA
  adpairs[2]=(int *)gcemalloc(sizeof(int)*4); // KI
  adpairs[3]=(int *)gcemalloc(sizeof(int)*8); // ACE
  adpairs[4]=(int *)gcemalloc(sizeof(int)*8); // NME
  numdihed=(int *)gcemalloc(sizeof(int)*5);  
  readdihedpairs(adpairs,numdihed);
  for (i=0;i<4;++i) numdihedtotal+=numdihed[i];

  return numdihedtotal;
}

double *CD(double *crd,int flagKOP) {
  int i,j,k,l;
  int nd;
  int numatom,numdihedtype,*numdihed,numdihedtotal;
  int **adpairs;

  double atom[4][3],dummy;
  double *dihed;
  double pi;
  FILE *logfile;

  numdihedtotal=0;
  pi=acos(-1.0);

  if (flagKOP=='P')  numdihedtype=1;
  else if(flagKOP=='O') numdihedtype=2;
  else numdihedtype=5;
  
  numatom=AP.NATOM;

  adpairs=(int **)gcemalloc(sizeof(int *)*5);
  adpairs[0]=(int *)gcemalloc(sizeof(int)*4); // PHI,PSI
  adpairs[1]=(int *)gcemalloc(sizeof(int)*4); // OMEGA
  adpairs[2]=(int *)gcemalloc(sizeof(int)*4); // KI
  adpairs[3]=(int *)gcemalloc(sizeof(int)*8); // ACE
  adpairs[4]=(int *)gcemalloc(sizeof(int)*8); // NME
  numdihed=(int *)gcemalloc(sizeof(int)*5);  
  readdihedpairs(adpairs,numdihed);
  for (i=0;i<4;++i) numdihedtotal+=numdihed[i];
  //  dihed=(double *)gcerealloc(dihed,sizeof(double)*numdihedtotal);
  dihed=(double *)gcemalloc(sizeof(double)*(numdihedtotal));
 
  nd=0;
  for (i=0;i<numdihedtype;++i) {
    for (j=0;j<numdihed[i];++j) {
      for (k=0;k<4;++k) 
	for (l=0;l<3;++l)
	  atom[k][l]=crd[(adpairs[i][j*4+k])*3+l];
      
      dihed[nd]=pick_dihed(atom[0],atom[1],atom[2],atom[3],0,dummy);
      if (dihed[nd] > pi)  dihed[nd] -= 2.0*pi;
      ++nd;
    }
  }

  logfile=efopen("log_CD.txt","w");
  for (i=0;i<numdihedtype;++i) fprintf(logfile,"%d ",numdihed[i]);
  fprintf(logfile,"\n");
  for (i=0;i<numdihedtype;++i) {
    for (j=0;j<numdihed[i];++j) {
      for (k=0;k<4;++k) 
	fprintf(logfile,"%4d ",adpairs[i][j*4+k]+1);
      for (k=0;k<4;++k) 
	fprintf(logfile,"%s ",AP.IGRAPH[adpairs[i][j*4+k]-1]);
      fprintf(logfile,"\n");
    }
  }
  fclose(logfile);

  return dihed;
}





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

double pick_bond_leng(  double atom_i[3],double atom_j[3]) {
  int alpha;
  double len=0.0;

  for (alpha=0;alpha<3;++alpha) {
    len += (atom_j[alpha]-atom_i[alpha])*(atom_j[alpha]-atom_i[alpha]);
  }
  len=sqrt(len);
  return len;
}

int check_phi_psi_omega_kai1(int atom1,int atom2, int atom3, int atom4 ) {
  int i,j;

  if (strncmp(AP.IGRAPH[atom1],"C\0\0",3) == 0  && strncmp(AP.IGRAPH[atom2],"N\0\0",3) == 0  && strncmp(AP.IGRAPH[atom3],"CA\0",3) == 0 && strncmp(AP.IGRAPH[atom4],"C\0\0",3) ==0 )
    return phi;
  else if (strncmp(AP.IGRAPH[atom1],"C\0\0",3) == 0  && strncmp(AP.IGRAPH[atom2],"N\0\0",3) == 0  && strncmp(AP.IGRAPH[atom3],"CH3",3) == 0 && strncmp(AP.IGRAPH[atom4],"HH31",4) ==0 )
    return phi;
  else if (strncmp(AP.IGRAPH[atom1],"N\0\0",3) == 0 && strncmp(AP.IGRAPH[atom2],"CA\0",3) == 0 && strncmp(AP.IGRAPH[atom3],"C\0\0",3) == 0 && strncmp(AP.IGRAPH[atom4],"N\0\0",3) == 0)
    return psi;
  else if (strncmp(AP.IGRAPH[atom1],"HH31",4) == 0 && strncmp(AP.IGRAPH[atom2],"CH3",3) == 0 && strncmp(AP.IGRAPH[atom3],"C\0\0",3) == 0 && strncmp(AP.IGRAPH[atom4],"N\0\0",3) == 0)
    return psi;
  else if (strncmp(AP.IGRAPH[atom1],"CA\0",3) == 0 && strncmp(AP.IGRAPH[atom2],"C\0\0",3) == 0 && strncmp(AP.IGRAPH[atom3],"N\0\0",3) == 0 && strncmp(AP.IGRAPH[atom4],"CA\0",3) == 0)
    return omega;
  else if (strncmp(AP.IGRAPH[atom1],"CH\0",2) == 0 && strncmp(AP.IGRAPH[atom2],"C\0\0",3) == 0 && strncmp(AP.IGRAPH[atom3],"N\0\0",3) == 0 && strncmp(AP.IGRAPH[atom4],"CA\0",3) == 0)
    return omega;
  else if (strncmp(AP.IGRAPH[atom1],"CA\0",3) == 0 && strncmp(AP.IGRAPH[atom2],"C\0\0",3) == 0 && strncmp(AP.IGRAPH[atom3],"N\0\0",3) == 0 && strncmp(AP.IGRAPH[atom4],"CH\0",2) == 0)
    return omega;
  else if (strncmp(AP.IGRAPH[atom1],"N\0\0",3) == 0 && strncmp(AP.IGRAPH[atom2],"CA\0",3) == 0 && strncmp(AP.IGRAPH[atom3],"CB\0",3) == 0 && strncmp(AP.IGRAPH[atom4],"HB1",3) == 0)
    return kai1;
  else
    return 0;

}

int pick_atom_dihed_pairs(int *atom_dihed_pair,int flag) {
  int i,j;
  int numdihed;
  int numatom,dihedtype;
  FILE *logfile;

  numdihed=0;
  for (i=0;i<AP.MPHIA;++i) {
    if (AP.PA[i][0] >=0 && AP.PA[i][1] >= 0 && AP.PA[i][2] >= 0 && AP.PA[i][3] >=0) {
      dihedtype=check_phi_psi_omega_kai1(abs(AP.PA[i][0])/3,abs(AP.PA[i][1])/3,abs(AP.PA[i][2])/3,abs(AP.PA[i][3])/3);
      if (dihedtype==phi || dihedtype==psi) {
	if (flag == 'p' || flag == 'o' || flag == 'k') {
	  for (j=0;j<4;++j) {
	    atom_dihed_pair[numdihed*4+j]=abs(AP.PA[i][j])/3+1;
	  }
	  ++numdihed;
	}	
      }
      else if (dihedtype==omega) {
	if (flag == 'o' || flag == 'k') {
	  for (j=0;j<4;++j) {
	    atom_dihed_pair[numdihed*4+j]=abs(AP.PA[i][j])/3+1;
	  }
	  ++numdihed;
	}	
      }
    }
  }
  for (i=0;i<AP.NPHIH;++i) {
    if (AP.PH[i][0] >=0 && AP.PH[i][1] >= 0 && AP.PH[i][2] >= 0 && AP.PH[i][3] >=0) {
      dihedtype=check_phi_psi_omega_kai1(abs(AP.PH[i][0])/3,abs(AP.PH[i][1])/3,abs(AP.PH[i][2])/3,abs(AP.PH[i][3])/3);
      if (dihedtype==phi || dihedtype==psi) {
	if (flag == 'p' || flag == 'o' || flag == 'k') {
	  for (j=0;j<4;++j) {
	    atom_dihed_pair[numdihed*4+j]=abs(AP.PH[i][j])/3+1;
	  }
	  ++numdihed;
	}	
      }
      else if (dihedtype==omega) {
	if (flag == 'o' || flag == 'k') {
	  for (j=0;j<4;++j) {
	    atom_dihed_pair[numdihed*4+j]=abs(AP.PH[i][j])/3+1;
	  }
	  ++numdihed;
	}	
      }
    }
  }
  if (flag == 'k') {
    for (i=0;i<AP.NPHIH;++i) {
      if (AP.PH[i][0] >=0 && AP.PH[i][1] >= 0 && AP.PH[i][2] >= 0 && AP.PH[i][3] >=0) {
	dihedtype=check_phi_psi_omega_kai1(abs(AP.PH[i][0])/3,abs(AP.PH[i][1])/3,abs(AP.PH[i][2])/3,abs(AP.PH[i][3])/3);
	if (dihedtype==kai1 ) {
	  for (j=0;j<4;++j) {
	    atom_dihed_pair[numdihed*4+j]=abs(AP.PH[i][j])/3+1;
	  }
	  ++numdihed;
	}	
      }
    }
  }

  logfile=efopen("log.txt","w");
  fprintf(logfile,"%d\n",numdihed);
  for (i=0;i<numdihed;++i) {
    for (j=0;j<4;++j) {
      fprintf(logfile,"%d ",atom_dihed_pair[i*4+j]);
    }
    for (j=0;j<4;++j) {
      fprintf(logfile,"%s ",AP.IGRAPH[atom_dihed_pair[i*4+j]-1]);
    }
    fprintf(logfile,"\n");
  }
  fclose(logfile);

  return numdihed;
}

void pick_dihed_all(double *coord,double *dihed,int numdihed,int *atom_dihed_pair) {
  int i,j;
  double atom_i[3],atom_j[3],atom_k[3],atom_l[3];
  double theta;
  double pi;

  for (i=0;i<numdihed;++i) {
    for (j=0;j<3;++j) {
      atom_i[j]=coord[(atom_dihed_pair[i*4]-1)*3+j];
      atom_j[j]=coord[(atom_dihed_pair[i*4+1]-1)*3+j];
      atom_k[j]=coord[(atom_dihed_pair[i*4+2]-1)*3+j];
      atom_l[j]=coord[(atom_dihed_pair[i*4+3]-1)*3+j];
    }
    
    pi=acos(-1.0);
    theta=pick_dihed(atom_i,atom_j,atom_k,atom_l,0,0.0);
    if (theta > pi)
      dihed[i] = theta -2.0*pi;
    else
      dihed[i] = theta;
  }
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

double wraped_angle(double angle_rad, double pi) {
  double angle_rad_wrapped;

  if (angle_rad>pi)
    angle_rad_wrapped-=2.0*pi;
  else if (angle_rad<-1.0*pi)
    angle_rad_wrapped+=2.0*pi;
  else
    angle_rad_wrapped=angle_rad;

  return angle_rad_wrapped;
}
