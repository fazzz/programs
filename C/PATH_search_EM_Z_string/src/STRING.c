
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "STRING.h"
#include "CSI.h"
#include "TOPO.h"
#include "PT.h"
#include "MB.h"
#include "LA.h"
#include "FF.h"
#include "EF.h"

double z_string_sample(double *p_stat/*,double *p_evoluted,double *f*/,int numpoint,double dt) {
  int i;

  double *fe;
  double *p_evoluted;
  
  //  fe=(double *)gcemalloc(sizeof(double)*numpoint*2);
  p_evoluted=(double *)gcemalloc(sizeof(double)*numpoint*2);


  // fe=(double *)emalloc(sizeof(double)*numpoint*2);
  //  p_evoluted=(double *)emalloc(sizeof(double)*numpoint*2);
  
  for (i=0;i<numpoint*2;++i) fe[i] = 0.0;

  for (i=0;i<numpoint;++i) {
    fe[i*2]   =-4.0*p_stat[i*2]*(1.0-p_stat[i*2]*p_stat[i*2]-p_stat[i*2+1]*p_stat[i*2+1])
      -2.0*p_stat[i*2]*p_stat[i*2+1]*p_stat[i*2+1]/((p_stat[i*2]*p_stat[i*2]+p_stat[i*2+1]*p_stat[i*2+1])*(p_stat[i*2]*p_stat[i*2]+p_stat[i*2+1]*p_stat[i*2+1]));
    
    fe[i*2+1] =-4.0*(1.0-p_stat[i*2]*p_stat[i*2]-p_stat[i*2+1]*p_stat[i*2+1])*p_stat[i*2+1]
      +2.0*p_stat[i*2+1]/(p_stat[i*2]*p_stat[i*2]+p_stat[i*2+1]*p_stat[i*2+1])
      -2.0*p_stat[i*2+1]*p_stat[i*2+1]*p_stat[i*2+1]/((p_stat[i*2]*p_stat[i*2]+p_stat[i*2+1]*p_stat[i*2+1])*(p_stat[i*2]*p_stat[i*2]+p_stat[i*2+1]*p_stat[i*2+1]));

    fe[i*2]   =-4.0*p_stat[i*2]*(1.0-p_stat[i*2]*p_stat[i*2]-p_stat[i*2+1]*p_stat[i*2+1])
      -2.0*p_stat[i*2]*p_stat[i*2+1]*p_stat[i*2+1]/((p_stat[i*2]*p_stat[i*2]+p_stat[i*2+1]*p_stat[i*2+1])*(p_stat[i*2]*p_stat[i*2]+p_stat[i*2+1]*p_stat[i*2+1]));
    
    fe[i*2+1] =-4.0*(1.0-p_stat[i*2]*p_stat[i*2]-p_stat[i*2+1]*p_stat[i*2+1])*p_stat[i*2+1]
      +2.0*p_stat[i*2+1]/(p_stat[i*2]*p_stat[i*2]+p_stat[i*2+1]*p_stat[i*2+1])
      -2.0*p_stat[i*2+1]*p_stat[i*2+1]*p_stat[i*2+1]/((p_stat[i*2]*p_stat[i*2]+p_stat[i*2+1]*p_stat[i*2+1])*(p_stat[i*2]*p_stat[i*2]+p_stat[i*2+1]*p_stat[i*2+1]));


  }
  
  z_string_evolution(p_stat,p_evoluted,fe,numpoint,dt,2);
  z_string_interpolation(p_stat,p_evoluted,numpoint,2);

  /*********************/
  /* free(fe);	       */
  /* free(p_evoluted); */
  /*********************/
}

double z_string_FASYS(double *path,int numpoint,double dt, double *v, double kd[2], double n[2]) {
  int i,j,k;
  double *path_evoluted,**q;
  double **f,*f_path;
  double **fb,**fa,**fd;

  f_path=(double *)gcemalloc(sizeof(double)*5*3*numpoint);
  path_evoluted=(double *)gcemalloc(sizeof(double)*5*3*numpoint);
  q=(double **)gcemalloc(sizeof(double *)*5);
  f=(double **)gcemalloc(sizeof(double *)*5);
  for (i=0;i<5;++i) {
    q[i]=(double *)gcemalloc(sizeof(double)*3);
    f[i]=(double *)gcemalloc(sizeof(double)*3);
  }
  // debug
  /**************************************************/
  /* fb=(double **)gcemalloc(sizeof(double *)*5);   */
  /* fa=(double **)gcemalloc(sizeof(double *)*5);   */
  /* fd=(double **)gcemalloc(sizeof(double *)*5);   */
  /* for (i=0;i<5;++i) {			    */
  /*   fb[i]=(double *)gcemalloc(sizeof(double)*3); */
  /*   fa[i]=(double *)gcemalloc(sizeof(double)*3); */
  /*   fd[i]=(double *)gcemalloc(sizeof(double)*3); */
  /* }						    */
  /**************************************************/
  // debug

  for (i=0;i<numpoint;++i) {
    for (j=0;j<5;++j) for (k=0;k<3;++k) q[j][k]=path[i*5*3+j*3+k];
    z_string_FASYS_calcforce(q,f,kd,n,fb,fa,fd);
    v[i]=z_string_FASYS_calcpote(q,kd,n);
    for (j=0;j<5;++j) for (k=0;k<3;++k) f_path[i*5*3+j*3+k]=-1.0*f[j][k];
  }
  // for debug
  /***********************************/
  /* FILE *db;			     */
  /* db=efopen("db_z_FAYSYS","w");   */
  /* for (j=0;j<numpoint*5*3;++j) {  */
  /*   fprintf(db,"%e\n",f_path[j]); */
  /* }				     */
  /* fclose(db);		     */
  /***********************************/
  // for debug

  z_string_evolution(path,path_evoluted,f_path,numpoint,dt,5*3);
  z_string_interpolation(path,path_evoluted,numpoint,5*3);

}

double z_string_FASYS_calcpote(double **q, double kd[2],double n[2]) {
  int i;
  double pi;
  double v=0.0;
  double kb=100.0,ka=50.0;
  //  double kd[2],n[2];
  double l_eq=1.53,a_eq;
  pi=acos(-1.0);
  a_eq=111.0/180*pi;

  /**************/
  /* kd[0]=2.5; */
  /* kd[1]=2.5; */
  /* n[0]=3.0;  */
  /* n[1]=3.0;  */
  /**************/

  for (i=0;i<4;++i) v+=0.5*kb*(len(q[i],q[i+1])-l_eq)*(len(q[i],q[i+1])-l_eq);
  for (i=0;i<3;++i) v+=0.5*ka*(ang(q[i],q[i+1],q[i+2])-a_eq)*(ang(q[i],q[i+1],q[i+2])-a_eq);
  
  v+=0.5*kd[0]*(1.0+cos(n[0]*(dih(q[0],q[1],q[2],q[3]))))+0.5*kd[1]*(1.0+cos(n[1]*dih(q[1],q[2],q[3],q[4])));

  // for debug
  /*******************************************************************************************************************************/
  /* FILE *db;															 */
  /* db=efopen("db_z_FAYSYS_dp","a");												 */
  /* fprintf(db,"%e\n",0.5*kd[0]*(1.0+cos(n[0]*(dih(q[0],q[1],q[2],q[3]))))+0.5*kd[1]*(1.0+cos(n[1]*dih(q[1],q[2],q[3],q[4])))); */
  /* fclose(db);														 */
  /*******************************************************************************************************************************/
  // for debug

  return v;
}

double z_string_FASYS_calcforce(double **q,double **f, double kd[2], double n[2], double **fb2, double **fa2, double **fd) {
  int i,j,k;
  double pi;
  double f_temp,f_temp1,f_temp2;
  double lenij,lenkj,angijk,angjik,dihedang;
  double k_bond=100.0,k_angle=50.0;
  double len_eq=1.53,ang_eq;

  double vij[3],vkj[3],cosijk;

  double atom_i[3],atom_j[3],atom_k[3],atom_l[3];
  double cosdih=0.0,sindih=0.0;

  double /*kd[2],n[2],*/p[2];
  double n1[3],n2[3],ln1,ln2;
  double vki[3],vjl[3],vji[3],vik[3],vkl[3];
  double fa,fb[3],fc[3]; 
  double op1[3],op2[3],op3[3],op4[3],op5[3],op6[3];

  pi=acos(-1.0);
  ang_eq=111.0/180*pi;

  /**************/
  /* kd[0]=2.5; */
  /* kd[1]=2.5; */
  /* n[0]=3.0;  */
  /* n[1]=3.0;  */
  /**************/
  p[0]=0.0;
  p[1]=0.0;

  for (i=0;i<5;++i) for (j=0;j<3;++j) f[i][j] = 0.0;
  // debug
  for (i=0;i<5;++i) {
    for (j=0;j<3;++j) {
      fb2[i][j] = 0.0;
      fa2[i][j] = 0.0;
      fd[i][j] = 0.0;
    }
  }
  // debug

  // bonding
  for (i=0;i<4;++i) {
    for (j=0;j<3;++j) {
      atom_i[j]=q[i][j];
      atom_j[j]=q[i+1][j];
    }
    lenij = len(atom_i,atom_j);
    for (j=0;j<3;++j) {
      f_temp = -k_bond*(lenij-len_eq)*(q[i+1][j]-q[i][j])/lenij;
      f[i][j] += f_temp;
      f[i+1][j] += -f_temp;
      // debug
      fb2[i][j] += f_temp;
      fb2[i+1][j] += -f_temp;
      // debug
      /************************/
      /* f[i][j] += -f_temp;  */
      /* f[i+1][j] += f_temp; */
      /************************/
    }
  }
  
  // angle
  for (i=0;i<3;++i) {
    lenij = len(q[i],q[i+1]);
    lenkj = len(q[i+2],q[i+1]);
    for (j=0;j<3;++j) {
      vij[j]=q[i+1][j]-q[i][j];
      vkj[j]=q[i+1][j]-q[i+2][j];
    }
    cosijk=inprod(vij,vkj,3);
    cosijk=cosijk/lenij/lenkj;
    angijk = acos(cosijk);
    
    for (j=0;j<3;++j) {
      f_temp1 = -k_angle*(angijk-ang_eq)/(lenij*sin(angijk))*(vkj[j]/lenkj-cosijk*vij[j]/lenij);
      f_temp2 = -k_angle*(angijk-ang_eq)/(lenkj*sin(angijk))*(vij[j]/lenij-cosijk*vkj[j]/lenkj);
      f[i][j] += -f_temp1;
      f[i+2][j] += -f_temp2;
      f[i+1][j] += f_temp1+f_temp2;
      // debug
      fa2[i][j] += f_temp1;
      fa2[i+2][j] += -f_temp2;
      fa2[i+1][j] += -f_temp1+f_temp2;
      // debug
      /**********************************/
      /* f[i][j] += f_temp1;	        */
      /* f[i+2][j] += f_temp2;	        */
      /* f[i+1][j] += -f_temp1-f_temp2; */
      /**********************************/
    }
  }

  // dihed
  for (i=0;i<2;++i) {
    for (j=0;j<3;++j) {
      vij[j] = q[i+1][j]-q[i][j];
      vkj[j] = q[i+1][j]-q[i+2][j];
      vki[j] = q[i][j]-q[i+2][j];
      vjl[j] = q[i+3][j]-q[i+1][j];
      vji[j] = q[i][j]-q[i+1][j];
      vik[j] = q[i+2][j]-q[i][j];
      vkl[j] = q[i+3][j]-q[i+2][j];
    }
    
    outprod(vij,vkj,n1);
    outprod(vkj,vkl,n2);
    ln1=sqrt(inprod(n1,n1,3));
    ln2=sqrt(inprod(n2,n2,3));
    
    csdih(q[i],q[i+1],q[i+2],q[i+3],&cosdih,&sindih);
    
    if (n[i]==1) fa=-n[i]*kd[i]/**cos(p[i])*/;
    else if (n[i]==2) fa=-n[i]*kd[i]*2.0*cosdih/**cos(p[i])*/;
    else if (n[i]==3) fa=-n[i]*kd[i]*(-4.0*sindih*sindih+3.0)/**cos(p[i])*/;
    else if (n[i]==4) fa=-n[i]*kd[i]*4.0*(cosdih*(2.0*cosdih*cosdih-1.0))/**cos(p[i])*/;
    
    for (j=0;j<3;++j) fb[j]=(n2[j]/ln2-cosdih*n1[j]/ln1)/ln1;
    for (j=0;j<3;++j) fc[j]=(n1[j]/ln1-cosdih*n2[j]/ln2)/ln2;
    
    outprod(fb,vkj,op1);
    outprod(fc,vki,op2);
    outprod(fb,vik,op3);
    outprod(fb,vij,op4);
    outprod(fc,vjl,op5);
    outprod(fc,vkj,op6);
    
    for (j=0;j<3;++j) {
      /*************************************/
      /* f[i][j] += fa*op1[j];		   */
      /* f[i+1][j] += fa*(-op2[j]+op3[j]); */
      /* f[i+2][j] += fa*(-op4[j]+op5[j]); */
      /* f[i+3][j] += fa*op6[j];	   */
      /*************************************/
      f[i][j] -= fa*op1[j];
      f[i+1][j] -= fa*(-op2[j]+op3[j]);
      f[i+2][j] -= fa*(-op4[j]+op5[j]);
      f[i+3][j] -= fa*op6[j];
      // debug
      fd[i][j] -= fa*op1[j];
      fd[i+1][j] -= fa*(-op2[j]+op3[j]);
      fd[i+2][j] -= fa*(-op4[j]+op5[j]);
      fd[i+3][j] -= fa*op6[j];
      // debug
    }
    // for debug
    /*****************************************/
    /* FILE *db;			     */
    /* db=efopen("db_z_fa_FASYS","a");	     */
    /* fprintf(db,"%e\n",fa);		     */
    /* fclose(db);			     */
    /* FILE *db2;			     */
    /* db2=efopen("db_z_sd_cd_FASYS","a");   */
    /* fprintf(db2,"%e %e\n",sindih,cosdih); */
    /* fclose(db2);			     */
    /*****************************************/
    // for debug
    // for debug
  }
}

double z_string_Cartesian(double *path,double *path_evoluted,double *fe,int numpoint,int numatom,double dt) {
  z_string_evolution(path,path_evoluted,fe,numpoint,dt,numatom*3);
  z_string_interpolation(path,path_evoluted,numpoint,numatom*3);
}

double z_string_evolution(double *p, double *p_evoluted, double *fe, int numpoint,double dt, int numdimension) {
  int i;
  // for debug
  /*******************************************/
  /* int j;				     */
  /* FILE *db;				     */
  /* db=efopen("db_z","w");		     */
  /* for (j=0;j<numpoint*numdimension;++j) { */
  /*   fprintf(db,"%e\n",fe[j]);	     */
  /* }					     */
  /* fclose(db);			     */
  /*******************************************/
  // for debug

  for (i=0;i<numpoint*numdimension;++i) p_evoluted[i]=p[i]+dt*fe[i];

}

double z_string_interpolation(double *p, double *p_evoluted,int numpoint, int numdimension) {
  int i,j;
  double width,delta;
  double *s,*al,*al_evoluted;
  double *a,*b,*c,*d;
  double *x,*y;
  double *x_interpolated,*y_interpolated;

  for (i=1;i<numpoint*numdimension;++i) p[i]=p_evoluted[i];

  s=(double *)gcemalloc(sizeof(double)*numpoint);
  al=(double *)gcemalloc(sizeof(double)*numpoint);
  al_evoluted=(double *)gcemalloc(sizeof(double)*numpoint);
  
  /****************************************************************/
  /* a=(double *)gcemalloc(sizeof(double)*(numpoint-1));	  */
  /* b=(double *)gcemalloc(sizeof(double)*numpoint);		  */
  /* c=(double *)gcemalloc(sizeof(double)*(numpoint-1));	  */
  /* d=(double *)gcemalloc(sizeof(double)*numpoint);		  */
  /* 								  */
  /* x=(double *)gcemalloc(sizeof(double)*numpoint);		  */
  /* y=(double *)gcemalloc(sizeof(double)*numpoint);		  */
  /* x_interpolated=(double *)gcemalloc(sizeof(double)*numpoint); */
  /* y_interpolated=(double *)gcemalloc(sizeof(double)*numpoint); */
  /****************************************************************/
  
  s[0]=0.0;
  for (i=1;i<numpoint;++i) {
    delta=0.0;
    for (j=0;j<numdimension;++j)
      delta+=(p_evoluted[i*numdimension+j]-p_evoluted[(i-1)*numdimension+j])*(p_evoluted[i*numdimension+j]-p_evoluted[(i-1)*numdimension+j]);
    delta=sqrt(delta);
    s[i]=s[i-1]+delta;
  }
  for (i=0;i<numpoint;++i)  al_evoluted[i]=s[i]/s[numpoint-1];
  for (i=0;i<numpoint;++i)  al[i]=(double)i/(numpoint-1);
  
  for (i=0;i<numdimension;++i) {
    a=(double *)gcemalloc(sizeof(double)*(numpoint-1));
    b=(double *)gcemalloc(sizeof(double)*numpoint);
    c=(double *)gcemalloc(sizeof(double)*(numpoint-1));
    d=(double *)gcemalloc(sizeof(double)*numpoint);
  
    x=(double *)gcemalloc(sizeof(double)*numpoint);
    y=(double *)gcemalloc(sizeof(double)*numpoint);
    x_interpolated=(double *)gcemalloc(sizeof(double)*numpoint);
    y_interpolated=(double *)gcemalloc(sizeof(double)*numpoint);
  
    for (j=0;j<numpoint;++j) {
      x[j]=al_evoluted[j];
      y[j]=p_evoluted[j*numdimension+i];
    }
    width=1.0/(numpoint-1);
    for (j=0;j<numpoint;++j)  x_interpolated[j]=j*width;
  
    csi_calc_abcd(x,y,a,b,c,d,numpoint);
    csi_interpolation2(x,a,b,c,d,x_interpolated,y_interpolated,width,numpoint,numpoint);
  
    /******************************/
    /* for (j=0;j<numpoint;++j) { */
    /*   x_interpolated[j]=x[j];  */
    /*   y_interpolated[j]=y[j];  */
    /* }			  */
    /******************************/

    for (j=0;j<numpoint;++j) {
      al[j]=x_interpolated[j];
      p[j*numdimension+i]=y_interpolated[j];
    }
  }
}


 
