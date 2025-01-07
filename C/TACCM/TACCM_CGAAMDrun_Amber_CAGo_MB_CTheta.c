#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "TACCM_CGAAMDrun_Amber_CAGo_MB_CTheta.h"

#include <math.h>

#include "EF.h"
#include "FFL.h"

#include "TOPO.h"
#include "LA.h"

double TACCM_CTheta_Amber_CAGo_MB(double *crd,int numatom,double *theta, 
				  int numdihe, int **pairs_dih_AA,
				  int numangl, int **pairs_ang_AA,
				  int numbond, int **pairs_bon_AA, 
				  double pi){
  int i,j,k,l;
  int ii,jj,kk,ll;

  double lenij,lenkj;
  double cosijk,angijk;

  double m[3],n[3],m_n[3],n_n[3],lm,ln;
  double vij[3],vkj[3],vkl[3];
  double lkj;
  double vijvkj,vklvkj;

  double atom[4][3];
  double angl,dihed;

  //  printf("here is line 381 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");

  for (i=0;i<numbond;++i) {
    printf("%4d -%4d\n",pairs_bon_AA[i][0],pairs_bon_AA[i][1]);
  }
  
  printf("\n");
  for (i=0;i<numangl;++i) {
    printf("%4d -%4d-%4d\n",pairs_ang_AA[i][0],pairs_ang_AA[i][1],pairs_ang_AA[i][2]);
  }
  
  printf("\n");
  for (i=0;i<numdihe;++i) {
    printf("%4d -%4d -%4d -%4d\n"
  	   ,pairs_dih_AA[i][0],pairs_dih_AA[i][1]
  	   ,pairs_dih_AA[i][2],pairs_dih_AA[i][3]);
  }

  for (i=0;i</*1*/numbond/*-1*/;++i) {
    ii=pairs_bon_AA[i][0]/*-1*/;
    jj=pairs_bon_AA[i][1]/*-1*/;
    for (j=0;j<3;++j) {
      atom[0][j]=crd[ii*3+j];
      atom[1][j]=crd[jj*3+j];
    }
    lenij=0.0;
    for (j=0;j<3;++j) {
      lenij += (atom[0][j]-atom[1][j])*(atom[0][j]-atom[1][j]);
    }
    //    printf("%5.3d %8.5lf\n",i,lenij);
    lenij=sqrt(lenij);
    //    printf("i=%5d lenij=%8.5lf ii=%5d jj=%5d\n",i,lenij,ii,jj);
    theta[i]=lenij;
    /**********************/
    //    printf("%5.3d %8.5lf\n",i,lenij);
  }

  //  printf("here is line 397 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");
  //  printf("numangl=%5.3d\n",numangl);
  ////////////////////////////////////////////////////////////

  for (i=0;i</*1*/numangl/*-1*/;++i) {
    ii=pairs_ang_AA[i][0]/*-1*/;
    jj=pairs_ang_AA[i][1]/*-1*/;
    kk=pairs_ang_AA[i][2]/*-1*/;
    for (j=0;j<3;++j) {
      atom[0][j]=crd[ii*3+j];
      atom[1][j]=crd[jj*3+j];
      atom[2][j]=crd[kk*3+j];
    }
  
    lenij = len(atom[0],atom[1]);
    lenkj = len(atom[2],atom[1]);
    for (j=0;j<3;++j) {
      vij[j]=atom[1][j]-atom[0][j];
      vkj[j]=atom[1][j]-atom[2][j];
    }
    cosijk=inprod(vij,vkj,3);
    cosijk=cosijk/lenij/lenkj;
    if (cosijk>=1.0) angijk=0.0;
    else if (cosijk<=0.0) angijk=pi;
    else  angijk = acos(cosijk);
  
    theta[i+numbond]=angijk;
  }

  //  printf("here is line 426 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");
  ////////////////////////////////////////////////////////////

  /******************************************************************/
  /* for (i=0;i<numdihe;++i) {					    */
  /*   ii=pairs_dih_AA[i][0]/\*-1*\/;				    */
  /*   jj=pairs_dih_AA[i][1]/\*-1*\/;				    */
  /*   kk=pairs_dih_AA[i][2]/\*-1*\/;				    */
  /*   ll=pairs_dih_AA[i][3]/\*-1*\/;				    */
  /* 								    */
  /*   for (j=0;j<3;++j) {					    */
  /*     atom[0][j]=crd[ii*3+j];				    */
  /*     atom[1][j]=crd[jj*3+j];				    */
  /*     atom[2][j]=crd[kk*3+j];				    */
  /*     atom[3][j]=crd[ll*3+j];				    */
  /*   }							    */
  /* 								    */
  /*   for (j=0;j<3;++j) {					    */
  /*     vij[j] = atom[1][j]-atom[0][j];			    */
  /*     vkj[j] = atom[1][j]-atom[2][j];			    */
  /*     vkl[j] = atom[3][j]-atom[2][j];			    */
  /*   }							    */
  /*   //    lkj=sqrt(inprod(vkj,vkj,3));			    */
  /* 								    */
  /* /\*   outprod(vij,vkj,m);		      *\/		    */
  /* /\*   outprod(vkj,vkl,n);		      *\/		    */
  /* /\*   lm=sqrt(inprod(m,m,3));		      *\/	    */
  /* /\*   ln=sqrt(inprod(n,n,3));		      *\/	    */
  /* /\*   for (j=0;j<3;++j) {		      *\/		    */
  /* /\*     m_n[j]=m[j]/lm;		      *\/		    */
  /* /\*     n_n[j]=n[j]/ln;		      *\/		    */
  /* /\*   }				      *\/		    */
  /* /\* 					      *\/	    */
  /* /\*   dihed=inprod(m_n,n_n,3);		      *\/	    */
  /* /\*   if (dihed>=1.0)			      *\/	    */
  /* /\*     dihed=0.0;			      *\/		    */
  /* /\*   else if (dihed<=-1.0)		      *\/	    */
  /* /\*     dihed=pi;			      *\/		    */
  /* /\*   else				      *\/		    */
  /* /\*     dihed=acos(dihed);		      *\/		    */
  /* /\*   if (inprod(vij,n,3)>0) dihed=-dihed;   *\/		    */
  /* /\*   if (dihed<-1.0*pi) dihed=2.0*pi+dihed; *\/		    */
  /* /\*   if (dihed>pi) dihed=-2.0*pi+dihed;     *\/		    */
  /* /\* 					      *\/	    */
  /* /\*   theta[i+numbond+numangl]=dihed;	      *\/	    */
  /* }								    */
  /******************************************************************/

  //  printf("here is line 471 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");

  for (i=0;i<numdihe;++i) {
    ii=pairs_dih_AA[i][0]/*-1*/;
    jj=pairs_dih_AA[i][1]/*-1*/;
    kk=pairs_dih_AA[i][2]/*-1*/;
    ll=pairs_dih_AA[i][3]/*-1*/;

    for (j=0;j<3;++j) {
      atom[0][j]=crd[ii*3+j];
      atom[1][j]=crd[jj*3+j];
      atom[2][j]=crd[kk*3+j];
      atom[3][j]=crd[ll*3+j];
    }

    for (j=0;j<3;++j) {
      vij[j] = atom[1][j]-atom[0][j];
      vkj[j] = atom[1][j]-atom[2][j];
      vkl[j] = atom[3][j]-atom[2][j];
    }
    lkj=sqrt(inprod(vkj,vkj,3));
    
    outprod(vij,vkj,m);
    outprod(vkj,vkl,n);
    lm=sqrt(inprod(m,m,3));
    ln=sqrt(inprod(n,n,3));
    for (j=0;j<3;++j) {
      m_n[j]=m[j]/lm;
      n_n[j]=n[j]/ln;
    }
    
    dihed=inprod(m_n,n_n,3);
    if (dihed>=1.0)
      dihed=0.0;
    else if (dihed<=-1.0)
      dihed=pi;
    else
      dihed=acos(dihed);
    if (inprod(vij,n,3)>0) dihed=-dihed;
    if (dihed<-1.0*pi) dihed=2.0*pi+dihed;
    if (dihed>pi) dihed=-2.0*pi+dihed;

    theta[i+numbond+numangl]=dihed;
  }

  return 0.0;
}
