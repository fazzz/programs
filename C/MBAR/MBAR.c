
#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "MBAR.h"
#include "HIST.h"
#include "mymath.h"
#include "LA.h"

double BAR_ite(double F,double *ene_U0_C0,double *ene_U0_C1, double *ene_U1_C0,double *ene_U1_C1, int n0, int n1, double criteria_BAR){
  int i,j,ite;
  double C=0.0,sum0,sum1,F_n,dF;

  for (ite=0;dF>=criteria_BAR;++ite) {
    for (j = 0; j < n0; ++j) {
      sum0+=Fermi(ene_U0_C0[i]-ene_U0_C1[i]+C);
    }
    for (j = 0; j < n1; ++j) {
      sum1+=Fermi(ene_U1_C0[i]-ene_U1_C1[i]+C);
    }
    
    F_n=log(sum0/sum1)+C-log(n1/n0);
    C=F+log(n1/n0);
    dF=fabs(F_n-F);
    F=F_n;
  }
  return dF;
}

double MBAR_ite(double *F, double ***ene, int n_sim, int *n, double criteria_BAR, int MAXITE ){
  int i,j,k,ite,nn;
  double *F_n,*dF;
  double num,din,max_dF=1.0;
  FILE *logfile;

  F_n  = (double *)emalloc(sizeof(double)*n_sim);
  dF   = (double *)emalloc(sizeof(double)*n_sim);

  for (i = 0; i < n_sim; ++i) {
    F[i]=0.0;   	
  }


  for (ite = 0; max_dF > criteria_BAR && ite < MAXITE ; ++ite) {
    for (i = 0; i < n_sim; ++i) {
      num=0.0;
      for (j = 0; j < n_sim; ++j) {
	for (nn = 0; nn < n[j]; ++nn) {
	  din=0.0;
	  for (k = 0; k < n_sim; ++k) {	  
	    din+=n[k]*exp(F[k]-ene[k][j][nn]);
	  }
	  num+=exp(-ene[i][j][nn])/din;
	}
      }
      F_n[i]=-log(num);
    }

    max_dF=F_n[0]-F[0];
    if (max_dF<0.0) max_dF=-max_dF;
    for (i = 1; i < n_sim; ++i) {
      dF[i]=F_n[i]-F[i];
      if (dF[i]<0.0) dF[i]=-dF[i];
      if (max_dF<dF[i]) max_dF=dF[i];
      F[i]=F_n[i];
    }

    //    logfile=efopen("logMBAR.txt","a");
    //    fprintf(logfile,"ite=%d \n# fene \n",ite);
    //    for (i=0;i<n_sim;++i)
    //      fprintf(logfile,"%d %10.4lf\n",i,F[i]);
    //    fclose(logfile);
    
    logfile=efopen("logMBAR.txt","a");
    fprintf(logfile,"ite=%d %12.9e\n",ite,max_dF);
    fclose(logfile);
  }
    
  return 0.0;
}

double MBAR_ite_high_speed(double *F, double ***expene, int n_sim, int *n, double criteria_BAR, int MAXITE ){
  int i,j,k,l,m,ite,nn;
  double *expF,*expF_n,*dF;
  double num,din,max_dF=1.0;
  FILE *logfile;

  expF  = (double *)gcemalloc(sizeof(double)*n_sim);
  expF_n  = (double *)gcemalloc(sizeof(double)*n_sim);
  dF   = (double *)gcemalloc(sizeof(double)*n_sim);

  for (i = 0; i < n_sim; ++i) {
    F[i]=0.0;   	
    expF[i]=1.0;   	
  }

  //  logfile=efopen("logMBAR.txt","w");
  for (ite = 0; max_dF > criteria_BAR && ite < MAXITE ; ++ite) {
    for (i = 0; i < n_sim; ++i) {
      num=0.0;
      for (j = 0; j < n_sim; ++j) {
	for (nn = 0; nn < n[j]; ++nn) {
	  din=0.0;
	  for (k = 0; k < n_sim; ++k) {	  
	    din+=n[k]*expene[k][j][nn]*expF[k];
	  }
	  if (expene[i][j][nn]!=0.0)  {  // 11-12-26
	    num+=expene[i][j][nn]/din;
	    if (din==0.0)  {
	      printf("error din=0.0\n");
	    }
	  }
	}
      }
      expF_n[i]=1.0/num;
    }

    max_dF=log(expF_n[0]/expF[0]);
    if (max_dF<0.0) max_dF=-max_dF;
    for (i = 1; i < n_sim; ++i) {
      dF[i]=log(expF_n[i]/expF[i]);
      if (dF[i]<0.0) dF[i]=-dF[i];
      if (max_dF<dF[i]) max_dF=dF[i];
      expF[i]=expF_n[i];
    }

    //    logfile=efopen("logMBAR.txt","a");
    //    fprintf(logfile,"ite=%d \n# fene \n",ite);
    //    for (i=0;i<n_sim;++i)
    //      fprintf(logfile,"%d %10.4lf\n",i,F[i]);
    //    fclose(logfile);
    
    logfile=efopen("logMBAR.txt","a");
    fprintf(logfile,"ite=%d %12.9e\n",ite,max_dF);
    fclose(logfile);
  }
  //  fclose(logfile);

  for (i=0;i<n_sim;++i)
    F[i]=log(expF[i]);

    
  return 0.0;
}

double MBAR_ite_high_speed_2(double *expF, double ***expene, int n_sim, int *n, double criteria_BAR, int MAXITE ){
  int i,j,k,l,m,ite,nn;
  double *expF_n,*dF;
  double num,din,max_dF=1.0;
  FILE *logfile;

  expF_n  = (double *)gcemalloc(sizeof(double)*n_sim);
  dF   = (double *)gcemalloc(sizeof(double)*n_sim);

  for (i = 0; i < n_sim; ++i) {
    expF[i]=1.0;   	
  }

  //  logfile=efopen("logMBAR.txt","w");
  for (ite = 0; max_dF > criteria_BAR && ite < MAXITE ; ++ite) {
    for (i = 0; i < n_sim; ++i) {
      num=0.0;
      for (j = 0; j < n_sim; ++j) {
	for (nn = 0; nn < n[j]; ++nn) {
	  din=0.0;
	  for (k = 0; k < n_sim; ++k) {	  
	    din+=n[k]*expene[k][j][nn]*expF[k];
	  }
	  num+=expene[i][j][nn]/din;
	  if (din==0.0)  {
	    printf("error din=0.0\n");
	  }
	}
      }
      expF_n[i]=1.0/num;
    }

    max_dF=log(expF_n[0]/expF[0]);
    if (max_dF<0.0) max_dF=-max_dF;
    for (i = 1; i < n_sim; ++i) {
      dF[i]=log(expF_n[i]/expF[i]);
      if (dF[i]<0.0) dF[i]=-dF[i];
      if (max_dF<dF[i]) max_dF=dF[i];
      expF[i]=expF_n[i];
    }
    
    logfile=efopen("logMBAR.txt","a");
    fprintf(logfile,"ite=%d %12.9e\n",ite,max_dF);
    fclose(logfile);
  }

  return 0.0;
}


double *MBAR_AVE_oned(double *fene, double ***enek, int n_sim, int *n, double **od_data, double width,double *max,double *min, int *frame){
  int i,j,k,t;
  
  double *ef;
  double ***euk;

  double *hist;

  int n_total;
  double din;
  double normz_const=0.0;
  double sum;

  ef=(double *)gcemalloc(sizeof(double)*n_sim);
  euk=(double ***)gcemalloc(sizeof(double **)*n_sim);
  for (i=0;i<n_sim;++i) 
    euk[i]=(double **)gcemalloc(sizeof(double *)*n_sim);
  for (i=0;i<n_sim;++i) 
    for (j=0;j<n_sim;++j) 
      euk[i][j]=(double *)gcemalloc(sizeof(double)*n[j]);

  for (i=0;i<n_sim;++i) 
    for (j=0;j<n_sim;++j) 
      for (t=0;t<n[i];++t)
	euk[i][j][t]=exp(-1.0*enek[i][j][t]);

  for (i=0;i<n_sim;++i) ef[i]=exp(fene[i]);

  n_total=0; for (i=0;i<n_sim;++i)  n_total+=n[i];

  normz_const=0.0;
  for (i=0;i<n_sim;++i) {
    for (t=0;t<n[i];++t) {
      din=0.0;
      for (j=0;j<n_sim;++j)
	din+=n[j]*ef[j]*euk[j][i][t];
      normz_const+=1.0/din;
    }
  }

  *max=od_data[0][0];*min=od_data[0][0];
  for (i=0;i<n_sim;++i) {
    for (j=0;j<n[i];++j) {
      if (*max < od_data[i][j]) *max=od_data[i][j];
      if (*min > od_data[i][j]) *min=od_data[i][j];
    }
  }
  *frame=(int)((*max-*min)/width)+1;
  hist=(double *)gcemalloc(sizeof(double)*(*frame));
  for (i=0;i<*frame;++i) hist[i]=0.0;

  for (i=0;i<n_sim;++i) {
    for (t=0;t<n[i];++t) {
      din=0.0;
      for (k=0;k<n_sim;++k)
	din+=n[k]*ef[k]*euk[k][i][t];
      hist[((int)((od_data[i][t]-*min)/width))]+=1.0/din/normz_const;
    }
  }

  sum=0.0;
  for (i=0;i<*frame;++i)
   sum+=hist[i];

  for (i=0;i<*frame;++i)
    hist[i]/=sum;

  return hist;
}

double *MBAR_AVE_oned_multi_temp(double *fene, double ***enek,double *betaene, int n_sim, int *n, double **od_data, double width,double *max,double *min, int *frame){
  int i,j,k,t,t_total;
  
  double *ef;
  double ***euk;

  double *hist;

  int n_total;
  double din;
  double expbetaene;
  double normz_const=0.0;
  double sum;

  ef=(double *)gcemalloc(sizeof(double)*n_sim);
  euk=(double ***)gcemalloc(sizeof(double **)*n_sim);
  for (i=0;i<n_sim;++i) 
    euk[i]=(double **)gcemalloc(sizeof(double *)*n_sim);
  for (i=0;i<n_sim;++i) 
    for (j=0;j<n_sim;++j) 
      euk[i][j]=(double *)gcemalloc(sizeof(double)*n[j]);

  for (i=0;i<n_sim;++i) 
    for (j=0;j<n_sim;++j) 
      for (t=0;t<n[i];++t)
	euk[i][j][t]=enek[i][j][t];                   // 11-12-26
  //	euk[i][j][t]=exp(-1.0*enek[i][j][t]); // 11-12-26

  for (i=0;i<n_sim;++i) ef[i]=exp(fene[i]);

  n_total=0; for (i=0;i<n_sim;++i)  n_total+=n[i];

  normz_const=0.0;
  t_total=0;
  for (i=0;i<n_sim;++i) {
    for (t=0;t<n[i];++t) {
      din=0.0;
      for (j=0;j<n_sim;++j)
	din+=n[j]*ef[j]*euk[j][i][t];
      expbetaene=exp(-1.0*betaene[t_total]);
      normz_const+=1.0/din*expbetaene;
      ++t_total;
    }
  }

  *max=od_data[0][0];*min=od_data[0][0];
  for (i=0;i<n_sim;++i) {
    for (j=0;j<n[i];++j) {
      if (*max < od_data[i][j]) *max=od_data[i][j];
      if (*min > od_data[i][j]) *min=od_data[i][j];
    }
  }
  *frame=(int)round((*max-*min)/width)+1;
  hist=(double *)gcemalloc(sizeof(double)*(*frame));
  for (i=0;i<*frame;++i) hist[i]=0.0;

  t_total=0;
  for (i=0;i<n_sim;++i) {
    for (t=0;t<n[i];++t) {
      din=0.0;
      for (k=0;k<n_sim;++k)
	din+=n[k]*ef[k]*euk[k][i][t];
      expbetaene=exp(-1.0*betaene[t_total]);
      hist[((int)round((od_data[i][t]-*min)/width))]+=1.0/din/normz_const*expbetaene;
      ++t_total;
    }
  }

  sum=0.0;
  for (i=0;i<*frame;++i)
   sum+=hist[i];

  for (i=0;i<*frame;++i)
    hist[i]/=sum;

  return hist;
}

double *MBAR_AVE_oned_multi_temp_2(double *ef, double ***euk, int n_sim, int *n, double **od_data, double width,double *max,double *min, int *frame, int numk){
  int i,j,k,t,t_total;
  
  double *hist;

  int n_total;
  double din;
  double expbetaene;
  double normz_const=0.0;
  double sum;

  n_total=0; for (i=0;i<n_sim;++i)  n_total+=n[i];

  *max=od_data[0][0];*min=od_data[0][0];
  for (i=0;i<n_sim;++i) {
    for (j=0;j<n[i];++j) {
      if (*max < od_data[i][j]) *max=od_data[i][j];
      if (*min > od_data[i][j]) *min=od_data[i][j];
    }
  }
  *frame=(int)round((*max-*min)/width)+1;
  hist=(double *)gcemalloc(sizeof(double)*(*frame));
  for (i=0;i<*frame;++i) hist[i]=0.0;

  for (i=0;i<n_sim;++i) {
    for (t=0;t<n[i];++t) {
      din=0.0; for (j=0;j<n_sim;++j) din+=n[j]*ef[j]*euk[j][i][t];
      hist[((int)round((od_data[i][t]-*min)/width))]+=1.0/din*ef[numk];
    }
  }

  sum=0.0;
  for (i=0;i<*frame;++i) sum+=hist[i];

  for (i=0;i<*frame;++i) hist[i]/=sum;

  return hist;
}

double MBAR_AVEV_multi_temp(double *ef, double ***euk, int n_sim, int *n, double **od_data,int numk){
  int i,j,t;
  
  double AVEV=0.0;
  double din;

  for (i=0;i<n_sim;++i) {
    for (t=0;t<n[i];++t) {
      din=0.0; for (j=0;j<n_sim;++j) din+=n[j]*ef[j]*euk[j][i][t];
      AVEV+=1.0/din*od_data[i][t]*ef[numk]/**euk[i][numk][t]*/;
    }
  }

  return AVEV;
}

double MBAR_AVE_value_n_multi_temp(double *fene, double ***enek,double *betaene, int n_sim, int *n, double **data, double n_pow){
  int i,j,k,t,t_total;
  
  double *ef;
  double ***euk;

  double ave;

  int n_total;
  double din;
  double expbetaene;
  double normz_const=0.0;

  ef=(double *)gcemalloc(sizeof(double)*n_sim);
  euk=(double ***)gcemalloc(sizeof(double **)*n_sim);
  for (i=0;i<n_sim;++i) 
    euk[i]=(double **)gcemalloc(sizeof(double *)*n_sim);
  for (i=0;i<n_sim;++i) 
    for (j=0;j<n_sim;++j) 
      euk[i][j]=(double *)gcemalloc(sizeof(double)*n[j]);

  for (i=0;i<n_sim;++i) 
    for (j=0;j<n_sim;++j) 
      for (t=0;t<n[i];++t)
	euk[i][j][t]=exp(-1.0*enek[i][j][t]);

  for (i=0;i<n_sim;++i) ef[i]=exp(fene[i]);

  n_total=0; for (i=0;i<n_sim;++i)  n_total+=n[i];

  normz_const=0.0;
  t_total=0;
  for (i=0;i<n_sim;++i) {
    for (t=0;t<n[i];++t) {
      din=0.0;
      for (j=0;j<n_sim;++j)
	din+=n[j]*ef[j]*euk[j][i][t];
      expbetaene=exp(-1.0*betaene[t_total]);
      normz_const+=1.0/din*expbetaene;
      ++t_total;
    }
  }

  ave=0.0;
  t_total=0;
  for (i=0;i<n_sim;++i) {
    for (t=0;t<n[i];++t) {
      din=0.0;
      for (k=0;k<n_sim;++k)
	din+=n[k]*ef[k]*euk[k][i][t];
      expbetaene=exp(-1.0*betaene[t_total]);
      ave+=pow(data[i][t],n_pow)/din/normz_const*expbetaene;
      ++t_total;
    }
  }

  return ave;
}

double *MBAR_ACM0(int n_sim, int *n, double *W, double *WT){
  int i,j,k,l,t;
  int n_total=0;
  double din;
  double **ene;
  double *N;
  double *A,*AP;

  double *WTAP,*WTAPW;
  double *WN,*WNWT;

  double *sigma,*Usigma,*B;
  double *AAP,*AAPA;

  for (i=0;i<n_sim;++i) n_total+=n[i];

  N=gcemalloc(sizeof(double)*n_sim*n_sim);
  msetzero(N,n_sim);
  for (i=0;i<n_sim;++i) N[i*n_sim+i]=n[i];

  WN=(double *)gcemalloc(sizeof(double)*n_total*n_sim);
  WNWT=(double *)gcemalloc(sizeof(double)*n_total*n_total);

  mnmult(W,n_total,n_sim,N,n_sim,n_sim,WN);
  mnmult(WN,n_total,n_sim,WT,n_sim,n_total,WNWT);

  A=gcemalloc(sizeof(double)*n_total*n_total);
  AP=gcemalloc(sizeof(double)*n_total*n_total);

  for (i=0;i<n_total;++i) {
    for (j=0;j<n_total;++j) {
      if (i==j)
	A[i*n_total+j]=1.0-WNWT[i*n_total+j];
      else
	A[i*n_total+j]=-WNWT[i*n_total+j];
    }
  }

  MPginvm2(A,AP,n_total,n_total);

  // check for GIM calc
  /***************************************************/
  /* AAP=gcemalloc(sizeof(double)*n_total*n_total);  */
  /* AAPA=gcemalloc(sizeof(double)*n_total*n_total); */
  /* 						     */
  /* LA_mmult(A,AP,AAP,n_total);		     */
  /* LA_mmult(AAP,A,AAPA,n_total);		     */
  /***************************************************/

  WTAP=gcemalloc(sizeof(double)*n_sim*n_total);
  WTAPW=gcemalloc(sizeof(double)*n_sim*n_sim);

  mnmult(WT,n_sim,n_total,AP,n_total,n_total,WTAP);
  mnmult(WTAP,n_sim,n_total,W,n_total,n_sim,WTAPW);

  return WTAPW;

}

double *MBAR_ACM(int n_sim, int *n, double *W){
  int i,j,k,l,t;
  int n_total=0;
  double din;
  double **ene;
  double *N;
  double *U,*VT,*V,*sv,*S;
  double *VTN,*VTNV,*SVTNV,*SVTNVS;
  double *A,*AP;
  double *SAP,*SAPS;
  double *VSAPS,*VSAPSVT;
  double *c1,*c2;

  double *sigma,*Usigma,*B;
  double *AAP,*AAPA;

  for (i=0;i<n_sim;++i) n_total+=n[i];

  N=gcemalloc(sizeof(double)*n_sim*n_sim);
  msetzero(N,n_sim);
  for (i=0;i<n_sim;++i) N[i*n_sim+i]=n[i];

  U=gcemalloc(sizeof(double)*n_total*n_total);
  VT=gcemalloc(sizeof(double)*n_sim*n_sim);
  V=gcemalloc(sizeof(double)*n_sim*n_sim);
  sv=gcemalloc(sizeof(double)*n_sim);
  svd(W,n_total,n_sim,U,VT,sv);

  sigma=gcemalloc(sizeof(double)*n_total*n_sim);
  for (i=0;i<n_total*n_sim;++i) sigma[i]=0.0;
  for (i=0;i<n_sim;++i) {
    if (sv[i]!=0.0)
      sigma[i*n_sim+i]=sv[i];
  }
  Usigma=gcemalloc(sizeof(double)*n_total*n_sim);  
  B=gcemalloc(sizeof(double)*n_total*n_total);     

  mnmult(U,n_total,n_total,sigma,n_total,n_sim,Usigma);
  mnmult(Usigma,n_total,n_sim,VT,n_sim,n_sim,B);
 
  S=gcemalloc(sizeof(double)*n_sim*n_sim);
  msetzero(S,n_sim);
  for (i=0;i<n_sim;++i) S[i*n_sim+i]=sv[i];

  VTN=gcemalloc(sizeof(double)*n_sim*n_sim);
  VTNV=gcemalloc(sizeof(double)*n_sim*n_sim);

  mtrans(VT,V,n_sim);
  LA_mmult(VT,N,VTN,n_sim);
  LA_mmult(VTN,V,VTNV,n_sim);

  SVTNV=gcemalloc(sizeof(double)*n_sim*n_sim);
  SVTNVS=gcemalloc(sizeof(double)*n_sim*n_sim);

  LA_mmult(S,VTNV,SVTNV,n_sim);
  LA_mmult(SVTNV,S,SVTNVS,n_sim);

  A=gcemalloc(sizeof(double)*n_sim*n_sim);
  AP=gcemalloc(sizeof(double)*n_sim*n_sim);

  for (i=0;i<n_sim;++i) {
    for (j=0;j<n_sim;++j) {
      if (i==j)
	A[i*n_sim+j]=1.0-SVTNVS[i*n_sim+j];
      else
	A[i*n_sim+j]=-SVTNVS[i*n_sim+j];
    }
  }

  //  MPginvm(A,AP,n_sim,n_sim);
  MPginvm2(A,AP,n_sim,n_sim);

  SAP=gcemalloc(sizeof(double)*n_sim*n_sim);
  SAPS=gcemalloc(sizeof(double)*n_sim*n_sim);

  LA_mmult(S,AP,SAP,n_sim);
  LA_mmult(SAP,S,SAPS,n_sim);

  VSAPS=gcemalloc(sizeof(double)*n_sim*n_sim);
  VSAPSVT=gcemalloc(sizeof(double)*n_sim*n_sim);

  LA_mmult(V,SAPS,VSAPS,n_sim);
  LA_mmult(VSAPS,VT,VSAPSVT,n_sim);

  return VSAPSVT;

}

double *MBAR_ACM2(int n_sim, int *n, double *W, double *WT){
  int i,j,k,l,t;
  int n_total=0;
  double din;
  double *WTW,*WTWI,*N;
  double *A,*AI;

  for (i=0;i<n_sim;++i) n_total+=n[i];

  WTW=gcemalloc(sizeof(double)*n_sim*n_sim);
  WTWI=gcemalloc(sizeof(double)*n_sim*n_sim);  

  mnmult(WT,n_sim,n_total,W,n_total,n_sim,WTW);

  N=gcemalloc(sizeof(double)*n_sim*n_sim);
  msetzero(N,n_sim);
  for (i=0;i<n_sim;++i) N[i*n_sim+i]=n[i];

  A=gcemalloc(sizeof(double)*n_sim*n_sim);
  AI=gcemalloc(sizeof(double)*n_sim*n_sim);

  invm(WTW,WTWI,n_sim);

  for (i=0;i<n_sim;++i)
    for (j=0;j<n_sim;++j)
      A[i*n_sim+j]=WTWI[i*n_sim+j]-N[i*n_sim+j]+1.0/n_total;

  invm(A,AI,n_sim);

  return AI;

}

double *MBAR_ACM_histod(double *fene, double ***enek, int n_sim, int *n, double **od_data, double width,double *max,double *min, int *frame, double *W){
  int i,j,k,l,t,num;
  int j2,j3,t2;
  int n_total=0;
  double din,normz_const;
  double **ene;
  double *N;
  double *U,*VT,*V,*sv,*S;
  double *VTN,*VTNV,*SVTNV,*SVTNVS;
  double *A,*AP;
  double *SAP,*SAPS;
  double *VSAPS,*VSAPSVT;

  double sum,CA;
  double *hist_CA,*hist_error;

  for (i=0;i<n_sim;++i) n_total+=n[i];
  ene=gcemalloc(sizeof(double *)*n_sim);
  for (i=0;i<n_sim;++i)
    ene[i]=gcemalloc(sizeof(double)*n_total);

  for (i=0;i<n_sim;++i) {
    k=0;
    for (j=0;j<n_sim;++j) {
      for (t=0;t<n[i];++t) {
	ene[i][k]=enek[i][j][t];
	++k;
      }
    }
  }

  normz_const=0.0;
  for (i=0;i<n_sim;++i) {
    for (t=0;t<n[i];++t) {
      din=0.0;
      for (j=0;j<n_sim;++j)
	din+=n[j]*exp(fene[j]-enek[j][i][t]);
      normz_const+=1.0/din;
    }
  }

  for (j=0;j<n_sim;++j) {
    j2=-1;
    for (t=0;t<n[j];++t) {
      ++j2;
      din=0.0;
      for (k=0;k<n_sim;++k) {
	din+=n[k]*exp(fene[k]-ene[k][j]);
      }
      W[j2*(n_sim+2)+n_sim]=1.0/normz_const/din;
    }
  }

  *max=od_data[0][0];*min=od_data[0][0];
  for (i=0;i<n_sim;++i) {
    for (j=0;j<n[i];++j) {
      if (*max < od_data[i][j]) *max=od_data[i][j];
      if (*min > od_data[i][j]) *min=od_data[i][j];
    }
  }
  *frame=(int)((*max-*min)/width)+1;
  hist_CA=(double *)gcemalloc(sizeof(double)*(*frame));
  hist_error=(double *)gcemalloc(sizeof(double)*(*frame));
  for (i=0;i<*frame;++i) hist_CA[i]=0.0;
  for (i=0;i<*frame;++i) hist_error[i]=0.0;

  for (i=0;i<n_sim;++i) {
    for (t=0;t<n[i];++t) {
      din=0.0;
      hist_CA[((int)((od_data[i][t]-*min)/width))]+=1.0;
    }
  }

  for (i=0;i<n_sim;++i) {
    for (t=0;t<n[i];++t) {

      j3=-1;
      CA=0.0;
      for (j2=0;j2<n_sim;++j2) {
	for (t2=0;t2<n[j2];++t2) {
	  ++j3;
	  din=0.0;
	  for (k=0;k<n_sim;++k) {
	    din+=n[k]*exp(fene[k]-enek[k][j2][t2]);
	  }
	  CA+=1.0/din*hist_CA[((int)((od_data[j2][t2]-*min)/width))];
	}
      }
      for (j2=0;j2<n_sim;++j2) {
	for (t2=0;t2<n[j2];++t2) {
	  ++j3;
	  din=0.0;
	  for (k=0;k<n_sim;++k) {
	    din+=n[k]*exp(fene[k]-enek[k][j2][t2]);
	  }
	  W[j3*(n_sim+2)+n_sim+1]=1.0/din/CA*hist_CA[((int)((od_data[i][t]-*min)/width))];
	}
      }

      N=gcemalloc(sizeof(double)*(n_sim+2)*(n_sim+2));
      msetzero(N,(n_sim+2));
      for (k=0;k<n_sim;++k) N[k*((n_sim+2))+k]=n[k];
    
      U=gcemalloc(sizeof(double)*n_total*n_total);
      VT=gcemalloc(sizeof(double)*(n_sim+2)*(n_sim+2));
      V=gcemalloc(sizeof(double)*(n_sim+2)*(n_sim+2));
      sv=gcemalloc(sizeof(double)*(n_sim+2));
      svd(W,n_total,(n_sim+2),U,VT,sv);

      S=gcemalloc(sizeof(double)*(n_sim+2)*(n_sim+2));
      msetzero(S,(n_sim+2));
      for (k=0;k<n_sim;++k) S[k*((n_sim+2))+k]=sv[k];
      
      VTN=gcemalloc(sizeof(double)*(n_sim+2)*(n_sim+2));
      VTNV=gcemalloc(sizeof(double)*(n_sim+2)*(n_sim+2));
      
      mtrans(VT,V,(n_sim+2));
      LA_mmult(VT,N,VTN,(n_sim+2));
      LA_mmult(VTN,V,VTNV,(n_sim+2));
      
      SVTNV=gcemalloc(sizeof(double)*(n_sim+2)*(n_sim+2));
      SVTNVS=gcemalloc(sizeof(double)*(n_sim+2)*(n_sim+2));

      LA_mmult(S,VTNV,SVTNV,(n_sim+2));
      LA_mmult(SVTNV,S,SVTNVS,(n_sim+2));

      A=gcemalloc(sizeof(double)*(n_sim+2)*(n_sim+2));
      AP=gcemalloc(sizeof(double)*(n_sim+2)*(n_sim+2));
      
      for (k=0;k<(n_sim+2);++k) {
	for (l=0;l<(n_sim+2);++l) {
	  if (k==l)
	    A[k*(n_sim+2)+l]=1.0-SVTNVS[k*(n_sim+2)+l];
	  else
	    A[k*(n_sim+2)+l]=-SVTNVS[k*(n_sim+2)+l];
	}
      }
	
      MPginvm(A,AP,(n_sim+2),(n_sim+2));

      SAP=gcemalloc(sizeof(double)*(n_sim+2)*(n_sim+2));
      SAPS=gcemalloc(sizeof(double)*(n_sim+2)*(n_sim+2));
	
      LA_mmult(S,AP,SAP,(n_sim+2));
      LA_mmult(SAP,S,SAPS,(n_sim+2));
	
      VSAPS=gcemalloc(sizeof(double)*(n_sim+2)*(n_sim+2));
      VSAPSVT=gcemalloc(sizeof(double)*(n_sim+2)*(n_sim+2));
    
      LA_mmult(V,SAPS,VSAPS,(n_sim+2));
      LA_mmult(VSAPS,VT,VSAPSVT,(n_sim+2));

      hist_error[((int)((od_data[i][t]-*min)/width))]=
	(hist_CA[((int)((od_data[i][t]-*min)/width))]/normz_const)
	*(hist_CA[((int)((od_data[i][t]-*min)/width))]/normz_const)
	*(VSAPSVT[(n_sim+1)*((n_sim+2))+n_sim+1]
	  +VSAPSVT[(n_sim)*((n_sim+2))+n_sim]
	  -2.0*VSAPSVT[(n_sim+1)*((n_sim+2))+n_sim]);
    }
  }

  return hist_error;

}

double *MBAR_ACM2_histod(double *fene, double ***enek, int n_sim, int *n, double **od_data, double width,double *max,double *min, int *frame, double *W, double *WT){
  int i,j,k,l,t,num,nn;
  int j2,j3,t2;
  int n_total=0;
  int *N;
  double din,normz_const,normz_const_for_hg;
  double *c1,*c2;
  double **ene;
  double *WTW,*WTWI;
  double *A,*AI;
  double *covm;

  double sum;
  double *hist_CA,*hist_error;

  for (i=0;i<n_sim;++i) n_total+=n[i];
  ene=gcemalloc(sizeof(double *)*n_sim);
  for (i=0;i<n_sim;++i) ene[i]=gcemalloc(sizeof(double)*n_total);

  for (i=0;i<n_sim;++i) {
    k=0;
    for (j=0;j<n_sim;++j) {
      for (t=0;t<n[i];++t) {
	ene[i][k]=enek[i][j][t];
	++k;
      }
    }
  }

  normz_const=0.0;
  for (i=0;i<n_sim;++i) {
    for (t=0;t<n[i];++t) {
      din=0.0;
      for (j=0;j<n_sim;++j)
	din+=n[j]*exp(fene[j]-enek[j][i][t]);
      normz_const+=1.0/din;
    }
  }

  num=n_sim;
  j=-1;
  for (i=0;i<n_sim;++i) {
    for (t=0;t<n[i];++t) {
      ++j;
      din=0.0;
      for (k=0;k<n_sim;++k) {
	din+=n[k]*exp(fene[k]-enek[k][i][t]);
      }
      W[j*(n_sim+2)+n_sim]=1.0/normz_const/din;
    }
  }

  *max=od_data[0][0];*min=od_data[0][0];
  for (i=0;i<n_sim;++i) {
    for (j=0;j<n[i];++j) {
      if (*max < od_data[i][j]) *max=od_data[i][j];
      if (*min > od_data[i][j]) *min=od_data[i][j];
    }
  }
  *frame=(int)((*max-*min)/width)+1;
  hist_CA=(double *)gcemalloc(sizeof(double)*(*frame));
  hist_error=(double *)gcemalloc(sizeof(double)*(*frame));
  for (i=0;i<*frame;++i) hist_CA[i]=0.0;
  for (i=0;i<*frame;++i) hist_error[i]=0.0;

  for (i=0;i<n_sim;++i) {
    for (t=0;t<n[i];++t) {
      din=0.0;
      for (k=0;k<n_sim;++k)
	din+=n[k]*exp(fene[k]-enek[k][i][t]);
      hist_CA[((int)((od_data[i][t]-*min)/width))]+=1.0/din;
    }
  }

  num=n_sim+1;
  for (i=0;i<*frame;++i) {
    nn=-1;
    for (j=0;j<n_sim;++j) {
      for (t=0;t<n[j];++t) {
	++nn;
	din=0.0;for (k=0;k<n_sim;++k) din+=n[k]*exp(fene[k]-enek[k][j][t]);
	if (i==((int)((od_data[j][t]-*min)/width)))
	  W[nn*(n_sim+2)+num]=1.0/din/hist_CA[i];
	else
	  W[nn*(n_sim+2)+num]=0.0;
      }
    }

    N=(int *)gcemalloc(sizeof(int)*(n_sim+2));
    for (k=0;k<n_sim;++k) N[k]=n[k];
    N[k]=0;N[k++]=0;
    for (j=0;j<n_total;++j)
    	for (k=0;k<n_sim+2;++k)
	  WT[k*n_total+j]=W[j*(n_sim+2)+k];
    //    covm=MBAR_ACM(n_sim+2,N,W);
    covm=MBAR_ACM2(n_sim+2,N,W,WT);
      
    hist_error[i]=
      (hist_CA[i]/normz_const)*(hist_CA[i]/normz_const)
      *(covm[(n_sim+1)*((n_sim+2))+n_sim+1]
	+covm[(n_sim)*((n_sim+2))+n_sim]
	-2.0*covm[(n_sim+1)*((n_sim+2))+n_sim]);
  }

  normz_const_for_hg=0.0;
  for (i=0;i<*frame;++i) normz_const_for_hg+=hist_CA[i];

  for (i=0;i<*frame;++i) hist_CA[i]=hist_CA[i]/normz_const_for_hg;

  for (i=0;i<*frame;++i) hist_error[i]=hist_error[i]/normz_const_for_hg;

  return hist_error;

}

double *MBAR_setw(double *fene, double ***enek, int n_sim, int *n, double *W, double *WT){
  int i,j,k,l,t;
  int n_total;
  double din;

  n_total=0;for (i=0;i<n_sim;++i)n_total+=n[i];

  l=-1;
  for (i=0;i<n_sim;++i) {
    for (t=0;t<n[i];++t) {
      ++l;
      for (j=0;j<n_sim;++j) {
  	din=0.0;
  	for (k=0;k<n_sim;++k) din+=n[k]*exp(fene[k]-enek[k][i][t]);
  	W[l*n_sim+j]=exp(fene[j]-enek[j][i][t])/din;
  	WT[j*n_total+l] =W[l*n_sim+j];
      }
    }
  }

}

double *MBAR_AVE_twod(double *fene, double ***enek, int n_sim, int *n, double ***td_data, double *width,double *max, double *min, int *frame){
  int i,j,k,t;
  
  double *ef;
  double ***euk;

  double *hist;

  int n_total;
  double din;
  double normz_const=0.0;
  double sum;

  ef=(double *)gcemalloc(sizeof(double)*n_sim);
  euk=(double ***)gcemalloc(sizeof(double **)*n_sim);
  for (i=0;i<n_sim;++i) euk[i]=(double **)gcemalloc(sizeof(double *)*n_sim);
  for (i=0;i<n_sim;++i) for (j=0;j<n_sim;++j) euk[i][j]=(double *)gcemalloc(sizeof(double)*n[/*j*/i]);

  for (i=0;i<n_sim;++i) 
    for (j=0;j<n_sim;++j) 
      for (t=0;t<n[i];++t) 
	euk[i][j][t]=exp(-1.0*enek[i][j][t]);

  for (i=0;i<n_sim;++i) ef[i]=exp(fene[i]);

  n_total=0; for (i=0;i<n_sim;++i)  n_total+=n[i];

  normz_const=0.0;
  for (i=0;i<n_sim;++i) {
    for (t=0;t<n[i];++t) {
      din=0.0;
      for (j=0;j<n_sim;++j)
	din+=n[j]*ef[j]*euk[j][i][t];
      normz_const+=1.0/din;
    }
  }

  for (i=0;i<2;++i) {
    max[i]=td_data[i][0][0];
    min[i]=td_data[i][0][0];
  }
  for (i=0;i<2;++i) {
    for (j=0;j<n_sim;++j) {
      for (k=0;k<n[j];++k) {
	if (max[i] < td_data[i][j][k]) max[i]=td_data[i][j][k];
	if (min[i] > td_data[i][j][k]) min[i]=td_data[i][j][k];
      }
    }
  }

  for (i=0;i<2;++i) frame[i]=(int)((max[i]-min[i])/width[i])+1;
  hist=(double *)gcemalloc(sizeof(double)*(frame[0])*(frame[1]));
  for (i=0;i<frame[0]*frame[1];++i) hist[i]=0.0;

  for (i=0;i<n_sim;++i) {
    for (t=0;t<n[i];++t) {
      din=0.0;
      for (k=0;k<n_sim;++k) din+=n[k]*ef[k]*euk[k][i][t];
      hist[((int)((td_data[0][i][t]-min[0])/width[0]))*(frame[1])+((int)((td_data[1][i][t]-min[1])/width[1]))]+=1.0/din/normz_const/*/width[0]/width[1]*/;
    }
  }

  sum=0.0; for (i=0;i<frame[0]*frame[1];++i) sum+=hist[i];
  
  for (i=0;i<frame[0]*frame[1];++i)  hist[i]/=sum;

  //  for (i=0;i<frame[0]*frame[1];++i) hist[i]=hist[i]/width[0]/width[1];

  return hist;
}

double **MBAR_AVE_twod_multi_temp(double *ef, double ***euk, int n_sim, int *n, double ***td_data, double widthx, double widthy,double *maxx,double *minx,double *maxy,double *miny, int *framex, int *framey, int numk){
  int i,j,k,t,t_total;
  
  double **hist;

  int n_total;
  double din;
  double expbetaene;
  double normz_const=0.0;
  double sum;

  n_total=0; for (i=0;i<n_sim;++i)  n_total+=n[i];

  *maxx=td_data[0][0][0];*minx=td_data[0][0][0];
  *maxy=td_data[0][0][1];*miny=td_data[0][0][1];
  for (i=0;i<n_sim;++i) {
    for (j=0;j<n[i];++j) {
      if (*maxx < td_data[i][j][0]) *maxx=td_data[i][j][0];
      if (*minx > td_data[i][j][0]) *minx=td_data[i][j][0];

      if (*maxy < td_data[i][j][1]) *maxy=td_data[i][j][1];
      if (*miny > td_data[i][j][1]) *miny=td_data[i][j][1];
    }
  }
  *framex=(int)round((*maxx-*minx)/widthx)+1;
  *framey=(int)round((*maxy-*miny)/widthy)+1;
  hist=(double **)gcemalloc(sizeof(double *)*(*framex));
  for (i=0;i<*framex;++i) hist[i]=(double *)gcemalloc(sizeof(double)*(*framey));

  for (i=0;i<n_sim;++i) {
    for (t=0;t<n[i];++t) {
      din=0.0; for (j=0;j<n_sim;++j) din+=n[j]*ef[j]*euk[j][i][t];
      hist[((int)round((td_data[i][t][0]-*minx)/widthx))][((int)round((td_data[i][t][1]-*miny)/widthy))]+=1.0/din*ef[numk];
    }
  }

  sum=0.0; for (i=0;i<*framex;++i) for (j=0;j<*framey;++j) sum+=hist[i][j];

  for (i=0;i<*framex;++i) for (j=0;j<*framey;++j) hist[i][j]/=sum;

  return hist;
}

double MBAR_A(double *expF, double ***expene, int n_sim, int *n, double *A, double *dA){
  int i,j,k,l,m,nn;
  double n1,n2,din,num;

  n1=0.0;
  n2=0.0;
  for (i = 0; i < n_sim; ++i) {
    for (j = 0; j < n_sim; ++j) {
      for (nn = 0; nn < n[j]; ++nn) {
	din=0.0; for (k=0;k<n_sim;++k) din+=n[k]*expene[k][j][nn]*expF[k];
	n2+=log(1.0/din);
      }
    }
    n1+=n[i]*log(expF[i]);
  }
  *A=-n1-n2;

  for (i = 0; i < n_sim; ++i) {
    num=0.0;
    for (j = 0; j < n_sim; ++j) {
      for (nn = 0; nn < n[j]; ++nn) {
	din=0.0;
	for (k = 0; k < n_sim; ++k) {	  
	  din+=n[k]*expene[k][j][nn]*expF[k];
	}
	num+=expene[i][j][nn]/din;
      }
    }
    dA[i]=n[i]*(num-1.0/expF[i]);
  }

  return 0.0;
}

double *MBAR_AVE_twod_wmaxmin(double *fene, double ***enek, int n_sim, int *n, double ***td_data, double *width,double *max, double *min, int *frame){
  int i,j,k,t;
  
  double *ef;
  double ***euk;

  double *hist;

  int n_total;
  double din;
  double normz_const=0.0;
  double sum;
  double hist_min,hist_max;
  int index_min,index_max;
  int c;

  ef=(double *)gcemalloc(sizeof(double)*n_sim);
  euk=(double ***)gcemalloc(sizeof(double **)*n_sim);
  for (i=0;i<n_sim;++i) euk[i]=(double **)gcemalloc(sizeof(double *)*n_sim);
  for (i=0;i<n_sim;++i) for (j=0;j<n_sim;++j) euk[i][j]=(double *)gcemalloc(sizeof(double)*n[/*j*/i]);

  for (i=0;i<n_sim;++i) {
    for (j=0;j<n_sim;++j) { 
      for (t=0;t<n[i];++t) {
	if (enek[i][j][t] < 300)
	  euk[i][j][t]=exp(-1.0*enek[i][j][t]);
	else // debug 2012 05 27
	  euk[i][j][t]=exp(-1.0*0.0);
      }
    }
  }

  for (i=0;i<n_sim;++i) ef[i]=exp(fene[i]);

  n_total=0; for (i=0;i<n_sim;++i)  n_total+=n[i];

  normz_const=0.0;
  for (i=0;i<n_sim;++i) {
    for (t=0;t<n[i];++t) {
      din=0.0;
      for (j=0;j<n_sim;++j)
	din+=n[j]*ef[j]*euk[j][i][t];
      normz_const+=1.0/din;
    }
  }

  /************************************************************************/
  /* for (i=0;i<2;++i) {						  */
  /*   max[i]=td_data[i][0][0];						  */
  /*   min[i]=td_data[i][0][0];						  */
  /* }									  */
  /* for (i=0;i<2;++i) {						  */
  /*   for (j=0;j<n_sim;++j) {						  */
  /*     for (k=0;k<n[j];++k) {						  */
  /* 	if (max[i] < td_data[i][j][k]) max[i]=td_data[i][j][k];		  */
  /* 	if (min[i] > td_data[i][j][k]) min[i]=td_data[i][j][k];		  */
  /*     }								  */
  /*   }								  */
  /* }									  */
  /************************************************************************/

  for (i=0;i<2;++i) frame[i]=(int)((max[i]-min[i])/width[i])+1;
  hist=(double *)gcemalloc(sizeof(double)*(frame[0])*(frame[1]));
  for (i=0;i<frame[0]*frame[1];++i) hist[i]=0.0;

  for (i=0;i<n_sim;++i) {
    for (t=0;t<n[i];++t) {
      din=0.0;
      for (k=0;k<n_sim;++k) din+=n[k]*ef[k]*euk[k][i][t];
      hist[((int)((td_data[0][i][t]-min[0])/width[0]))*(frame[1])+((int)((td_data[1][i][t]-min[1])/width[1]))]+=1.0/din/normz_const/*/width[0]/width[1]*/;

      /***************************************************************************************************************************/
      /* if ((c=((int)((td_data[0][i][t]-min[0])/width[0]))*(frame[1])+((int)((td_data[1][i][t]-min[1])/width[1])))==171)	 */
      /* 	if (hist[171]>0.1)												 */
      /* 	  printf("yes\n");												 */
      /***************************************************************************************************************************/

    }
  }

  sum=0.0; for (i=0;i<frame[0]*frame[1];++i) sum+=hist[i];

  /***********************************************************************************************************/
  /* hist_min=0.0;											     */
  /* hist_max=0.0;											     */
  /* for (i=0;i</\*=*\/frame[0];++i) {									     */
  /*   for (j=0;j</\*=*\/frame[1];++j) {								     */
  /*     if ( (hist_min > hist[i*frame[1]+j] && hist[i*frame[1]+j]!=0) || hist_min==0.0 ) {		     */
  /* 	hist_min=hist[i*frame[1]+j];									     */
  /* 	index_min=i*frame[1]+j;										     */
  /*     }												     */
  /*     if ( hist_max < hist[i*frame[1]+j] && hist[i*frame[1]+j]!=0 && hist[i*frame[1]+j] < 1000 ) {	     */
  /* 	hist_max=hist[i*frame[1]+j];									     */
  /* 	index_max=i*frame[1]+j;										     */
  /*     }												     */
  /*   }												     */
  /* }													     */
  /***********************************************************************************************************/
  
  for (i=0;i<frame[0]*frame[1];++i)  hist[i]/=sum;

  for (i=0;i<frame[0]*frame[1];++i) hist[i]=hist[i]/width[0]/width[1];

  hist_min=0.0;
  hist_max=0.0;
  for (i=0;i</*=*/frame[0];++i) {
    for (j=0;j</*=*/frame[1];++j) {
      if ( (hist_min > hist[i*frame[1]+j] && hist[i*frame[1]+j]!=0) || hist_min==0.0 ) {
  	hist_min=hist[i*frame[1]+j];
  	index_min=i*frame[1]+j;
      }
      if ( hist_max < hist[i*frame[1]+j] && hist[i*frame[1]+j]!=0 && hist[i*frame[1]+j] < 1000 ) {
  	hist_max=hist[i*frame[1]+j];
  	index_max=i*frame[1]+j;
      }
    }
  }

  return hist;
}
