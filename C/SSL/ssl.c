
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SSL.h"
#include "LA.h"
#include "efunc.h"

double epsilon=/*0.0001*/0.00001;
double epsilon_subg=/*0.0001*/0.00001;

void pickupS(double *COVM,int numv, int numcolum,double *S);
void pickupW(double *Sigma,int numv, int numcolum,double *W);

void setLambda(double lambda,double *l,int numv,int numcolum,double *Lambda);
void setSigma(double *omega,int numv, int numcolum,double *Sigma);
int jugement_conversion(double *Sigma, double *Sigmap,double *COVM,int numv);

/*****************************************/
/* input	                         */
/* S 1x(M-1),  W (M-1)*(M-1), rou 1x1	 */
/*                    ||                 */
/*                    \/                 */
/* output                                */
/* beta 1x(M-1)                          */
/*****************************************/
int ssl_subgraalg(double *S, double rou, double *W, int numv, double *beta) {
  int i,j,k,numite,numcolum,judge;
  double A,*betaold;
  FILE *ctest;

  ctest=efopen("subgraalg.txt","a");  
  for (i=0;i<numcolum;++i)
    beta[i]=0.0;

  judge =1;
  numite=0;
  numcolum=numv-1;
  betaold=ecalloc(numcolum,sizeof(double));

  while (judge != 0 || numite < 20) {
    for (i=0;i<numcolum;++i) {
      A = S[i];
      for (j=0;j<numcolum;++j) {
	if (j!=i)
	  A-=W[i*numcolum+j]*beta[j];
      }
      if (A > rou)
	beta[i]=(A-rou)/W[i*numcolum+i];
      else if (A >= -rou && A <= rou)
	beta[i]=0.0;
      else
	beta[i]=(A+rou)/W[i*numcolum+i];
    }
    ++numite;
    fprintf(ctest,"%d : ",numite);
    for (i=0;i<numcolum;++i)
      fprintf(ctest,"%e : ",beta[i]);
    fprintf(ctest,"\n");
    judge = 0;
    for (i=0;i<numcolum;++i) {
      if (fabs(beta[i]-betaold[i])>epsilon_subg)
    	judge = 1;
    }
    for (i=0;i<numcolum;++i)
      betaold[i]=beta[i];
  }

  free(betaold);
  fclose(ctest);
  return 0;
}


/***************************************************/
/* input	                                   */
/* S 1x(M-1), W (M-1)*(M-1),  rou 1x1,  sigma 1x1  */
/*                        ||                       */
/*                        \/                       */
/* output                                          */
/* lambda 1x1, l 1x(M-1), omega 1x(M-1)            */
/***************************************************/
double ssl_gralasso(double *S, double rou, int numv, double *W, double sigma, double lambda, double *l, double *omega) {
  int i;
  double *beta,*Wl;

  beta = (double *)ecalloc(numv-1,sizeof(double));
  Wl = (double *)ecalloc(numv-1,sizeof(double));
  ssl_subgraalg(S, rou, W, numv, beta);
  
  lambda = 1.0/(sigma-vtmvmult(beta,W,beta,numv-1));
  for (i=0;i<numv-1;++i)
    l[i]=-lambda*beta[i];
  mvmult(W,l,Wl,numv-1);
  for (i=0;i<numv-1;++i)
    omega[i]=-Wl[i]/lambda;

  free(beta);
  free(Wl);
  return lambda;
}

/*******************/
/* input	   */
/* LS MxM, rou 1x1 */
/*      ||         */
/*      \/         */
/* output          */
/* Lambda MxM      */
/*******************/
int ssl_gralassomain(double *COVM, double rou, int numv, double *Lambda, double *Sigma ) {
  int i,j,k;
  int numite;
  double sigma,lambda;
  double *W,*S,*l,*omega,*Sigmap,*unit;

  Sigmap=(double *)ecalloc(numv*numv,sizeof(double));
  unit=(double *)ecalloc(numv*numv,sizeof(double));
  W=(double *)ecalloc((numv-1)*(numv-1),sizeof(double));
  S=(double *)ecalloc((numv-1),sizeof(double));
  l=(double *)ecalloc((numv-1),sizeof(double));
  omega=(double *)ecalloc((numv-1),sizeof(double));


  for (i=0;i<numv;++i)
    for (j=i;j<numv;++j)
      if (i==j)
	Sigma[i*numv+i]=COVM[i*numv+i]+rou;
      else
	Sigma[i*numv+j]=COVM[i*numv+j];

  for (j=0;j<numv;++j)
    for (k=0;k<numv;++k)
      Sigmap[j*numv+k]=Sigma[j*numv+k];

  numite=0;
  do {
    for (j=0;j<numv;++j)
      for (k=0;k<numv;++k)
	Sigmap[j*numv+k]=Sigma[j*numv+k];
    for (i=0;i<numv;++i) {
      pickupS(COVM,numv,i,S);
      pickupW(Sigma,numv,i,W);
      lambda=ssl_gralasso(S, rou, numv, W, Sigma[i*numv+i], lambda, l, omega);
      setLambda(lambda,l,numv,i,Lambda);
      setSigma(omega,numv,i,Sigma);
    }
    ++numite;
  }while (jugement_conversion(Sigma,Sigmap,COVM,numv)!=0 || numite < 20);

  for (i=0;i<numv;++i) {
    for (j=i;j<numv;++j) {
      Lambda[j*numv+i]=Lambda[i*numv+j];
      Sigma[j*numv+i]=Sigma[i*numv+j];
    }
  }

  mmult(Lambda,Sigma,unit,numv);

  free(Sigmap);
  free(S);
  free(l);
  free(omega);
  free(W);

  return numite;
}

void pickupS(double *COVM,int numv, int numcolum,double *S) {
  int i;

  for (i=0;i<numcolum;++i)
    S[i]=COVM[i*numv+numcolum];
  for (i=numcolum;i<numv-1;++i)
    S[i]=COVM[numcolum*numv+i+1];
}

void pickupW(double *Sigma,int numv, int numcolum,double *W) {
  int i,j;

  for (i=0;i<numcolum;++i)
    for (j=i;j<numcolum;++j)
      W[i*(numv-1)+j]=Sigma[i*numv+j];
  for (i=0;i<numcolum;++i)
    for (j=numcolum;j<numv-1;++j)
      W[i*(numv-1)+j]=Sigma[i*numv+j+1];
  for (i=numcolum;i<numv-1;++i)
    for (j=i;j<numv-1;++j)
      W[i*(numv-1)+j]=Sigma[(i+1)*numv+(j+1)];
}

void setLambda(double lambda,double *l,int numv,int numcolum,double *Lambda){
  int i,j;

  for (i=0;i<numcolum;++i)
    Lambda[i*numv+numcolum]=l[i];
  for (i=numcolum;i<numv-1;++i)
    Lambda[numcolum*numv+i+1]=l[i];

  Lambda[numcolum*numv+numcolum]=lambda;
}

void setSigma(double *omega,int numv, int numcolum,double *Sigma) {
  int i,j;

  for (i=0;i<numcolum;++i)
    Sigma[i*numv+numcolum]=omega[i];
  for (i=numcolum;i<numv-1;++i)
    Sigma[numcolum*numv+i+1]=omega[i];
}

int jugement_conversion(double *Sigma, double *Sigmap,double *COVM,int numv){
  int i,j,numcolum,num=0;
  double det=0.0,aveSdiag=0.0;
  FILE *ctest;

  ctest=efopen("glasso.txt","w");

  numcolum=numv-1;
  num=0;
  for (i=0;i<numv;++i) {
    for (j=i+1;j<numv;++j) {
      ++num;
      det = (num*det+fabs(Sigma[i*numv+j]-Sigmap[i*numv+j]))/(num+1);
    }
  }
  num=0;
  for (i=0;i<numv;++i) {
    for (j=i+1;j<numv;++j) {
      ++num;
      aveSdiag=(num*aveSdiag+fabs(COVM[i*numv+j]))/(num+1);
    }
  }

  fprintf(ctest,"%e \n",det);
  fclose(ctest);

  if (det < epsilon*aveSdiag)
    return 0;
  else
    return 1;
}

double ssl_covm(double *timeseries, int nums, int numv, double *COVM){
  int i,j,k;

  for (j=0;j<numv;++j)
    for (k=j;k<numv;++k)
      COVM[j*numv+k]=0.0;
  for (i=0;i<nums;++i)
    for (j=0;j<numv;++j)
      for (k=j;k<numv;++k)
	COVM[j*numv+k]+=timeseries[i*numv+j]*timeseries[i*numv+k];
  for (j=0;j<numv;++j)
    for (k=j;k<numv;++k)
      COVM[j*numv+k]=COVM[j*numv+k]/nums;

  return 0.0;
}

double ssl_ave(double *timeseries, int nums, int numv, double *AVE){
  int i,j;

  for (i=0;i<nums;++i)
    for (j=0;j<numv;++j)
      AVE[j]=(i*AVE[j]+timeseries[i*numv+j])/(i+1);

  return 0.0;
}


double ssl_normalize(double *timeseries, int nums, int numv, double *timeseriesnorm) {
  int i,j;
  double *ave,*var;
  FILE *output;

  ave=(double *)ecalloc(numv,sizeof(double));
  var=(double *)ecalloc(numv,sizeof(double));

  ssl_ave(timeseries,nums,numv,ave);
  for (i=0;i<nums;++i)
    for (j=0;j<numv;++j)
      timeseriesnorm[i*numv+j]=timeseries[i*numv+j]-ave[j];
  for (i=0;i<numv;++i)
    var[i]=0.0;
  for (i=0;i<nums;++i) {
    for (j=0;j<numv;++j) {
      var[j]+=timeseriesnorm[i*numv+j]*timeseriesnorm[i*numv+j];
    }
  }
  for (j=0;j<numv;++j)
    var[j]=sqrt(var[j]/nums);
  for (i=0;i<nums;++i)
    for (j=0;j<numv;++j)
      timeseriesnorm[i*numv+j]=timeseriesnorm[i*numv+j]/var[j];

  for (i=0;i<numv;++i)
    var[i]=0.0;
  for (i=0;i<nums;++i) {
    for (j=0;j<numv;++j) {
      var[j]+=timeseriesnorm[i*numv+j]*timeseriesnorm[i*numv+j];
    }
  }
  for (j=0;j<numv;++j)
    var[j]=sqrt(var[j]/nums);

  output=efopen("datanorm.txt","w");
  for (i=0;i<nums;++i) {
    fprintf(output,"%d ",i+1);
    for (j=0;j<numv;++j) {
      fprintf(output,"%e ",timeseriesnorm[i*numv+j]);
    }
    fprintf(output,"\n");
  }
  fclose(output);

}

 double ssl_KLdiv(double *Lambda, double *Sigma, double *Lambdap, double *Sigmap, int numv,  double *KLdiv, int topten[10], double vtopten[10]) {
  int i,j,k;
  double *WA,*WB,*omegaA,*omegaB,*lA,*lB,*lAB,*lBA,lambdaA,lambdaB,sigmaA,sigmaB,dAB,dBA,sumKLdiv;

  WA=(double *)ecalloc((numv-1)*(numv-1),sizeof(double));
  WB=(double *)ecalloc((numv-1)*(numv-1),sizeof(double));
  omegaA=(double *)ecalloc((numv-1),sizeof(double));
  omegaB=(double *)ecalloc((numv-1),sizeof(double));
  lA=(double *)ecalloc((numv-1),sizeof(double));
  lB=(double *)ecalloc((numv-1),sizeof(double));
  lBA=(double *)ecalloc((numv-1),sizeof(double));
  lAB=(double *)ecalloc((numv-1),sizeof(double));

  for (i=0;i<numv;++i) {
    pickupS(Sigma,numv,i,omegaA);
    pickupW(Sigma,numv,i,WA);
    pickupS(Sigmap,numv,i,omegaB);
    pickupW(Sigmap,numv,i,WB);
    pickupS(Lambda,numv,i,lA);
    pickupS(Lambdap,numv,i,lB);
    lambdaA=Lambda[i*numv+i];
    lambdaB=Lambdap[i*numv+i];
    sigmaA=Sigma[i*numv+i];
    sigmaB=Sigmap[i*numv+i];

    for (j=0;j<numv-1;++j) {
      lBA[j] = lB[j]-lA[j];
      lAB[j] = lA[j]-lB[j];
    }

    dAB = inprod(omegaA,lBA,numv-1)+0.5*(vtmvmult(lB,WA,lB,numv-1)/lambdaB-vtmvmult(lA,WA,lA,numv-1)/lambdaA)+0.5*(log(lambdaA/lambdaB)+sigmaA*(lambdaB-lambdaA));
    dBA = inprod(omegaB,lAB,numv-1)+0.5*(vtmvmult(lA,WB,lA,numv-1)/lambdaA-vtmvmult(lB,WB,lB,numv-1)/lambdaB)+0.5*(log(lambdaB/lambdaA)+sigmaB*(lambdaA-lambdaB));
    
    if (dAB >= dBA)
      KLdiv[i] = dAB;
    else
      KLdiv[i] = dBA;
    //  KLdiv[i] = dBA;
  }

  for (i=0;i<10;++i) {
    vtopten[i] = 0.0;
    topten[i]=0;
  }

  sumKLdiv=0.0;
  for (i=0;i<numv;++i) {
    sumKLdiv += KLdiv[i];
    for (j=0;j<9;++j) {
      if (vtopten[j]>=KLdiv[i]) {
	;
      }
      else {
	for (k=8;k>=j;--k) {
	  vtopten[k+1]=vtopten[k];
	  topten[k+1]=topten[k];
	}
	vtopten[j]=KLdiv[i];
	topten[j]=i;
	break;
      }
    }
  }
  return sumKLdiv;
}
