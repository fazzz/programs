#include <stdio.h>
#include <math.h>

#include "f2c.h"
#include "clapack.h"

#include "ABA_multimer.h"
#include "gener.h"
#include "force.h"
#include "MD.h"
#include "UmsSan.h"

void Calc_acc_debug(double *eig_dihed);

void mmmult(double *m1, double *m2, double *m1m2, int n);
void mvmult(double *m, double *v, double *mv, int n);
void mtrans(double *m1,double *m2, int n);
void msetIni(double *m, int n);
void msetzero(double *m, int n);
int invm(double mat[7][7],double inv[7][7]);
int invm6(double mat[6][6],double inv[6][6]);
int invm3(double mat[3][3],double inv[3][3]);

// 次のステップの 2 面角を計算を行う関数
void calc_dd_theta_cycle(int pflag,int cflag,
			 double *ele, double *ALJ, double *BLJ,
			 double *p_e, double *p_1_4_e, 
			 double *p_LJ,double *p_1_4_LJ,
			 double *f_e, double *f_1_4_e, 
			 double *f_LJ,double *f_1_4_LJ,
			 int numnb, int *indexnb,
			 int num14, int *index14,double *eig,double *eig_14,double *eig_dihed) {
  int i;
  int nNumClut;
  int nNumClutOrigBranch = 0;
  double add[6];

  // 現在の構造でのタンパク質のポテンシャルエネルギー、
  // タンパク質に及ぼす力を計算
  //	calc_force();
  if (cflag == ON) {
    //  printf("yes-d1\n");
    calc_force2(pflag,
		ele, ALJ, BLJ,
		p_e, p_1_4_e, 
		p_LJ,p_1_4_LJ,
		f_e,f_1_4_e, 
		f_LJ,f_1_4_LJ,
		numnb, indexnb,
		num14, index14,eig,eig_14,eig_dihed);
    //  printf("yes-d2\n");
  }
  //	Calc_dihed_Umbrella_Potential();
  
  // 現在の構造でのタンパク質の ABI 、Bias Force を計算 ( 末端 -> 始点 )
  for (nNumClut=prot.DOF/*-1*/;
       nNumClut>=0;
       --nNumClut) {
    // 現在の構造でのタンパク質の ABI を計算
    calc_ABA_cycle(nNumClut);
    // 現在の構造でのタンパク質の bias force を計算
    //		calc_bias_force_cycle(nnumclut);
  }
  // 現在の構造でのタンパク質の abi 、bias force を計算 ( 末端 -> 始点 )
  for (nNumClut=prot.DOF/*-1*/;
       nNumClut>=0;
       --nNumClut){
    // 現在の構造でのタンパク質の abi を計算
    //		calc_aba_cycle(nnumclut);
    // 現在の構造でのタンパク質の bias force を計算
    calc_Bias_force_cycle(nNumClut);
  }

  // 剛体 0 番目の spatial acceleration の設定
  if (TermMoveMode==OFF) {
    for(i=0; i<6; ++i){
      clust[0].sp_acc[i] = 0.0;
    }
  }
  else {
  }
  
  if (TermMoveMode==OFF && ABI[0].hingmat!=6 ) {
    nNumClut= 1;
  }
  else {
    nNumClut = 0;
  }
  // 現在の構造でのタンパク質の spatial acceleration を計算 ( 始点 -> 末端 )
  for(/*nNumClut = 1*//*0*/; nNumClut < prot.DOF; ++nNumClut) {
    if (clust[nNumClut].num_branch > 1) {
      nNumClutOrigBranch = nNumClut;
    }
    
    // 予測子加速度の計算を行う
    //		set_predict_acc(nnumclut, nnumclutorigbranch);
    // spatial acceleration の計算を行う
    calc_sp_acc_cycle(nNumClut);
  }
  if (MODE == NVT) {
    if (TermMoveMode==OFF && ABI[0].hingmat!=6) {
      nNumClut = 1;
    }
    else {
      nNumClut = 0;
    }
    for(/* nNumClut = 1/\*0*\/ */; nNumClut < prot.DOF; ++nNumClut) {
      if (GearMODE == ON) {
	s_NVT = s_NVT_predict[0];
	dot_s_NVT = s_NVT_predict[1]/deltat;
      }
      if (TermMoveMode==ON && nNumClut==0 && TermMoveMode2==12) {
	add[0]=0.0;
	add[1]=0.0;
	add[2]=0.0;
	add[3]=(-vel_Term[2]*vel_Term[3+1]+vel_Term[1]*vel_Term[3+2]);
	add[4]=( vel_Term[2]*vel_Term[3+0]-vel_Term[0]*vel_Term[3+2]);
	add[5]=(-vel_Term[1]*vel_Term[3+0]+vel_Term[0]*vel_Term[3+1]);

	for (i=0;i<6;++i) {
	  acc_Term2[i]-=dot_s_NVT/s_NVT*vel_Term[i];
	  acc_Term[i]=acc_Term2[i]-add[i];
	}
      }
      clust[nNumClut].dddihedang[0] -= dot_s_NVT/s_NVT/*xi_dummy**//*xi_dummy**/*clust[nNumClut].ddihedang[0];
    }
  }

  if (TermMoveMode==ON) {
    /*********************/
    /* Calc_acc_debug(); */
    /*********************/
  }

}

void Calc_acc_debug(double *eig_dihed){
  int i,j;
  /******************************************/
  /* double *Mat,*InvMat;		    */
  /* double *MatPsi,*MatInt;		    */
  /* double *force,*biasf,*T,*theta;	    */
  /* double MatHingPsiIintPsiTHingT2[7][7]; */
  /******************************************/
  double *In0,*In1,*Psi,*PsiIn1,*In1PsiT,*PsiT,*PsiIn1PsiT,*Inv2,*T,*theta,*PsiIn1a,*In1a,*a,*b1,*Psib1,*Ha,*In1Ha;
  double *T2,*theta2;
  double *T3,*theta3;
  double *domega,*S0,*dS,*RS0,*R;
  double *omega,*omegaRS0,*omegaomegaRS0;
  double *sx,*omegasx,*omegasxv,*v;
  double dotS[3];
  double Mat[7][7],Inv[7][7],Mat2[6][6],InvIn1[6][6],Mat3[3][3],InvIn3[3][3],m;
  double q[4],dotq[4],temp[4][4],temp2[4];

  In0=malloc(sizeof(double)*6*6);
  In1=malloc(sizeof(double)*6*6);
  Psi=malloc(sizeof(double)*6*6);
  PsiT=malloc(sizeof(double)*6*6);
  PsiIn1=malloc(sizeof(double)*6*6);
  In1PsiT=malloc(sizeof(double)*6*6);
  PsiIn1PsiT=malloc(sizeof(double)*6*6);
  Inv2=malloc(sizeof(double)*7*7);
  T=malloc(sizeof(double)*6);
  theta=malloc(sizeof(double)*6);
  a=malloc(sizeof(double)*6);
  In1a=malloc(sizeof(double)*6);
  PsiIn1a=malloc(sizeof(double)*6);
  b1=malloc(sizeof(double)*6);
  Psib1=malloc(sizeof(double)*6);
  Ha=malloc(sizeof(double)*6);
  In1Ha=malloc(sizeof(double)*6);
  T2=malloc(sizeof(double)*6);
  theta2=malloc(sizeof(double)*6);
  T3=malloc(sizeof(double)*3);
  theta3=malloc(sizeof(double)*6);
  domega=malloc(sizeof(double)*3*3);
  S0=malloc(sizeof(double)*3);
  RS0=malloc(sizeof(double)*3);
  R=malloc(sizeof(double)*3*3);
  dS=malloc(sizeof(double)*3);
  omega=malloc(sizeof(double)*3*3);
  omegaRS0=malloc(sizeof(double)*3);
  omegaomegaRS0=malloc(sizeof(double)*3);
  sx=malloc(sizeof(double)*3*3);
  omegasx=malloc(sizeof(double)*3*3);
  omegasxv=malloc(sizeof(double)*3);
  v=malloc(sizeof(double)*3);

  for(i=0;i<6;++i) {
    for(j=0;j<6;++j) {
      In0[i*6+j]=clust[0].InertiaMatrix[i][j];
      In1[i*6+j]=clust[1].InertiaMatrix[i][j];
      Psi[i*6+j]=clust[1].TransMatrix[0][i][j];
    }
  }
  
  mmmult(Psi,In1,PsiIn1,6);
  mtrans(Psi,PsiT,6);
  mmmult(In1,PsiT,In1PsiT,6);
  mmmult(PsiIn1,PsiT,PsiIn1PsiT,6);

  for(i=0;i<6;++i) {
    Mat[6][i]=PsiIn1[i*6+2];
    Mat[i][6]=In1PsiT[2*6+i];
    for(j=0;j<6;++j) {
      Mat[i][j]=In0[i*6+j]+PsiIn1PsiT[i*6+j];
    }
  }
  Mat[6][6]=In1[2*6+2];
  invm(Mat,Inv);
  for(i=0;i<7;++i) {
    for(j=0;j<7;++j) {
      Inv2[i*7+j]=Inv[i][j];
    }
  }
  for (i=0;i<6;++i) {
    a[i]=clust[1].Coriolis_acc[i];
  }
  for (i=0;i<3;++i) {
    Ha[i]=clust[1].Coriolis_acc[i];
  }
  for (i=3;i<6;++i) {
    Ha[i]=0.0;
  }
  for (i=0;i<6/*3*/;++i) {
    b1[i]=clust[1].Coriolis_b[i];
  }
  mvmult(PsiIn1,a,PsiIn1a,6);
  mvmult(Psi,b1,Psib1,6);
  for(i=0;i<6;++i) {
    T[i]=-PsiIn1a[i]-Psib1[i]-clust[0].Coriolis_b[i];
  }
  mvmult(In1,Ha,In1Ha,6);
  mvmult(In1,a,In1a,6);
  if (dhstopflag!=ON) {
    Calc_dihed_Potential_Force(eig_dihed);
  }
  T[6]=-clust[1].f_c.f_dihed-In1a[2]-clust[1].Coriolis_b[2];
  mvmult(Inv2,T,theta,7);

  for (i=0;i<3/*6*/;++i)
    T2[i]=-clust[0].Coriolis_b[i];

  domega[0]= 0.0;  domega[1]=-theta3[2];  domega[2]= theta3[1];
  domega[3]= theta3[2];  domega[4]= 0.0;  domega[5]=-theta3[0];
  domega[6]=-theta3[1];  domega[7]= theta3[0];  domega[8]= 0.0;
  omega[0]= 0.0;  omega[1]=-vel_Term[2];  omega[2]= vel_Term[1];
  omega[3]= vel_Term[2];  omega[4]= 0.0;  omega[5]=-vel_Term[0];
  omega[6]=-vel_Term[1];  omega[7]= vel_Term[0];  omega[8]= 0.0;
  sx[0]= 0.0;  sx[1]=-clust[0].qCOM[2];  sx[2]= clust[0].qCOM[1];
  sx[3]= clust[0].qCOM[2];  sx[4]= 0.0;  sx[5]=-clust[0].qCOM[0];
  sx[6]=-clust[0].qCOM[1];  sx[7]= clust[0].qCOM[0];  sx[8]= 0.0;
  for (i=0;i<3;++i)
    v[i]=vel_Term[i+3];
  m=0.0;
  for (i=0;i<8;++i)
    m+=clust[0].mass_clust[i];

  mmmult(omega,sx,omegasx,3);
  mvmult(omegasx,v,omegasxv,3);

  for (i=0;i<3;++i)
    T2[i]-=m*omegasxv[i];

  /***********************************/
  /* for (i=0;i<6;++i)		     */
  /*   T2[i]=clust[0].Coriolis_b[i]; */
  /***********************************/

  for (i=0;i<6;++i)
    for (j=0;j<6;++j)
      Mat2[i][j]=clust[0].InertiaMatrix[i][j];
  /***********************/
  /* for (i=0;i<3;++i)	 */
  /*   for (j=3;j<6;++j) */
  /*     Mat2[i][j]=0.0; */
  /***********************/
  invm6(Mat2,InvIn1);
  mvmult(InvIn1,T2,theta2,6);

  for (i=0;i<3;++i)
    T3[i]=-clust[0].Coriolis_b[i];
  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      Mat3[i][j]=clust[0].InertiaMatrix[i][j];
  invm3(Mat3,InvIn3);
  mvmult(InvIn3,T3,theta3,3);
  for (i=0;i<3;++i)
    S0[i]=clust[0].qCOM[i];

  if (nNumStep==1) {
    for (i=0;i<3;++i)
      q[i]=0.0;
    q[3]=1.0;
  }

  temp[0][0]=-q[2];
  temp[0][1]=-q[3];
  temp[0][2]= q[1];
  temp[0][3]= q[0];
  temp[1][0]= q[3];
  temp[1][1]=-q[2];
  temp[1][2]=-q[0];
  temp[1][3]= q[1];
  temp[2][0]= q[0];
  temp[2][1]= q[1];
  temp[2][2]= q[3];
  temp[2][3]= q[2];
  temp[3][0]=-q[1];
  temp[3][1]= q[0];
  temp[3][2]=-q[2];
  temp[3][3]= q[3];

  for (i=0;i<3;++i)
    temp2[i]=vel_Term[i];
  temp2[3]=0.0;
  for (i=0;i<3;++i)
    dotq[i]=0.0;
  for (i=0;i<4;++i)
    for (j=0;j<4;++j)
      dotq[i]+=0.5*temp[i][j]*temp2[j];

  for(i=0;i<4;++i)
    q[i]+=dotq[i]*deltat;

  R[0] = q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3];
  R[1*3+1] = q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3];
  R[2*3+2] = q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];
  
  R[1] = 2.0*(q[1]*q[2]-q[0]*q[3]);
  R[1*3+0] = 2.0*(q[2]*q[1]+q[0]*q[3]);
  
  R[2] = 2.0*(q[1]*q[3]+q[0]*q[2]);
  R[2*3+0] = 2.0*(q[3]*q[1]-q[0]*q[2]);
  
  R[1*3+2] = 2.0*(q[2]*q[3]-q[0]*q[1]);
  R[2*3+1] = 2.0*(q[3]*q[2]+q[0]*q[1]);

  mvmult(domega,S0,dS,3);
  mvmult(R,S0,RS0,3);
  mvmult(omega,RS0,omegaRS0,3);
  mvmult(omega,omegaRS0,omegaomegaRS0,3);
  for (i=0;i<3;++i)
    theta3[i+3]=-dS[i]-omegaomegaRS0[i];
  for (i=0;i<3;++i)
    dotS[i]=clust[0].sp_velo[i+3]+omegaRS0[i];
   
  if (TermMoveMode2 == 1) {
    for (i=0;i<6;++i) {
      acc_Term[i]=theta[i];
    }
    clust[1].dddihedang[0]=theta[6];
  }
  else if (TermMoveMode2 == 2) {
    for (i=0;i<6;++i) {
      acc_Term[i]=theta2[i];
    }
  }
  else if (TermMoveMode2 == 3) {
    for (i=0;i<6;++i) {
      acc_Term[i]=theta3[i];
    }

    /***********************/
    /* for (i=3;i<6;++i) { */
    /*   acc_Term[i]=0.0;  */
    /* }		   */
    /***********************/

  }


  free(In0);
  free(In1);
  free(Psi);
  free(PsiT);
  free(PsiIn1);
  free(PsiIn1PsiT);
  free(a);
  free(In1a);
  free(PsiIn1a);
  free(b1);
  free(Psib1);
  free(Ha);
  free(In1Ha);
  free(T3);
  free(theta3);
  free(domega);
  free(S0);
  free(dS);
  free(omega);
  free(RS0);
  free(omegaRS0);
  free(omegaomegaRS0);
  free(sx);
  free(omegasx);
  free(omegasxv);
  free(v);

  /********************************************************************************/
  /* Mat=malloc(sizeof(double)*12*12);						  */
  /* InvMat=malloc(sizeof(double)*7*7);						  */
  /* MatPsi=malloc(sizeof(double)*12*12);					  */
  /* MatInt=malloc(sizeof(double)*12*12);					  */
  /* force=malloc(sizeof(double)*12);						  */
  /* biasf=malloc(sizeof(double)*12);						  */
  /* T=malloc(sizeof(double)*12);						  */
  /* theta=malloc(sizeof(double)*7);						  */
  /* 										  */
  /* msetIni(MatPsi,12);							  */
  /* msetzero(MatInt,12);							  */
  /* 										  */
  /* for (i=0;i<6;++i) {							  */
  /*   for (j=0;j<6;++j) {							  */
  /*     MatPsi[(i+6)*12+j]=clust[1].TransMatrix[0][i][j];			  */
  /*     MatInt[i*12+j]=clust[0].InertiaMatrix[i][j];				  */
  /*     MatInt[(i+6)*12+(j+6)]=clust[1].InertiaMatrix[i][j];			  */
  /*   }									  */
  /* }										  */
  /* 										  */
  /* force[8]=-clust[1].f_c.f_dihed;						  */
  /* for (i=0;i<3;++i) {							  */
  /*   biasf[i]=clust[0].Coriolis_b[i]-clust[0].f_c.sp_f_clust[0].N_clust[i];	  */
  /*   biasf[3+i]=clust[0].Coriolis_b[3+i]-clust[0].f_c.sp_f_clust[0].f_clust[i]; */
  /*   biasf[6+i]=clust[1].Coriolis_b[i]-clust[1].f_c.sp_f_clust[0].N_clust[i];	  */
  /*   biasf[9+i]=clust[1].Coriolis_b[3+i]-clust[1].f_c.sp_f_clust[0].f_clust[i]; */
  /* }										  */
  /* /\*****************************\/						  */
  /* /\* mvmult(MatHing,force,12); *\/						  */
  /* mvmult(MatPsi,biasf,12);							  */
  /* /\* mvmult(MatHing,biasf,12); *\/						  */
  /* for (i=0;i<12;++i) {							  */
  /*   T[i]=-biasf[i];								  */
  /* }										  */
  /* T[8]+=clust[1].f_c.f_dihed;						  */
  /* /\*****************************\/						  */
  /* 										  */
  /* for (i=0;i<12;++i)								  */
  /*   for (j=0;j<12;++j)							  */
  /*     Mat[i*12+j]=MatInt[i*12+j];						  */
  /* mmmult(MatPsi,Mat,12);							  */
  /* mtrans(MatPsi,12);								  */
  /* mmmult(Mat,MatPsi,12);							  */
  /* //  for (i=0;i<12;++i)							  */
  /* //    for (j=0;j<12;++j)							  */
  /* //      test2[i][j]=MatHing[i*12+j];					  */
  /* //  invm(test2);								  */
  /* for (i=0;i<12;++i)								  */
  /*   for (j=0;j<12;++j)							  */
  /*     Mat[i*12+j]=MatPsi[i*12+j];						  */
  /* //  mmmult(MatHing,Mat,12);						  */
  /* //  mmmult(Mat,MatHing,12);						  */
  /* //  for (i=0;i<12;++i)							  */
  /* //    for (j=0;j<12;++j)							  */
  /* //    Mat[i*12+j]=MatHing[i*12+j];						  */
  /* for (i=0;i<6;++i)								  */
  /*   for (j=0;j<6;++j)							  */
  /*     MatHingPsiIintPsiTHingT2[i][j]=Mat[i*12+j];				  */
  /* for (j=0;j<7;++j) {							  */
  /*   MatHingPsiIintPsiTHingT2[6][i]=Mat[6*12+j];				  */
  /*   MatHingPsiIintPsiTHingT2[i][6]=Mat[j*12+6];				  */
  /* }										  */
  /* MatHingPsiIintPsiTHingT2[6][6]=Mat[6*12+6];				  */
  /* 										  */
  /* invm(MatHingPsiIintPsiTHingT2);						  */
  /* for (i=0;i<7;++i)								  */
  /*   for (j=0;j<7;++j)							  */
  /*     InvMat[i*7+j]=MatHingPsiIintPsiTHingT2[i][j];				  */
  /* 										  */
  /* for (i=0;i<6;++i)								  */
  /*   theta[i]=T[i];								  */
  /* theta[6]=T[8];								  */
  /* 										  */
  /* mvmult(InvMat,theta,7);							  */
  /* //  for (i=0;i<12;++i) {							  */
  /* //    theta[i]=T[i];							  */
  /* //  }									  */
  /* 										  */
  /* if (TermMoveMode2 == 1) {							  */
  /*   for (i=0;i<6;++i) {							  */
  /*     acc_Term[i]=theta[i];							  */
  /*   }									  */
  /*   clust[1].dddihedang[0]=theta[6];						  */
  /* }										  */
  /* 										  */
  /* free(MatPsi);								  */
  /* free(MatInt);								  */
  /* free(Mat);									  */
  /* free(force);								  */
  /* free(biasf);								  */
  /* free(theta);								  */
  /* free(T);									  */
  /********************************************************************************/

}

void mmmult(double *m1, double *m2, double *m1m2, int n){
  int i,j,k;
  
  for (i=0;i<n;++i)
    for (j=0;j<n;++j)
      m1m2[i*n+j] = 0.0;

  for (i=0;i<n;++i)
    for (j=0;j<n;++j)
      for (k=0;k<n;++k)
	m1m2[i*n+j] += m1[i*n+k]*m2[k*n+j];

}

void mvmult(double *m, double *v, double *mv, int n){
  int i,j,k;

  
  for (i=0;i<n;++i)
      mv[i] = 0.0;

  for (i=0;i<n;++i)
    for (j=0;j<n;++j)
	mv[i] += m[i*n+j]*v[j];

}

void mtrans(double *m1,double *m2, int n){
  int i,j;


  for (i=0;i<n;++i)
    for (j=0;j<n;++j)
      m2[i*n+j] = m1[j*n+i];

}

void msetIni(double *m, int n) {
  int i,j;


  for (i=0;i<n;++i)
    for (j=0;j<n;++j)
      m[i*n+j] = 0.0;

  for (i=0;i<n;++i)
      m[i*n+i] = 1.0;

}

void msetzero(double *m, int n) {
  int i,j;


  for (i=0;i<n;++i)
    for (j=0;j<n;++j)
      m[i*n+j] = 0.0;

}

int invm(double mat[7][7],double inv[7][7]) {
  int i,j,k;
  double A[49];
  double test[7][7];
  static long int m=7,n=7,lda=7,info,piv[7],lwork=7;
  static double work[7];

  k=0;
  for(i=0;i<7;++i) {
    for(j=0;j<7;++j) {
      A[k]=mat[j][i];
      ++k;
    }
  }

  dgetrf_(&m,&n,A,&lda,piv,&info);

  dgetri_(&n,A,&lda,piv,work,&lwork,&info);

  k=0;
  for(i=0;i<7;++i) {
    for(j=0;j<7;++j) {
      inv[j][i]=A[k];
      ++k;
    }
  }

  for(i=0;i<7;++i)
    for(j=0;j<7;++j)
      test[i][j] = 0.0;

  for(i=0;i<7;++i)
    for(j=0;j<7;++j)
      for(k=0;k<7;++k)
	test[i][j]+=mat[i][k]*inv[k][j];

  return 0;

}

int invm6(double mat[6][6],double inv[6][6]) {
  int i,j,k;
  double A[49];
  double test[6][6];
  static long int m=6,n=6,lda=6,info,piv[6],lwork=6;
  static double work[6];

  k=0;
  for(i=0;i<6;++i) {
    for(j=0;j<6;++j) {
      A[k]=mat[j][i];
      ++k;
    }
  }

  dgetrf_(&m,&n,A,&lda,piv,&info);

  dgetri_(&n,A,&lda,piv,work,&lwork,&info);

  k=0;
  for(i=0;i<6;++i) {
    for(j=0;j<6;++j) {
      inv[j][i]=A[k];
      ++k;
    }
  }

  for(i=0;i<6;++i)
    for(j=0;j<6;++j)
      test[i][j] = 0.0;

  for(i=0;i<6;++i)
    for(j=0;j<6;++j)
      for(k=0;k<6;++k)
	test[i][j]+=mat[i][k]*inv[k][j];

  return 0;

}

int invm3(double mat[3][3],double inv[3][3]) {
  int i,j,k;
  double A[9];
  double test[3][3];
  static long int m=3,n=3,lda=3,info,piv[3],lwork=3;
  static double work[3];

  k=0;
  for(i=0;i<3;++i) {
    for(j=0;j<3;++j) {
      A[k]=mat[j][i];
      ++k;
    }
  }

  dgetrf_(&m,&n,A,&lda,piv,&info);

  dgetri_(&n,A,&lda,piv,work,&lwork,&info);

  k=0;
  for(i=0;i<3;++i) {
    for(j=0;j<3;++j) {
      inv[j][i]=A[k];
      ++k;
    }
  }

  for(i=0;i<3;++i)
    for(j=0;j<3;++j)
      test[i][j] = 0.0;

  for(i=0;i<3;++i)
    for(j=0;j<3;++j)
      for(k=0;k<3;++k)
	test[i][j]+=mat[i][k]*inv[k][j];

  return 0;

}
