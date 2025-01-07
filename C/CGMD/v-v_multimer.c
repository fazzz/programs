#include <stdio.h>
#include <math.h>

#include "gener.h"
#include "ABA_multimer.h"
#include "MD.h"
#include "physics.h"
#include "EF.h"

/*******************************/
/* double GearsConstant[6];    */
/* double GearsConstant2[6];   */
/* double GearsConstant3[5];   */
/* double Telar_Matrix[6][6];  */
/* double Telar_Matrix2[5][5]; */
/*******************************/
//double delta_s_NVT/*, acc_s_NVT*/;
//double predict_dihedang[6];

void Gear_six(int nNumClut);

// Velocity-Verlet 法を使った時間進展を行う関数
// 次のステップの 2 面角を計算
int velocity_verlet_step1(int nNumStep, int pflag,
			  double *ele, double *ALJ, double *BLJ,
			  double *p_e, double *p_1_4_e, 
			  double *p_LJ,double *p_1_4_LJ,
			  double *f_e, double *f_1_4_e, 
			  double *f_LJ,double *f_1_4_LJ,
			  int numnb, int *indexnb,
			  int num14, int *index14){
  int i,j,cflag;
  int nNumClt,nNumClut,nNumClutOrigBranch;
  FILE *outtest, *outtest2, *outtest3, *outtest4;
  
  double *deltathetat/*[MAXDOF]*/;
  double *dd_thetat_deltat/*[MAXDOF]*/;
  double *valueofforce/*[MAXDOF]*/;

  double /*delta_s_NVT,*/delta_xi_NVT/*,acc_s_NVT*/,inertiaofvertialvariable,inertiaofvertialvariable2;

  deltathetat=(double *)gcemalloc(sizeof(double)*prot.DOF);
  dd_thetat_deltat=(double *)gcemalloc(sizeof(double)*prot.DOF);
  valueofforce=(double *)gcemalloc(sizeof(double)*prot.DOF);

  // 次のステップの 2 面角加速度を計算
  /************************************************/
  /* if (GearMODE == OFF) {			  */
  /*   cflag=ON;				  */
  /*   calc_dd_theta_cycle(pflag,cflag,		  */
  /* 			ele, ALJ, BLJ,		  */
  /* 			p_e, p_1_4_e, 		  */
  /* 			p_LJ,p_1_4_LJ,		  */
  /* 			f_e,f_1_4_e, 		  */
  /* 			f_LJ,f_1_4_LJ,		  */
  /* 			numnb,indexnb,		  */
  /* 			num14,index14);		  */
  /* }						  */
  /************************************************/

  ///	if ((outtest=fopen("dddihedang.out","a")) == NULL)
  //	{
  //		printf("in\n");
  //		exit(1);
  //	}
  //	if ((outtest2=fopen("force2.out","a")) == NULL)
  //	{
  //		exit(1);
  //	}
  //	if ((outtest3=fopen("eata.out","a")) == NULL)
  //	{
  //		exit(1);
  //	}
  //	if ((outtest4=fopen("ddihedang.out","a")) == NULL)
  //	{
  //		exit(1);
  //	}
  
  //	for (nNumClt=1;nNumClt<prot.DOF;++nNumClt)
  //	{
  //		valueofforce[nNumClt] = clust[nNumClt].f_c.sp_f_clust[0].f_clust[0]*clust[nNumClt].f_c.sp_f_clust[0].f_clust[0]
  //		                 +clust[nNumClt].f_c.sp_f_clust[0].f_clust[1]*clust[nNumClt].f_c.sp_f_clust[0].f_clust[1]
  //		                 +clust[nNumClt].f_c.sp_f_clust[0].f_clust[2]*clust[nNumClt].f_c.sp_f_clust[0].f_clust[2];
  //	}

  GearsConstant[0] = 3.0/16.0;
  GearsConstant[1] = 251.0/360.0;
  GearsConstant[2] = 1.0;
  GearsConstant[3] = 11.0/18.0;
  GearsConstant[4] = 1.0/6.0;
  GearsConstant[5] = 1.0/60.0;
  
  GearsConstant2[0] = 3.0/20.0;
  GearsConstant2[1] = 251.0/360.0;
  GearsConstant2[2] = 1.0;
  GearsConstant2[3] = 11.0/18.0;
  GearsConstant2[4] = 1.0/6.0;
  GearsConstant2[5] = 1.0/60.0;
  
  GearsConstant3[0] = 251.0/720.0;
  GearsConstant3[1] = 1.0;
  GearsConstant3[2] = 11.0/12.0;
  GearsConstant3[3] = 1.0/3.0;
  GearsConstant3[4] = 1.0/24.0;
	
  for (i=0;i<6;++i){
    for (j=0;j<6;++j){
      if (i != j) {
	Telar_Matrix[i][j] = 0.0;
      }
      else{
	Telar_Matrix[i][i] = 1.0;
      }
    }
  }
	
  for (i=0;i<5;++i){
    for (j=0;j<5;++j){
      if (i != j){
	Telar_Matrix2[i][j] = 0.0;
      }
      else {
	Telar_Matrix2[i][i] = 1.0;
      }
    }
  }
	
  Telar_Matrix[0][1] = 1.0;
  Telar_Matrix[0][2] = 1.0;
  Telar_Matrix[0][3] = 1.0;
  Telar_Matrix[0][4] = 1.0;
  Telar_Matrix[0][5] = 1.0;
  Telar_Matrix[1][2] = 2.0;
  Telar_Matrix[1][3] = 3.0;
  Telar_Matrix[1][4] = 4.0;
  Telar_Matrix[1][5] = 5.0;
  Telar_Matrix[2][3] = 3.0;
  Telar_Matrix[2][4] = 6.0;
  Telar_Matrix[2][5] = 10.0;
  Telar_Matrix[3][4] = 4.0;
  Telar_Matrix[3][5] = 10.0;
  Telar_Matrix[4][5] = 5.0;
  
  Telar_Matrix2[0][1] = 1.0;
  Telar_Matrix2[0][2] = 1.0;
  Telar_Matrix2[0][3] = 1.0;
  Telar_Matrix2[0][4] = 1.0;
  Telar_Matrix2[1][2] = 2.0;
  Telar_Matrix2[1][3] = 3.0;
  Telar_Matrix2[1][4] = 4.0;
  Telar_Matrix2[2][3] = 3.0;
  Telar_Matrix2[2][4] = 6.0;
  Telar_Matrix2[3][4] = 4.0;


  if (MODE==NVT) {
    for (i=0;i<6;++i){
      s_NVT_predict[i] = 0.0;
    }
	  
    for (i=0;i<6;++i){
      for (j=0;j<6;++j){
	s_NVT_predict[i] += Telar_Matrix[i][j]*s_NVT_correct[j];
      }
    }

    s_NVT = s_NVT_predict[0];
    dot_s_NVT = s_NVT_predict[1]/deltat;

  }


  // 次のステップの 2 面角を計算
  if (ABI[0].hingmat!=6 ) {
    nNumClt= 1;
  }
  else {
    nNumClt = 0;
  }
  for(/*nNumClt = 1*/; nNumClt < prot.DOF; ++nNumClt) {
    if (GearMODE == ON)	{
      if (ABI[nNumClt].hingmat!=6)
	Gear(nNumClt);
      else 
	Gear_six(nNumClt);
    }
  }

  return nNumStep;
}

// Velocity-Verlet 法を使った時間進展
// 次のステップの 2 面角速度を計算
void velocity_verlet_step2(int pflag,
			   double *ele, double *ALJ, double *BLJ,
			   double *p_e, double *p_1_4_e, 
			   double *p_LJ,double *p_1_4_LJ,
			   double *f_e, double *f_1_4_e, 
			   double *f_LJ,double *f_1_4_LJ,
			   int numnb, int *indexnb,
			   int num14, int *index14,
			   double vel_Term[3],double *eig,double *eig_14,double *eig_dihed)

{
  int i,cflag;
  int nNumClt;
  FILE *outtest, *outtest2;
  double /*delta_s_NVT,*/delta_xi_NVT/*,acc_s_NVT*/,inertiaofvertialvariable,inertiaofvertialvariable2,kinetic_ene_t;
  //04_10
  int degF;
	
  GearsConstant[0] = 3.0/16.0;
  GearsConstant[1] = 251.0/360.0;
  GearsConstant[2] = 1.0;
  GearsConstant[3] = 11.0/18.0;
  GearsConstant[4] = 1.0/6.0;
  GearsConstant[5] = 1.0/60.0;
	
  GearsConstant2[0] = 3.0/20.0;
  GearsConstant2[1] = 251.0/360.0;
  GearsConstant2[2] = 1.0;
  GearsConstant2[3] = 11.0/18.0;
  GearsConstant2[4] = 1.0/6.0;
  GearsConstant2[5] = 1.0/60.0;
	
  GearsConstant3[0] = 251.0/720.0;
  GearsConstant3[1] = 1.0;
  GearsConstant3[2] = 11.0/12.0;
  GearsConstant3[3] = 1.0/3.0;
  GearsConstant3[4] = 1.0/24.0;

	

//	pre_dyn();
  if (MODE == NVT) {
    s_NVT = s_NVT_predict[0];
    dot_s_NVT = s_NVT_predict[1]/deltat;
    pre_dyn();

    CalcT(vel_Term);
    inertiaofvertialvariable2=2.0*Energy_kinetic8/q_NVT*s_NVT*(4.18407*100.0);

    if (TermMoveMode2==12)
      degF=prot.DOF-1+6;
    else
      degF=prot.DOF-1;

    acc_s_NVT = inertiaofvertialvariable2-/*(prot.DOF-1)*//*degF*/DOFOFPROT*k_B_kcm*/*1.38065e-23*/T_Kelvin/q_NVT*s_NVT
      *4.18407*100.0
      /*/1.660539e-27/1.0e-20*1.0e-24*/
      +dot_s_NVT*dot_s_NVT/s_NVT;

    /*************************************************************************/
    /* acc_s_NVT = inertiaofvertialvariable2-1.0/Energy_kinetic8/q_NVT*s_NVT */
    /*   +dot_s_NVT*dot_s_NVT/s_NVT;					     */
    /*************************************************************************/


  }
	// 次のステップの 2 面角を計算
  cflag = ON;
  //  printf("yes-v1\n");
  calc_dd_theta_cycle(pflag,cflag,
		      ele, ALJ, BLJ,
		      p_e, p_1_4_e, 
		      p_LJ,p_1_4_LJ,
		      f_e,f_1_4_e, 
		      f_LJ,f_1_4_LJ,
		      numnb,indexnb,
		      num14,index14,eig,eig_14,eig_dihed);
  //  printf("yes-v2\n");

  if (ABI[0].hingmat!=6 ) {
    nNumClt= 1;
  }
  else {
    nNumClt = 0;
  }
  for(/*nNumClt = 1*/; nNumClt < prot.DOF; ++nNumClt) {
    if (ABI[nNumClt].hingmat!=6)
      Gear2(nNumClt);
    else 
      Gear2_six(nNumClt);
  }
  if (MODE==NVT) {
    for (i=0;i<6;++i){
      s_NVT_correct[i] = 0.0;
    }

    for (i=0;i<6;++i){
      s_NVT_correct[i]
	= s_NVT_predict[i]+GearsConstant[i]*(0.5*deltat*deltat*acc_s_NVT-s_NVT_predict[2]);
    }
    s_NVT = s_NVT_correct[0];
    dot_s_NVT = s_NVT_correct[1]/deltat;

  }
  
}

void Gear(int nNumClut) {
  int i,j;
  
  for (i=0;i<6;++i) {
    clust[nNumClut].predict_dihedang[i] = 0.0;
  }

  clust[nNumClut].correct_dihedang[0] = clust[nNumClut].dihedang[0];
  
  for (i=0;i<6;++i) {
    for (j=0;j<6;++j) {
      clust[nNumClut].predict_dihedang[i] += Telar_Matrix[i][j]
	*clust[nNumClut].correct_dihedang[j];
    }
  }

  clust[nNumClut].now_deltadihedang[0] = 0.0;
  
  for (i=1;i<6;++i){
    //		clust[nNumClut].now_deltadihedang[0] += clust[nNumClut].predict_dihedang[i]
    //												/*clust[nNumClut].correct_dihedang[i]*/;
  }
  for (j=1;j<6;++j){
    clust[nNumClut].now_deltadihedang[0] += Telar_Matrix[0][j]
    *clust[nNumClut].correct_dihedang[j];
  }
  //		clust[nNumClut].now_deltadihedang[0] = clust[nNumClut].predict_dihedang[0]
  //											  -clust[nNumClut].correct_dihedang[0];
  
  clust[nNumClut].ddihedang[0] =   clust[nNumClut].predict_dihedang[1]/deltat;


  
}

void Gear2(int nNumClut) {
  int i;
  double delta_dddihedang;
  FILE *debug;

	//04_10

  GearsConstant[0] = 3.0/16.0;
  GearsConstant[1] = 251.0/360.0;
  GearsConstant[2] = 1.0;
  GearsConstant[3] = 11.0/18.0;
  GearsConstant[4] = 1.0/6.0;
  GearsConstant[5] = 1.0/60.0;
  
  GearsConstant2[0] = 3.0/20.0;
  GearsConstant2[1] = 251.0/360.0;
  GearsConstant2[2] = 1.0;
  GearsConstant2[3] = 11.0/18.0;
  GearsConstant2[4] = 1.0/6.0;
  GearsConstant2[5] = 1.0/60.0;
  
  GearsConstant3[0] = 251.0/720.0;
  GearsConstant3[1] = 1.0;
  GearsConstant3[2] = 11.0/12.0;
  GearsConstant3[3] = 1.0/3.0;
  GearsConstant3[4] = 1.0/24.0;

  
  delta_dddihedang =
    0.5*deltat*deltat*clust[nNumClut].dddihedang[0]
    -clust[nNumClut].predict_dihedang[2];
  
  for (i=0;i<6;++i){
    clust[nNumClut].correct_dihedang[i] = 0.0;
  }

//	clust[nNumClut].correct_dihedang[0] = clust[nNumClut].dihedang[0];
  
  for (i=0;i<6;++i){
    clust[nNumClut].correct_dihedang[i]
      = clust[nNumClut].predict_dihedang[i]+GearsConstant[i]*delta_dddihedang;
  }
  
  clust[nNumClut].now_deltadihedang[0] = 0.0;

  /*********************************************************************************/
  /* debug=fopen("debug_force.txt","a");					   */
  /* fprintf(debug,"%d %d %e \n",nNumStep,nNumClut,clust[nNumClut].dddihedang[0]); */
  /* fclose(debug);								   */
  /*********************************************************************************/

  //	clust[nNumClut].now_deltadihedang[0] = /*clust[nNumClut].correct_dihedang[0]-clust[nNumClut].predict_dihedang[0]*/
  //										   GearsConstant[0]*delta_dddihedang;
  //	clust[nNumClut].now_deltadihedang[0] = clust[nNumClut].correct_dihedang[0]-clust[nNumClut].predict_dihedang[0];
  //										   GearsConstant[0]*delta_dddihedang;
  clust[nNumClut].now_deltadihedang[0] = GearsConstant[0]*delta_dddihedang;
  clust[nNumClut].ddihedang[0] = clust[nNumClut].correct_dihedang[1]/deltat;

  //	if (MODE == NVT)
  //	{
  //		verlet_of_vartial_variable(nNumStep);
  //	}

  
}

void Gear_six(int nNumClut) {
  int i,j,k;
  
  for (i=0;i<6;++i) {
    for (j=0;j<6;++j) {
      clust[nNumClut].predict_dihedang_six[i][j] = 0.0;
    }
  }
  
  /***********************************************************************************/
  /* for (i=0;i<6;++i) {							     */
  /*   clust[nNumClut].correct_dihedang_six[i][0] = clust[nNumClut].dihedang_six[i]; */
  /* }										     */
  /***********************************************************************************/
    
  for (i=0;i<6;++i) {
    for (j=0;j<6;++j) {
      for (k=0;k<6;++k) {
	clust[nNumClut].predict_dihedang_six[i][j] += Telar_Matrix[j][k]
	  *clust[nNumClut].correct_dihedang_six[i][k];
      }
    }
  }

  for (i=0;i<6;++i)
    clust[nNumClut].now_deltadihedang_six[i] = 0.0;
  
  for (i=0;i<6;++i){
    for (j=1;j<6;++j){
      clust[nNumClut].now_deltadihedang_six[i] += Telar_Matrix[0][j]
	*clust[nNumClut].correct_dihedang_six[i][j];
    }
  }

  for (i=0;i<6;++i)
    clust[nNumClut].ddihedang_six[i] = clust[nNumClut].predict_dihedang_six[i][1]/deltat;
  
}

void Gear2_six(int nNumClut) {
  int i,j;
  double delta_dddihedang[6];
  FILE *debug;

	//04_10

  GearsConstant[0] = 3.0/16.0;
  GearsConstant[1] = 251.0/360.0;
  GearsConstant[2] = 1.0;
  GearsConstant[3] = 11.0/18.0;
  GearsConstant[4] = 1.0/6.0;
  GearsConstant[5] = 1.0/60.0;
  
  GearsConstant2[0] = 3.0/20.0;
  GearsConstant2[1] = 251.0/360.0;
  GearsConstant2[2] = 1.0;
  GearsConstant2[3] = 11.0/18.0;
  GearsConstant2[4] = 1.0/6.0;
  GearsConstant2[5] = 1.0/60.0;
  
  GearsConstant3[0] = 251.0/720.0;
  GearsConstant3[1] = 1.0;
  GearsConstant3[2] = 11.0/12.0;
  GearsConstant3[3] = 1.0/3.0;
  GearsConstant3[4] = 1.0/24.0;

  
  for (i=0;i<6;++i)
    delta_dddihedang[i] =
      0.5*deltat*deltat*clust[nNumClut].dddihedang_six[i]
      -clust[nNumClut].predict_dihedang_six[i][2];
  
  for (i=0;i<6;++i){
    for (j=0;j<6;++j){
      clust[nNumClut].correct_dihedang_six[i][j] = 0.0;
    }
  }

//	clust[nNumClut].correct_dihedang[0] = clust[nNumClut].dihedang[0];
  
  for (i=0;i<6;++i){
    for (j=0;j<6;++j){
      clust[nNumClut].correct_dihedang_six[i][j]
	= clust[nNumClut].predict_dihedang_six[i][j]+GearsConstant[j]*delta_dddihedang[i];
    }
  }
  
  //  clust[nNumClut].now_deltadihedang[0] = 0.0;

  for (i=0;i<6;++i){
    clust[nNumClut].now_deltadihedang_six[i] = GearsConstant[0]*delta_dddihedang[i];
  }

  /*********************************************************************************/
  /* debug=fopen("debug_force.txt","a");					   */
  /* fprintf(debug,"%d %d %e \n",nNumStep,nNumClut,clust[nNumClut].dddihedang[0]); */
  /* fclose(debug);								   */
  /*********************************************************************************/

  //	clust[nNumClut].now_deltadihedang[0] = /*clust[nNumClut].correct_dihedang[0]-clust[nNumClut].predict_dihedang[0]*/
  //										   GearsConstant[0]*delta_dddihedang;
  //	clust[nNumClut].now_deltadihedang[0] = clust[nNumClut].correct_dihedang[0]-clust[nNumClut].predict_dihedang[0];
  //										   GearsConstant[0]*delta_dddihedang;
  //  clust[nNumClut].now_deltadihedang[0] = GearsConstant[0]*delta_dddihedang;
  for (i=0;i<6;++i)
    clust[nNumClut].ddihedang_six[i] = clust[nNumClut].correct_dihedang_six[i][1]/deltat;

  //	if (MODE == NVT)
  //	{
  //		verlet_of_vartial_variable(nNumStep);
  //	}

  
}

int verlet_of_vartial_variable(int nNumStep)
{
////	double xi_dummy;
//////	dot_dot_s_NVT = 1.0/(s_NVT*tau_NVT*deltat*tau_NVT*deltat)*(T_Kelvin_Now/T_Kelvin-1.0-1.0/(prot.DOF-1));
//////	dot_s_NVT += deltat*0.5*(dot_dot_s_NVT+dot_dot_s_NVT_old);
//////	s_NVT += deltat*dot_s_NVT+deltat*deltat*0.5*;
//////	dot_dot_s_NVT_old = dot_dot_s_NVT;
//	dot_dot_s_NVT = s_NVT/(tau_NVT*deltat*tau_NVT*deltat)*(T_Kelvin_Now/T_Kelvin-1.0-1.0/(prot.DOF-1))+dot_s_NVT/s_NVT;
//	dot_s_NVT += deltat*0.5*(dot_dot_s_NVT+dot_dot_s_NVT_old);
//	s_NVT += deltat*dot_s_NVT+deltat*deltat*0.5*dot_dot_s_NVT;
//	dot_dot_s_NVT_old = dot_dot_s_NVT;
//	xi_dummy = dot_s_NVT/s_NVT;
//	xi_NVT += deltat/(tau_NVT*deltat*tau_NVT*deltat)*(T_Kelvin_Now/T_Kelvin-1.0-1.0/(prot.DOF-1));
//////	s_NVT += deltat*xi_NVT/*+deltat*deltat/(tau_NVT*deltat*tau_NVT*deltat)*(T_Kelvin_Now/T_Kelvin-1.0)*/;

  return nNumStep;
}

/****************************************************************************************************/
/* int modified_velocity_verlet_1(int nNumStep, int pflag,					    */
/* 			       double *ele, double *ALJ, double *BLJ,				    */
/* 			       double *p_e, double *p_1_4_e, 					    */
/* 			       double *p_LJ,double *p_1_4_LJ,					    */
/* 			       double *f_e, double *f_1_4_e, 					    */
/* 			       double *f_LJ,double *f_1_4_LJ,					    */
/* 			       int numnb, int *indexnb,						    */
/* 			       int num14, int *index14,						    */
/* 			       int numDOF) {							    */
/*   int nNumClut;										    */
/* 												    */
/*   for (nNumClut=0;nNumClut<numDOF;++nNumClut) {						    */
/*     clust[nNumClut].ddihedang[0] = 1.5*vel_bh[nNumClut] - 0.5*vel_boh[nNumClut];		    */
/*   }												    */
/* 												    */
/*   pre_dyn();											    */
/*   calc_dd_theta_cycle(pflag,	    								    */
/* 		      ele, ALJ, BLJ,								    */
/* 		      p_e, p_1_4_e, 								    */
/* 		      p_LJ,p_1_4_LJ,								    */
/* 		      f_e,f_1_4_e,								    */
/* 		      f_LJ,f_1_4_LJ,								    */
/* 		      numnb,indexnb,								    */
/* 		      num14,index14);								    */
/* 												    */
/*   for (nNumClut=0;nNumClut<numDOF;++nNumClut) {						    */
/*     clust[nNumClut].ddihedang[0] = vel_bh[nNumClut]+deltat*clust[nNumClut].dddihedang[0];	    */
/*     vel_bh[nNumClut] += deltat*clust[nNumClut].dddihedang[0];				    */
/*     clust[nNumClut].now_deltadihedang[0] = 0.5*vel_bh[nNumClut];				    */
/*   }												    */
/* 												    */
/* }												    */
/****************************************************************************************************/

