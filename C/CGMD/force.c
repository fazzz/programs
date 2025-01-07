#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "MD.h"
#include "force.h"
#include "BD.h"
#include "EF.h"

//void calc_dot_Pseduo_TransMatrix(int nNumClutI,
//						         int nNumClutK,
//						         double mat[4][4]);

// 現在の構造でのタンパク質のポテンシャルエネルギー、
// タンパク質に及ぼす力の計算を行う関数
void calc_force(int pflag)
{
	int i,ii,i_c,alpha,alpha2,alpha3,j,k;

	int nNumClut;
	int nNumClutI;
	int nNumClutJ;
	int nNumAtom;

	int nNumOrigc;
	int nNumAbsoc;

	int origin_of_his_branch = 0;

	int nNumClutOfParent;
	int nNmAtomOfCN_1_A;

	double **q_c/*[MAXNAC][3]*/;

	double qq[3], qq2[3];

	double f_clust[100][4];
	double N_clust[100][3];
	double f_i_c[3];
	double f_vect[400];

	double RotatnNumtonNumMiOn[3][3];

	double mat[4][4];
	double dummy[4][4];
	double dummy_xoord[4];
	double dummy_xoord2[4];
	double dotTransQOne[10][20][4];
	double dotTransQTwo[10][20][4];
	double dotTransQ[10][20][4];
	double dotTransQ_mat[10][40];

	FILE *out;

	q_c=(double **)gcemalloc(sizeof(double *)*prot.num_atom);
	for (i=0;i<prot.num_atom;++i)
	  q_c[i]=(double *)gcemalloc(sizeof(double)*3);


	// spatial_force の初期化を行う
	for(nNumClut = 0; nNumClut < prot.DOF; ++nNumClut)
	{
		for(alpha=0;alpha<3;++alpha)
		{
			clust[nNumClut].f_c.sp_f_clust[0].f_clust[alpha] = 0.0;
			clust[nNumClut].f_c.sp_f_clust[0].N_clust[alpha] = 0.0;
			f_clust[nNumClut][alpha] = 0.0;
			N_clust[nNumClut][alpha] = 0.0;
		}
	}

	for(nNumClut = 0;nNumClut < prot.DOF; ++nNumClut)
	{
		clust[nNumClut].f_c.f_dihed = 0.0;
		for(i=0;i < clust[nNumClut].num_atom_clust;++i)
		{
			for(alpha=0;alpha < 3; ++alpha)
			{
				clust[nNumClut].f_c.f_elesta[i][alpha] = 0.0;
				clust[nNumClut].f_c.f_L_J[i][alpha] = 0.0;
				clust[nNumClut].f_c.f_1_4_elesta[i][alpha] = 0.0;
				clust[nNumClut].f_c.f_1_4_L_J[i][alpha] = 0.0;
			}
		}
	}

	/*******************************/
        /* Calc_vector_atoms();	       */
        /*******************************/

	// 2 面角相互作用の計算
//	Calc_dihed_Potential_Force();


	// VDW 相互作用の計算
	if (pflag!=AMBERMODE) {
	  Calc_L_J_PotentialandForce();
	}
	else {
	  //	  Calc_L_J_PotentialandForce2();
	}

	if (restflag==ON) {
	  // 制限をかける
	  //	  Calc_restraintForce();
	}
	// 静電相互作用の計算
//	Calc_ele_sta_PotentialandForce();

//	for (nNumAtom=0;nNumAtom<prot.num_atom;++nNumAtom)
//	{
//		for (alpha=0;alpha<4;++alpha)
//		{
//			f_clust[nNumAtom][alpha] = 0.0;
//		}
//	}

//	nNumAtom=0;
	for (nNumClut=/*1*/0;nNumClut</*=*/prot.DOF/*-1*/;++nNumClut)
	{
		for(i_c=0;i_c < clust[nNumClut].num_atom_clust;++i_c)
		{
			for (alpha=0;alpha<3;++alpha)
			{
				f_clust[nNumClut/*nNumAtom*/][alpha]
			              +=
							   clust[nNumClut].f_c.f_L_J[i_c][alpha]
						 	 + clust[nNumClut].f_c.f_elesta[i_c][alpha]
							 + clust[nNumClut].f_c.f_1_4_L_J[i_c][alpha]
						 	 + clust[nNumClut].f_c.f_1_4_elesta[i_c][alpha]
							;
			}
			f_clust[nNumClut/*nNumAtom*/][3] = 0.0;
//			++nNumAtom;
		}
	}
/////////////////////////////////////////////////////////////////////////////////////////

	for (nNumClut=1;nNumClut<=prot.DOF-1;++nNumClut)
	{
		for (nNumClutJ=1;nNumClutJ<=prot.DOF-1;++nNumClutJ)
//		for (nNumAtom=0;nNumAtom<prot.num_atom;++nNumAtom)
		{
			for (i=0;i<4;++i)
			{
				dotTransQ[nNumClut][nNumClutJ/*nNumAtom*/][i] = 0.0;
			}
		}
	}

//	nNumClutOfParent = clust[nNumClut].nNumClutOfParent-1;

	for (nNumClutI=1;nNumClutI<=prot.DOF-1;++nNumClutI)
	{
		for (nNumClutJ=1;nNumClutJ<=nNumClutI;++nNumClutJ)
		{
			for (i=0;i<4;++i)
			{
				for (j=0;j<4;++j)
				{
					mat[i][j] = 0.0;
				}
			}
//			calc_dot_Pseduo_TransMatrix(nNumClutI,nNumClutJ,mat);
//			calc_dot_Pseduo_TransMatrix(nNumClutJ,nNumClutI,mat);
//			calc_dot_Pseduo_TransMatrix(nNumClutI,nNumClutI,mat);
			calc_dot_Pseduo_TransMatrix(nNumClutJ,nNumClutJ,mat);

//			for (nNumAtom=clust[nNumClutI].origin_xoord_a+1;
//				 nNumAtom<clust[nNumClutI].num_xoord_a;
//				 ++nNumAtom)
//			{
				for (i=0;i<3;++i)
				{
//					dummy_xoord[i] = clust[nNumClutI].xoord_clust/*[0]*/[nNumAtom][i];
//					dummy_xoord[i] = clust[nNumClutJ].qCOM[i];
					dummy_xoord[i] = clust[nNumClutI].qCOM[i];
				}
				dummy_xoord[3] = /*0.0*/1.0;

				if (nNumClutI != nNumClutJ)
				{
//					for(nNumClut=nNumClutJ+1;nNumClut<=nNumClutI;++nNumClut)
					for(nNumClut=nNumClutI;nNumClut>nNumClutJ;--nNumClut)
					{
						for (i=0;i<4;++i)
						{
							dummy_xoord2[i] = dummy_xoord[i];
						}

						for (i=0;i<4;++i)
						{
							dummy_xoord[i] = 0.0;
						}

						for (i=0;i<4;++i)
						{
							for (j=0;j<4;++j)
							{
								dummy_xoord[i] += clust[nNumClut].PsedoTransMatrix[i][j]*dummy_xoord2[j];
							}
						}
					}
				}

					for (i=0;i<4;++i)
					{
						for (j=0;j<4;++j)
						{
//							nNmAtomOfCN_1_A = prot.num_atom - clust[nNumClutI].num_xoord_a;
							dotTransQ/*Two*/[nNumClutI][nNumClutJ/*nNumAtom+nNmAtomOfCN_1_A*/][i]
											+= mat[i][j]
											  *dummy_xoord[j];
						}
					}
//			}
		}
	}

/////////////////////////////////////////////////////////////////////////////////////////
	for(nNumClut=1;nNumClut<=prot.DOF-1;++nNumClut)
	{
		Patrix[nNumClut] = 0.0;
	}

	nNumAtom = 0;
	for(nNumClut=1;nNumClut<=prot.DOF-1;++nNumClut)
	{
		for(alpha=0;alpha < /*4*/3;++alpha)
		{
//			f_vect[nNumAtom] = 0.0;
			++nNumAtom;
			f_vect[nNumAtom] = f_clust[nNumClut][alpha];
		}
	}

	for(nNumClutJ=1;nNumClutJ<=prot.DOF-1;++nNumClutJ)
	{
		nNumAtom = 0;
		for (nNumClut=1;nNumClut<=prot.DOF-1;++nNumClut)
		{
			for(alpha=0;alpha < /*4*/3;++alpha)
			{
//				dotTransQ_mat[nNumClut][nNumAtom] = 0.0;
				++nNumAtom;
//				dotTransQ_mat[nNumClut][nNumAtom] = dotTransQ[nNumClut][nNumClutJ][alpha];
				dotTransQ_mat[nNumClutJ][nNumAtom] = dotTransQ[nNumClut][nNumClutJ][alpha];
			}
		}
	}

	for(nNumClut=1;nNumClut<=prot.DOF-1;++nNumClut)
	{
		for (nNumAtom=1;nNumAtom<=(prot.DOF-1)/**4*/*3;++nNumAtom)
		{
			Patrix[nNumClut] += f_vect[nNumAtom]
								   *dotTransQ_mat[nNumClut][nNumAtom];
		}
	}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
	for(nNumClut=1;nNumClut<=prot.DOF-1;++nNumClut)
	{
		clust[nNumClut].f_c.f_dihed = 0.0;
	}
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

	for(nNumClut=1;nNumClut<=prot.DOF-1;++nNumClut)
	{
		clust[nNumClut].f_c.f_dihed /*-*/+= Patrix[nNumClut];
//		clust[nNumClut].f_c.f_dihed = -1.0*clust[nNumClut].f_c.f_dihed;
//		clust[nNumClut].f_c.f_dihed /*+*/-= Patrix[nNumClut];
	}
}

//void calc_dot_Pseduo_TransMatrix(int nNumClutI,
//						         int nNumClutK,
//						         double mat[4][4])
//{
//	int i,j,k;
//	int nNumClut;
//	double dummy_matrix[4][4], dummy_matrix2[4][4];
//
//	for (i=0;i<4;++i)
//	{
//		for (j=0;j<4;++j)
//		{
//			dummy_matrix[i][j] = ident[i][j];
//		}
//	}
//
//	for (nNumClut=1;nNumClut<=nNumClutI;++nNumClut)
//	{
//		for (i=0;i<4;++i)
//		{
//			for (j=0;j<4;++j)
//			{
//				dummy_matrix2[i][j] = dummy_matrix[i][j];
//			}
//		}
//
//		for (i=0;i<4;++i)
//		{
//			for (j=0;j<4;++j)
//			{
//				dummy_matrix[i][j] = 0.0;
//			}
//		}
//
//		for (i=0;i<4;++i)
//		{
//			for (j=0;j<4;++j)
//			{
//				for (k=0;k<4;++k)
//				{
//					dummy_matrix[i][j] += dummy_matrix2[i][k]
//										 *clust[nNumClut].PsedoTransMatrix[k][j];
//				}
//			}
//		}
//	}
//
//	for (i=0;i<4;++i)
//	{
//		for (j=0;j<4;++j)
//		{
//			dummy_matrix2[i][j] = dummy_matrix[i][j];
//		}
//	}
//
//	for (i=0;i<4;++i)
//	{
//		for (j=0;j<4;++j)
//		{
//			dummy_matrix[i][j] = 0.0;
//		}
//	}
//
//	for (i=0;i<4;++i)
//	{
//		for (j=0;j<4;++j)
//		{
//			for (k=0;k<4;++k)
//			{
//				dummy_matrix[i][j] += dummy_matrix2[i][k]
//									  *delta_matrix[k][j];
//			}
//		}
//	}
//
//	for (nNumClut=nNumClutI+1;nNumClut<=nNumClutK;++nNumClut)
//	{
//		for (i=0;i<4;++i)
//		{
//			for (j=0;j<4;++j)
//			{
//				dummy_matrix2[i][j] = dummy_matrix[i][j];
//			}
//		}
//
//		for (i=0;i<4;++i)
//		{
//			for (j=0;j<4;++j)
//			{
//				dummy_matrix[i][j] = 0.0;
//			}
//		}
//
//		for (i=0;i<4;++i)
//		{
//			for (j=0;j<4;++j)
//			{
//				for (k=0;k<4;++k)
//				{
//					dummy_matrix[i][j] += dummy_matrix2[i][k]
//										 *clust[nNumClut].PsedoTransMatrix[k][j];
//				}
//			}
//		}
//	}
//
//	for (i=0;i<4;++i)
//	{
//		for (j=0;j<4;++j)
//		{
//			mat[i][j] = dummy_matrix[i][j];
//		}
//	}
//}

// 現在の構造でのタンパク質のポテンシャルエネルギー、
// タンパク質に及ぼす力の計算を行う関数
void calc_force2(int pflag,
		 double *ele, double *ALJ, double *BLJ,
		 double *p_e, double *p_1_4_e, 
		 double *p_LJ,double *p_1_4_LJ,
		 double *f_e, double *f_1_4_e, 
		 double *f_LJ,double *f_1_4_LJ,
		 int numnb, int *indexnb,
		 int num14, int *index14,double *eigen,double *eig_14,double *eig_dihed)
{
  int i,ii,i_c,alpha,alpha2,alpha3,j;

	int num;
	int nNumClut;

	int nNumOrigc;
	int nNumAbsoc;

	int origin_of_this_branch = 0;

	double **q_c/*[MAXNAC][3]*/;

	double qq[3], qq2[3];

	double **f_clust/*[MAXDOF][3]*/;
	double **N_clust/*[MAXDOF][3]*/;
	double ***f_clust3/*[MAXDOF][20][3]*/;
	double ***f_clust4/*[MAXDOF][20][3]*/;
	double f_i_c[3];

	double RotatnNumtonNumMiOn[3][3];

	double *cord;

	FILE *out;

	int num_a_prot;
	int tnum_atom_clust,tDOF=prot.DOF;

	q_c=(double **)gcemalloc(sizeof(double *)*prot.num_atom);
	for (i=0;i<prot.num_atom;++i)
	  q_c[i]=(double *)gcemalloc(sizeof(double)*3);
	f_clust=(double **)gcemalloc(sizeof(double *)*prot.DOF);
	N_clust=(double **)gcemalloc(sizeof(double *)*prot.DOF);
	for (i=0;i<prot.DOF;++i) {
	  f_clust[i]=(double *)gcemalloc(sizeof(double)*3);
	  N_clust[i]=(double *)gcemalloc(sizeof(double)*3);
	}
	f_clust3=(double **)gcemalloc(sizeof(double *)*prot.DOF);
	f_clust4=(double **)gcemalloc(sizeof(double *)*prot.DOF);
	for (i=0;i<prot.DOF;++i) {
	  f_clust3[i]=(double **)gcemalloc(sizeof(double *)*20);
	  f_clust4[i]=(double **)gcemalloc(sizeof(double *)*20);
	  for (j=0;j<20;++j) {
	  f_clust3[i][j]=(double *)gcemalloc(sizeof(double)*20);
	  f_clust4[i][j]=(double *)gcemalloc(sizeof(double)*20);
	  }
	}


	// spatial_force の初期化を行う
	for(nNumClut = 0; nNumClut < prot.DOF; ++nNumClut)
	{
		for(alpha=0;alpha<3;++alpha)
		{
			clust[nNumClut].f_c.sp_f_clust[0].f_clust[alpha] = 0.0;
			clust[nNumClut].f_c.sp_f_clust[0].N_clust[alpha] = 0.0;
			f_clust[nNumClut][alpha] = 0.0;
			N_clust[nNumClut][alpha] = 0.0;
		}
	}

	for(nNumClut = 0;nNumClut < prot.DOF; ++nNumClut)
	{
		for(i=0;i < clust[nNumClut].num_atom_clust;++i)
		{
			for(alpha=0;alpha < 3; ++alpha)
			{
				clust[nNumClut].f_c.f_elesta[i][alpha] = 0.0;
				clust[nNumClut].f_c.f_L_J[i][alpha] = 0.0;
				clust[nNumClut].f_c.f_1_4_elesta[i][alpha] = 0.0;
				clust[nNumClut].f_c.f_1_4_L_J[i][alpha] = 0.0;
				f_clust3[nNumClut][i][alpha] = 0.0;
				f_clust4[nNumClut][i][alpha] = 0.0;
			}
		}
	}

	/*******************************/
        /* Calc_vector_atoms();	       */
        /*******************************/

	if (dhstopflag==2) {
	  Calc_dihed_Potential_Force_for_db();
	}
	else if (dhstopflag!=ON) {
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
	  // 2 面角相互作用の計算
	  //  printf("yes-f1\n");
	  Calc_dihed_Potential_Force(eig_dihed);
	  //  printf("yes-f2\n");
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
	}
	
	if(nbstopflag!=ON) {
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
	// VDW 相互作用の計算
	if (pflag!=AMBERMODE) {
	  Calc_L_J_PotentialandForce();
	}
	else {
	  if ((cord=(double *)malloc(sizeof(double)*prot.num_atom*3))==NULL) {
	    printf("error of allocation cord\n");
	  }
	  for (i=0;i<prot.num_atom;++i) {
	    for (alpha=0;alpha<3;++alpha) {
	      cord[i*3+alpha]=prot.coord[i][alpha];
	    }
	  }
	  //  printf("yes-f3\n");
	  Calc_L_J_PotentialandForce2(ele, ALJ, BLJ,
				      p_e,  p_1_4_e, 
				      p_LJ, p_1_4_LJ,
				      f_e, f_1_4_e, 
				      f_LJ,f_1_4_LJ,
				      numnb, indexnb,
				      num14, index14, cord,eigen,eig_14);
	  //  printf("yes-f4\n");

	  /******************************************************************************************/
          /* num_a_prot=0;									    */
	  /* for(nNumClut = 0;nNumClut < tDOF; ++nNumClut) {					    */
	  /*   tnum_atom_clust = clust[nNumClut].num_atom_clust;				    */
	  /*   for(i=0;i < tnum_atom_clust;++i) {						    */
	  /*     potential_pro.p_elesta[num_a_prot] = p_e[num_a_prot];				    */
	  /*     potential_pro.p_1_4_elesta[num_a_prot] = p_1_4_e[num_a_prot];			    */
	  /*     potential_pro.p_L_J[num_a_prot] = p_LJ[num_a_prot];				    */
	  /*     potential_pro.p_1_4_L_J[num_a_prot] = p_1_4_LJ[num_a_prot];			    */
	  /*     for(alpha=0; alpha<3; ++alpha) {						    */
	  /* 	clust[nNumClut].f_c.f_L_J[i][alpha] = f_LJ[num_a_prot*3+alpha];			    */
	  /* 	clust[nNumClut].f_c.f_1_4_L_J[i][alpha] = f_1_4_LJ[num_a_prot*3+alpha];		    */
	  /* 	clust[nNumClut].f_c.f_elesta[i][alpha] = f_e[num_a_prot*3+alpha];		    */
	  /* 	clust[nNumClut].f_c.f_1_4_elesta[i][alpha] = f_1_4_e[num_a_prot*3+alpha];	    */
	  /*     }										    */
	  /*     ++num_a_prot;									    */
	  /*   }										    */
	  /* }											    */
          /******************************************************************************************/
	  free(cord);

	}
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
	}
	if (restflag==ON) {
	  // 制限をかける
	  Calc_restraintforce();
	}


	// 静電相互作用の計算
//	Calc_ele_sta_PotentialandForce();


	for(nNumClut=0;nNumClut<prot.DOF;++nNumClut)
	{
		for(i_c=0;i_c<clust[nNumClut].num_atom_clust;++i_c)
		{
			for(alpha=0;alpha < 3;++alpha)
			{
				f_clust[nNumClut][alpha]
				              +=
								   clust[nNumClut].f_c.f_L_J[i_c][alpha]
							 	 + clust[nNumClut].f_c.f_elesta[i_c][alpha]
								 + clust[nNumClut].f_c.f_1_4_L_J[i_c][alpha]
							 	 + clust[nNumClut].f_c.f_1_4_elesta[i_c][alpha]
								;
							  
				f_clust3[nNumClut][i_c][alpha]
				              =
								   clust[nNumClut].f_c.f_L_J[i_c][alpha]
							 	 + clust[nNumClut].f_c.f_elesta[i_c][alpha]
								 + clust[nNumClut].f_c.f_1_4_L_J[i_c][alpha]
							 	 + clust[nNumClut].f_c.f_1_4_elesta[i_c][alpha];

			}
		}
	}

//	out = fopen("force.out","a");
	for(nNumClut=0;nNumClut<prot.DOF;++nNumClut)
	{
		for (alpha=0;alpha<3;++alpha)
		{
			for (alpha2=0;alpha2<3;++alpha2)
			{
					// Rot_0_to_N * Rot_n-1_to_0
					clust[nNumClut].f_c.sp_f_clust[0].f_clust[alpha]
					+= 	clust[nNumClut].trans_A_to_CN[0][alpha][alpha2]
				       *f_clust[nNumClut][alpha2];
			}
		}
	}

	for(nNumClut=0;nNumClut<prot.DOF;++nNumClut)
	{
		for(i_c=0;i_c<clust[nNumClut].num_atom_clust;++i_c)
		{
			for (alpha=0;alpha<3;++alpha)
			{
				for (alpha2=0;alpha2<3;++alpha2)
				{
						// Rot_0_to_N * Rot_n-1_to_0
						f_clust4[nNumClut][i_c][alpha]
						+= 	clust[nNumClut].trans_A_to_CN[0][alpha][alpha2]
					       *f_clust3[nNumClut][i_c][alpha2];
				}
			}
		}
	}

	for(nNumClut=0;nNumClut<prot.DOF;++nNumClut)
	{
	  num = 0/*0410*//*clust[nNumClut].origin_xoord_a*/;
		for(i_c=0;i_c<clust[nNumClut].num_atom_clust;++i_c)
		{
			clust[nNumClut].f_c.sp_f_clust[0].N_clust[0]
			+= -clust[nNumClut].xoord_clust/*[0]*/[i_c+num][2]*f_clust4[nNumClut][i_c][1]
			   +clust[nNumClut].xoord_clust/*[0]*/[i_c+num][1]*f_clust4[nNumClut][i_c][2];

			clust[nNumClut].f_c.sp_f_clust[0].N_clust[1] 
			+= clust[nNumClut].xoord_clust/*[0]*/[i_c+num][2]*f_clust4[nNumClut][i_c][0]
			  -clust[nNumClut].xoord_clust/*[0]*/[i_c+num][0]*f_clust4[nNumClut][i_c][2];

			clust[nNumClut].f_c.sp_f_clust[0].N_clust[2]
			+= -clust[nNumClut].xoord_clust/*[0]*/[i_c+num][1]*f_clust4[nNumClut][i_c][0]
			   +clust[nNumClut].xoord_clust/*[0]*/[i_c+num][0]*f_clust4[nNumClut][i_c][1];

		}
	}

//	for(nNumClut=0;nNumClut<prot.DOF;++nNumClut)
//	{
//		nNumOrigc = clust[nNumClut].origin_atom_a - 1;
//
//		for(i_c=0;i_c < clust[nNumClut].num_atom_clust;++i_c)
//		{
//			nNumAbsoc = clust[nNumClut].origin_atom_a - 1 + i_c;
//
//			for(alpha=0;alpha<3;++alpha)
//			{
///*(7_19)*/		qq/*[i_c]*/[alpha] = prot.coord[nNumAbsoc][alpha]
//				                   -prot.coord[nNumOrigc][alpha]; /* A */
//			}
//
//			for (alpha=0;alpha<3;++alpha)
//			{
//				qq2[alpha] = 0.0;
//			}
//
//			for (alpha=0;alpha<3;++alpha)
//			{
//				for (alpha2=0;alpha2<3;++alpha2)
//				{
//					// Rot_0_to_N * Rot_n-1_to_0
//					qq2[alpha]
//						+= clust[nNumClut].trans_A_to_CN[0][alpha][alpha2]
//					      *qq[alpha2];
//				}
//			}
//
//			for (alpha=0;alpha<3;++alpha)
//			{
//				f_i_c[alpha] = 0.0;
//			}
//
//			for (alpha=0;alpha<3;++alpha)
//			{
//				for (alpha2=0;alpha2<3;++alpha2)
//				{
//					// Rot_0_to_N * Rot_n-1_to_0
//					f_i_c[alpha]
//						+= clust[nNumClut].trans_A_to_CN[0][alpha][alpha2]
//					  *(
//					  	 clust[nNumClut].f_c.f_L_J[i_c][alpha2]
//					    +clust[nNumClut].f_c.f_elesta[i_c][alpha2]
//					    +clust[nNumClut].f_c.f_1_4_L_J[i_c][alpha2]
//					    +clust[nNumClut].f_c.f_1_4_elesta[i_c][alpha2]
//					  	)
//						;
//				}
//			}
//
//			clust[nNumClut].f_c.sp_f_clust[0].N_clust[0] /* N*m */
//			+= 	/*(*/  -(qq2/*_c[i_c]*/[2]/* A *//**1.0e-10*//* m/A */)
//				  *f_i_c[1]/* N */
//			      +(qq2/*_c[i_c]*/[1]/* A *//**1.0e-10*//* m/A */)
//				  *f_i_c[2]/* N */
//				/*)*2.0*/;
//
//			clust[nNumClut].f_c.sp_f_clust[0].N_clust[1] /* N*m */
//			+= 	/*(*/   (qq2/*_c[i_c]*/[2]/* A *//**1.0e-10*//* m/A */)
//				   *f_i_c[0]/* N */
//			      -(qq2/*_c[i_c]*/[0]/* A *//**1.0e-10*//* m/A */)
//				   *f_i_c[2]/* N */
//				/*)*2.0*/;
//
//			clust[nNumClut].f_c.sp_f_clust[0].N_clust[2] /* N*m */ 
//			+= 	/*(*/   -(qq2/*_c[i_c]*/[1]/* A *//**1.0e-10*//* m/A */)
//				   *f_i_c[0]/* N */
//			       +(qq2/*_c[i_c]*/[0]/* A *//**1.0e-10*//* m/A */)
//				   *f_i_c[1]/* N */
//				/*)*2.0*/;

//	for(nNumClut=0;nNumClut<prot.DOF;++nNumClut)
//	{
//			clust[nNumClut].f_c.sp_f_clust[0].N_clust[0] /* N*m */
//			+= 	/*(*/  -(clust[nNumClut].qCOM/*_c[i_c]*/[2]/* A *//**1.0e-10*//* m/A */)
//				  *clust[nNumClut].f_c.sp_f_clust[0].f_clust[1]/* N */
//			      +(clust[nNumClut].qCOM/*_c[i_c]*/[1]/* A *//**1.0e-10*//* m/A */)
//				  *clust[nNumClut].f_c.sp_f_clust[0].f_clust[2]/* N */
//				/*)*2.0*/;
//
//			clust[nNumClut].f_c.sp_f_clust[0].N_clust[1] /* N*m */
//			+= 	/*(*/   (clust[nNumClut].qCOM/*_c[i_c]*/[2]/* A *//**1.0e-10*//* m/A */)
//				   *clust[nNumClut].f_c.sp_f_clust[0].f_clust[0]/* N */
//			      -(clust[nNumClut].qCOM/*_c[i_c]*/[0]/* A *//**1.0e-10*//* m/A */)
//				   *clust[nNumClut].f_c.sp_f_clust[0].f_clust[2]/* N */
//				/*)*2.0*/;
//
//			clust[nNumClut].f_c.sp_f_clust[0].N_clust[2] /* N*m */ 
//			+= 	/*(*/   -(clust[nNumClut].qCOM/*_c[i_c]*/[1]/* A *//**1.0e-10*//* m/A */)
//				   *clust[nNumClut].f_c.sp_f_clust[0].f_clust[0]/* N */
//			       +(clust[nNumClut].qCOM/*_c[i_c]*/[0]/* A *//**1.0e-10*//* m/A */)
//				   *clust[nNumClut].f_c.sp_f_clust[0].f_clust[1]/* N */
//				/*)*2.0*/;

//			N_clust[nNumClut][0] /* N*m */
//			+= 	  -(q/*_c[i_c]*/[2]/* A */*1.0e-10/* m/A */)
//				  *f_i_c[1]/* N */
//			      +(q/*_c[i_c]*/[1]/* A */*1.0e-10/* m/A */)
//				  *f_i_c[2]/* N */;

//			N_clust[nNumClut][1] /* N*m */
//			+= 	   (q/*_c[i_c]*/[2]/* A */*1.0e-10/* m/A */)
//				   *f_i_c[0]/* N */
//			      -(q/*_c[i_c]*/[0]/* A */*1.0e-10/* m/A */)
//				   *f_i_c[2]/* N */;

//			N_clust[nNumClut][2] /* N*m */ 
//			+= 	   -(q/*_c[i_c]*/[1]/* A */*1.0e-10/* m/A */)
//					*f_i_c[0]/* N */
//			       +(q/*_c[i_c]*/[0]/* A */*1.0e-10/* m/A */)
//					*f_i_c[1]/* N */;
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//		}
//	}
	// LD、BD のとき揺動の計算を行う
	if (DYNMODE == LD || DYNMODE == BD)
	{
		Calc_Brownian();
	}
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
}

