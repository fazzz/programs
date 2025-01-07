#include <stdio.h>
#include <math.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "MD.h"
#include "force.h"
#include "physics.h"
#include "math.h"

double FPIeata_o;
double FPIeata_o2;

double Calc_stat_elect_pote_1_4(int nNumClut, int nNumClut2, int i_c, int ii_c);

double Calc_L_J_P_1_4(int nNumClut, int nNumClut2, int i_c,int ii_c);

// VDW 相互作用の計算を行う関数
double Calc_L_J_PotentialandForce(void)
{
	int i,j,ii;

	int alpha;

	int nNumClut, nNumClut2;

	int num_non_bond_clust;
	double eata;
	double eata_L_J[MAXA];
	double gamma_L_J[MAXA];
	double P_L_J[MAXA];
	double dihed_ang[10];
	double P_e_s[MAXA];
	double PI2;

	num_a_prot = 0;
	NUM_A_PROT = 0;

	FPIeata_o = 4.0*PI*eata_o/1.602177e-19/1.602177e-19*1.660539e-27*1.0e-24/1.0e-30;/*F^-1*m=m^3*kg*s^-4*A^-2*/
//	FPIeata_o2 = 4.0*PI*eata_o;/*F^-1*m=m^3*kg*s^-4*A^-2*/
	PI2 =  3.1415926536;
	eata = 1.0/(2.99792458e8*2.99792458e8*4*PI*1.0e-7);/*F^-1*m=m^3*kg*s^-4*A^-2*/
	FPIeata_o2 = 4.0*PI*eata_o;/*F^-1*m=m^3*kg*s^-4*A^-2*/

	// 静電相互作用ポテンシャルエネルギーの初期化
	for (i=0;i<prot.num_atom;++i)
	{
		potential_pro.p_elesta[i] = 0.0;
	}

	potential_pro.p_elestat=0.0;

	// 静電相互作用ポテンシャルエネルギーの初期化
	for (i=0;i<prot.num_atom;++i)
	{
		potential_pro.p_1_4_elesta[i] = 0.0;
	}

	for(i=0; i<prot.num_atom; ++i)
	{
		P_e_s[i] = 0.0;
	}


	potential_pro.p_L_Jt=0.0;
	for(nNumClut = 0;nNumClut < prot.DOF; ++nNumClut)
	{
		for(i=0;i < clust[nNumClut].num_atom_clust;++i)
		{
			eata_L_J[num_a_prot] = clust[nNumClut].f_p_clust.eata_L_J_p[i];
			gamma_L_J[num_a_prot] = clust[nNumClut].f_p_clust.gamma_L_J_p[i];
			for(alpha=0; alpha<3; ++alpha)
			{
				clust[nNumClut].f_c.f_L_J[i][alpha]/*N*/ = 0.0;
			}
			++num_a_prot;
		}
	}

	// VDW 相互作用ポテンシャルエネルギーの初期化
	for(i=0; i<prot.num_atom; ++i)
	{
		potential_pro.p_L_J[i] = 0.0;
	}

	// VDW 相互作用ポテンシャルエネルギーの初期化
	for(i=0; i<prot.num_atom; ++i)
	{
		potential_pro.p_1_4_L_J[i] = 0.0;
	}

	num_a_prot = 0;
	for(nNumClut = 0;nNumClut < prot.DOF; ++nNumClut)
	{
		for(i=0;i < clust[nNumClut].num_atom_clust;++i)
		{
			NUM_A_PROT = 0;
			for(nNumClut2 = 0;nNumClut2 < prot.DOF; ++nNumClut2)
			{
				for(j=0;j < clust[nNumClut2].num_atom_clust;++j)
				{
//					if (num_a_prot < NUM_A_PROT)
//					{
					if (which_calc_non_bonding_pote(nNumClut, i, NUM_A_PROT)
					                                                    == 1)
					{
///////////////////////////////////////////////////////////////////////////////
						potential_pro.p_elesta[num_a_prot]/*J*/ 
							  // 静電相互作用ポテンシャルエネルギーの計算
							  += Calc_stat_elect_pote(nNumClut, nNumClut2, i, /*ii*/j)
												 /**1.0e-3*/;
////					potential_pro.p_elestat += Calc_stat_elect_pote(nNumClut, nNumClut2, i, ii);
							// 静電相互作用力の計算を行う
							Calc_stat_elect_f(nNumClut, nNumClut2, i, /*ii*/j);
								/*N=J/m*/
////							Calc_stat_elect_f(nNumClut2, nNumClut, ii, i);
												 /*N=J/m*/
///////////////////////////////////////////////////////////////////////////////

//						printf("%d v.s. %d \n",num_a_prot+1, NUM_A_PROT+1);
						// VDW 相互作用ポテンシャルエネルギーの計算
						potential_pro.p_L_J[num_a_prot]/*J*/
						 += Calc_L_J_P(nNumClut, nNumClut2, i, j)/**1.0e-3*/;
						// VDW 相互作用力の計算
						Calc_L_J_F(nNumClut,i,j);
//						Calc_L_J_F(nNumClut,nNumClut2,i,j);
					}
					else if(which_calc_1_4_non_bonding_pote(nNumClut,
					                                         i, NUM_A_PROT)==1)
					{
///////////////////////////////////////////////////////////////////////////////
							potential_pro.p_1_4_elesta[num_a_prot]/*J*/ 
							  // 静電相互作用ポテンシャルエネルギーの計算
							  += /*1.0/1.2**/Calc_stat_elect_pote_1_4(nNumClut, nNumClut2, i, /*ii*/j)
												 /**1.0e-3*/;
							// 静電相互作用力の計算を行う
							Calc_stat_elect_f_1_4(nNumClut, nNumClut2, i, /*ii*/j);
												 /*N=J/m*/
////						Calc_stat_elect_f_1_4(nNumClut2, nNumClut, ii, i);
												 /*N=J/m*/
///////////////////////////////////////////////////////////////////////////////

						// VDW 相互作用ポテンシャルエネルギーの計算
						potential_pro.p_1_4_L_J[num_a_prot]/*J*/
						 += /*1.0/2.0**/Calc_L_J_P_1_4(nNumClut, nNumClut2, i, j)/**1.0e-3*/;
						// VDW 相互作用力の計算
						Calc_L_J_F_1_4(nNumClut,i,j);
//						Calc_L_J_F_1_4(nNumClut, nNumClut2, i,j);
					}
//					}
					++NUM_A_PROT;
				}
			}
			++num_a_prot;
		}
	}

///////////////////////////////////////////////////////////////////////////////

/*	for (i=0;i<prot.num_atom;++i)
	{
		potential_pro.p_elestat += potential_pro.p_elesta[i]/2.0;
	}

	for (i=0;i<prot.num_atom;++i)
	{
		potential_pro.p_1_4_elestat += potential_pro.p_1_4_elesta[i]/2.0;
	}

///////////////////////////////////////////////////////////////////////////////

	for (i=0;i<prot.num_atom;++i)
	{
		potential_pro.p_L_Jt += potential_pro.p_L_J[i]/2.0;
	}

	for (i=0;i<prot.num_atom;++i)
	{
		potential_pro.p_1_4_L_Jt += potential_pro.p_1_4_L_J[i]/2.0;
	}*/
}
///////////////////////////////////////////////////////////////////////////////

// 静電相互作用ポテンシャルエネルギーの計算を行う関数
double Calc_stat_elect_pote(int nNumClut, int nNumClut2, int i_c, int ii_c)
{
	FILE *out;
	double elect_pote = 0.0;

	if(len_q[num_a_prot][NUM_A_PROT]>0)
	{
		elect_pote /*J*/ = clust[nNumClut].f_p_clust.e_f[i_c]/*e*/
		                  *1.602177e-19/18.2223/*C=s*A*/
				 	      *clust[nNumClut2].f_p_clust.e_f[ii_c]/*e*/
				          *1.602177e-19/18.2223/*C=s*A*/
				          /FPIeata_o2/*m^3*kg*s^-4*A^-2*/
					      /(len_q[num_a_prot][NUM_A_PROT]/*A*/*1.0e-10/*m/A*/);
	}
	else
	{
		elect_pote = 0.0; /*J*/
		printf("error: length of %d th atom and %d th atom is not reasnoable\n", 
									                      num_a_prot, NUM_A_PROT);
	}

	elect_pote = elect_pote/**2.3889e-4/*Kcal/J*///
						   /4.184000/1000.0
						   *6.022142e23/*/mol*/
//						   *1.11268e-13
						   ;
//	out = fopen("testtest.out","a");
//	fprintf(out, "%d %d %lf %lf \n", num_a_prot+1, NUM_A_PROT+1, elect_pote, 1.0/(len_q[num_a_prot][NUM_A_PROT]*len_q[num_a_prot][NUM_A_PROT]));
//	fclose(out);

	return elect_pote; /*J*/

}

// 静電相互作用ポテンシャルエネルギーの計算を行う関数
double Calc_stat_elect_pote_1_4(int nNumClut, int nNumClut2, int i_c, int ii_c)
{
	double elect_pote = 0.0;

	if(len_q[num_a_prot][NUM_A_PROT]>0)
	{
		elect_pote /*J*/ = clust[nNumClut].f_p_clust.e_f[i_c]/*e*/
		                  *1.602177e-19/18.2223/*C=s*A*/
				 	      *clust[nNumClut2].f_p_clust.e_f[ii_c]/*e*/
				          *1.602177e-19/18.2223/*C=s*A*/
				          /FPIeata_o2/*m^3*kg*s^-4*A^-2*/
					      /(len_q[num_a_prot][NUM_A_PROT]/*A*/*1.0e-10/*m/A*/);
	}
	else
	{
		elect_pote = 0.0; /*J*/
		printf("error: length of %d th atom and %d th atom is not reasnoable\n", 
									                      num_a_prot, NUM_A_PROT);
	}

	elect_pote = elect_pote/*2.3889e-4*//*Kcal/J*/
							/4.184000/1000.0
						   *6.022142e23/1.2/*/mol*/;

	return elect_pote; /*J*/

}

// 静電相互作用力の計算を行う関数
void Calc_stat_elect_f(int nNumClut, int nNumClut2, int i_c, int ii_c)
{
	int alpha;

	double Th_len_q;

	FILE *out;

	Th_len_q = pow/*er*/(len_q[num_a_prot][NUM_A_PROT]/**1.0e-10*/,3);

//	out = fopen("FelestatA.out", "a");
//	fprintf(out,"%e \n",( clust[nNumClut].f_p_clust.e_f[i_c]/*e*/*1.602e-19
//				      *clust[nNumClut2].f_p_clust.e_f[ii_c]/*e*/*1.602e-19 )
//				     /FPIeata_o/*m^3*kg*s^-4*A^-2*/
//				     /Th_len_q/*A^3*/);
//	fclose(out);

	for(alpha=0; alpha<3; ++alpha)
	{
		clust[nNumClut].f_c.f_elesta[i_c][alpha]/*N*/
		/*-*/-= /*-(7_9)*/( clust[nNumClut].f_p_clust.e_f[i_c]/*e*//**1.602e-19*//18.2223
				      *clust[nNumClut2].f_p_clust.e_f[ii_c]/*e*//**1.602e-19*//18.2223 )
				     /FPIeata_o/*m^3*kg*s^-4*A^-2*/
				     /Th_len_q/*A^3*/
				     *(q[num_a_prot][NUM_A_PROT][alpha]/*A*//**1.0e-10*//*m/A*/)/**10e-2*/;
	}
}

// 静電相互作用力の計算を行う関数
void Calc_stat_elect_f_1_4(int nNumClut, int nNumClut2, int i_c, int ii_c)
{
	int alpha;

	double len;
	double Th_len_q;

	len = len_q[num_a_prot][NUM_A_PROT]/**1.0e-10*/;
	Th_len_q = pow/*er*/(len,3);

	for(alpha=0; alpha<3; ++alpha)
	{
		clust[nNumClut].f_c.f_1_4_elesta[i_c][alpha]/*N*/
		/*-*/-= /*-(7_9)*/( clust[nNumClut].f_p_clust.e_f[i_c]/18.2223/*e*//**1.602e-19*/
				           *clust[nNumClut2].f_p_clust.e_f[ii_c]/18.2223/*e*//**1.602e-19*/ )
						  /FPIeata_o/*m^3*kg*s^-4*A^-2*/
				          /Th_len_q/*m^3*/
				          *(q[num_a_prot][NUM_A_PROT][alpha]/*A*//**1.0e-10*//*m/A*/)/**10e-2*/
					      /1.2;
	}
}

//// 静電相互作用力の計算を行う関数
//void Calc_stat_elect_f(int nNumClut, int nNumClut2, int i_c, int ii_c)
//{
//	int alpha;
//
//	double Th_len_q;
//
//	for(alpha=0; alpha<3; ++alpha)
//	{
//		Th_len_q = power(len_q[num_a_prot][NUM_A_PROT]*1.0e-10,3);
//		clust[nNumClut].f_c.f_elesta[i_c][alpha]/*N*/
//				 += -( clust[nNumClut].f_p_clust.e_f[i_c]/*e*/*1.602e-19
//				      *clust[nNumClut2].f_p_clust.e_f[ii_c]/*e*/*1.602e-19 )
//				     /FPIeata_o/*m^3*kg*s^-4*A^-2*/
//				     /Th_len_q/*A^3*/
//				     *(q[num_a_prot][NUM_A_PROT][alpha]/*A*/*1.0e-10/*m/A*/);
//	}
//}

//// 静電相互作用力の計算を行う関数
//void Calc_stat_elect_f_1_4(int nNumClut, int nNumClut2, int i_c, int ii_c)
//{
//	int alpha;
//
//	double Th_len_q;
//
//	for(alpha=0; alpha<3; ++alpha)
//	{
//		Th_len_q = power(len_q[num_a_prot][NUM_A_PROT]*1.0e-10,3);
//		clust[nNumClut].f_c.f_1_4_elesta[i_c][alpha]/*N*/
//				 	+= -( clust[nNumClut].f_p_clust.e_f[i_c]/*e*/*1.602e-19
//				       *clust[nNumClut2].f_p_clust.e_f[ii_c]/*e*/*1.602e-19 )
//				       /FPIeata_o/*m^3*kg*s^-4*A^-2*/
//				       /Th_len_q/*A^3*/
//				       *(q[num_a_prot][NUM_A_PROT][alpha]/*A*/*1.0e-10/*m/A*/)
//					   /1.2;
//	}
//}

///////////////////////////////////////////////////////////////////////////////

// VDW 相互作用ポテンシャルエネルギーの計算を行う関数
double Calc_L_J_P(int nNumClut, int nNumClut2, int i_c,int ii_c)
{
	double P_L_J;
	double P_L_J_12;
	double P_L_J_6;
	double gammma_len;
	double one_len = 0.0;
	double A_i_j;
	double B_i_j;
	double l;

	int i,j,index;

	FILE *out;

//	eata_i_j/*J*/ =    ( clust[nNumClut].f_p_clust.eata_L_J_p[i_c]/*kcal/mol*/
//						   /6.022142e23/*mol*/
//						   *4.18e3/*J/kcal*/ 
//						)
//				     * ( clust[nNumClut2].f_p_clust.eata_L_J_p[ii_c]/*kcal/mol*/
//						   /6.022142e23/*mol*/
//						   *4.18e3/*J/kcal*/ 
//				 		);

//	if (eata_i_j > 0)
//	{
//		eata_i_j/*J*/ = sqrt(eata_i_j/*J^2*/);
//	}
//	else
//	{
//		eata_i_j/*J*/ = 0.0;
//		printf("error: parameter eata is not reasnable\n");
//	}

//	ganmma_i_j/*A*/
//	 = /*0.5*/(   clust[nNumClut].f_p_clust.gamma_L_J_p[i_c]
//				+ clust[nNumClut2].f_p_clust.gamma_L_J_p[ii_c]  )/*A*/;

//	A_i_j = eata_i_j*power(ganmma_i_j,12);
//	B_i_j = eata_i_j*power(ganmma_i_j,6);

	i = prot.L_J_parm.atomtype[num_a_prot]-1;
	j = prot.L_J_parm.atomtype[NUM_A_PROT]-1;

	index = prot.L_J_parm.atomtypeIndex[i][j]-1;

	A_i_j = prot.L_J_parm.A[index];
	B_i_j = prot.L_J_parm.B[index];

//	A_i_j = eata_i_j*power(ganmma_i_j,12);
//	B_i_j = eata_i_j*power(ganmma_i_j,6)*2.0;

//	A_i_j = A_i_j/6.022142e23/*mol*/*4.18e3/*J/kcal*/;
//	B_i_j = B_i_j/6.022142e23/*mol*/*4.18e3/*J/kcal*/;

//	if (len_q[num_a_prot][NUM_A_PROT] > 0)
//	{
//		gammma_len/*A*/ = ganmma_i_j/*A*/
//		                 /len_q[num_a_prot][NUM_A_PROT]/*A*/;
//	}
//	else
//	{
//		gammma_len = 0.0;
//		printf("error2: length of %d th atom and %d th atom is not reasnable\n"
//		                                             , num_a_prot, NUM_A_PROT);
//	}

	if (len_q[num_a_prot][NUM_A_PROT] <= 0)
	{
//		printf("error2: length of %d th atom and %d th atom is not reasnable\n"
//		                                             , num_a_prot, NUM_A_PROT);
	}

//	if (len_q[num_a_prot][NUM_A_PROT] > 3*ganmma_i_j)
//	{
//		P_L_J = 0.0;
//	}
//	else
//	{
//		P_L_J_12 = power(gammma_len,12);
//		P_L_J_6 = power(gammma_len,6);

//		P_L_J /*J*/ = eata_i_j/*J*/*(P_L_J_12-P_L_J_6); 
//	}

		one_len = 1.0/len_q[num_a_prot][NUM_A_PROT];

		P_L_J_12 = A_i_j*pow/*er*/(one_len,12);
		P_L_J_6  = B_i_j*pow/*er*/(one_len,6);

//		l = power(one_len,6);
		P_L_J /*J*/ = 
					P_L_J_12
					- P_L_J_6
					; 

//	out = fopen("testLJ","a");
//	fprintf(out,"%d %d %lf %lf %lf %lf\n", num_a_prot, NUM_A_PROT, P_L_J, power(one_len,6)*1.0e5, A_i_j, B_i_j);
//	fclose(out);
	return P_L_J; /*J*/
}

// VDW 相互作用ポテンシャルエネルギーの計算を行う関数
double Calc_L_J_P_1_4(int nNumClut, int nNumClut2, int i_c,int ii_c)
{
	double P_L_J;
	double P_L_J_12;
	double P_L_J_6;
	double gammma_len;
	double one_len = 0.0;
	double A_i_j;
	double B_i_j;
	double l;

	int i,j,index;

	FILE *out;

//	eata_i_j/*J*/ =    ( clust[nNumClut].f_p_clust.eata_L_J_p[i_c]/*kcal/mol*/
//						   /6.022142e23/*mol*/
//						   *4.18e3/*J/kcal*/ 
//						)
//				     * ( clust[nNumClut2].f_p_clust.eata_L_J_p[ii_c]/*kcal/mol*/
//						   /6.022142e23/*mol*/
//						   *4.18e3/*J/kcal*/ 
//				 		);

//	if (eata_i_j > 0)
//	{
//		eata_i_j/*J*/ = sqrt(eata_i_j/*J^2*/);
//	}
//	else
//	{
//		eata_i_j/*J*/ = 0.0;
//		printf("error: parameter eata is not reasnable\n");
//	}

//	ganmma_i_j/*A*/
//	 = /*0.5*/(   clust[nNumClut].f_p_clust.gamma_L_J_p[i_c]
//				+ clust[nNumClut2].f_p_clust.gamma_L_J_p[ii_c]  )/*A*/;

//	A_i_j = eata_i_j*power(ganmma_i_j,12);
//	B_i_j = eata_i_j*power(ganmma_i_j,6);

	i = prot.L_J_parm.atomtype[num_a_prot]-1;
	j = prot.L_J_parm.atomtype[NUM_A_PROT]-1;

	index = prot.L_J_parm.atomtypeIndex[i][j]-1;

	A_i_j = prot.L_J_parm.A[index];
	B_i_j = prot.L_J_parm.B[index];

//	A_i_j = eata_i_j*power(ganmma_i_j,12);
//	B_i_j = eata_i_j*power(ganmma_i_j,6)*2.0;

//	A_i_j = A_i_j/6.022142e23/*mol*/*4.18e3/*J/kcal*/;
//	B_i_j = B_i_j/6.022142e23/*mol*/*4.18e3/*J/kcal*/;

//	if (len_q[num_a_prot][NUM_A_PROT] > 0)
//	{
//		gammma_len/*A*/ = ganmma_i_j/*A*/
//		                 /len_q[num_a_prot][NUM_A_PROT]/*A*/;
//	}
//	else
//	{
//		gammma_len = 0.0;
//		printf("error2: length of %d th atom and %d th atom is not reasnable\n"
//		                                             , num_a_prot, NUM_A_PROT);
//	}

	if (len_q[num_a_prot][NUM_A_PROT] <= 0)
	{
//		printf("error2: length of %d th atom and %d th atom is not reasnable\n"
//		                                             , num_a_prot, NUM_A_PROT);
	}

//	if (len_q[num_a_prot][NUM_A_PROT] > 3*ganmma_i_j)
//	{
//		P_L_J = 0.0;
//	}
//	else
//	{
//		P_L_J_12 = power(gammma_len,12);
//		P_L_J_6 = power(gammma_len,6);

//		P_L_J /*J*/ = eata_i_j/*J*/*(P_L_J_12-P_L_J_6); 
//	}

		one_len = 1.0/len_q[num_a_prot][NUM_A_PROT];

		P_L_J_12 = A_i_j*pow/*er*/(one_len,12);
		P_L_J_6  = B_i_j*pow/*er*/(one_len,6);

		l = pow/*er*/(one_len,6);
		P_L_J /*J*/ = (P_L_J_12 - P_L_J_6)/2.0; 

//	out = fopen("testLJ","a");
//	fprintf(out,"%d %d %lf %lf %lf %lf\n", num_a_prot, NUM_A_PROT, P_L_J, power(one_len,6)*1.0e5, A_i_j, B_i_j);
//	fclose(out);
	return P_L_J; /*J*/
}

// VDW 相互作用力の計算を行う関数
void Calc_L_J_F(int nNumClut, int i_c, int j_c)
//void Calc_L_J_F(int nNumClut,int nNumClut2, int i_c, int j_c)
{
	int alpha,ii,jj,index;

	double F_L_J_A = 0.0;
	double F_L_J_A_12 = 0.0;
	double F_L_J_A_6 = 0.0;
	double F_L_J_k[3];
	double len;
	double A_i_j,B_i_j;

	FILE *out;

	for(alpha=0; alpha<3; ++alpha)
	{
		F_L_J_k[alpha] = 0.0;
	}

	ii = prot.L_J_parm.atomtype[num_a_prot]-1;
	jj = prot.L_J_parm.atomtype[NUM_A_PROT]-1;

	index = prot.L_J_parm.atomtypeIndex[ii][jj]-1;

//	A_i_j = prot.L_J_parm.A[index]/6.022142e23/*mol*/*4.18e3/*J/kcal*/*1.0e-150;
//	B_i_j = prot.L_J_parm.B[index]/6.022142e23/*mol*/*4.18e3/*J/kcal*/*1.0e-60;

	A_i_j = prot.L_J_parm.A[index]/6.022142e23/*mol*////2.3889e-4/*Kcal/J*//**1.0e-150;*/
					*4.184000*1000.0
				   /1.660539e-27/1.0e-20*1.0e-24/**1.0e120*/;
	B_i_j = prot.L_J_parm.B[index]/6.022142e23/*mol*////2.3889e-4/*Kcal/J*//**1.0e-60;*/
					*4.184000*1000.0
				   /1.660539e-27/1.0e-20*1.0e-24/**1.0e60*/;

	len = len_q[num_a_prot][NUM_A_PROT]/*A*//**1.0e-10*/;

	F_L_J_A_12/*m^-2*/ = /*2.0**//*24.0*/12.0*A_i_j/pow/*er*/(len,14);

	F_L_J_A_6/*m^-2*/ = 6.0*B_i_j/pow/*er*/(len,8);

	F_L_J_A/*J*m^-2*/ = /*24.0**/ /*(*/     
											F_L_J_A_12 
									 /*-*/- F_L_J_A_6/*)*/ /*m^-2*/
										;

//	out = fopen("FLJA.out", "a");
//	fprintf(out,"%e \n",F_L_J_A);
//	fclose(out);

	for(alpha=0;alpha<3;++alpha)
	{
		F_L_J_k[alpha]/*N=J/m*/
		             = F_L_J_A/*J*m^-2*/
					   *q[num_a_prot][NUM_A_PROT][alpha]/*A*/
					   /**1.0e-10*//*m/A*/;

		clust[nNumClut].f_c.f_L_J[i_c][alpha]/*N=J/m*/
		            /*+*/-= F_L_J_k[alpha]/*N=J/m*/;
	}
}

// VDW 相互作用力の計算を行う関数
void Calc_L_J_F_1_4(int nNumClut, int i_c, int j_c)
//void Calc_L_J_F_1_4(int nNumClut,int nNumClut2, int i_c, int j_c)
{
	int alpha,ii,jj,index;

	double F_L_J_A_1_4 = 0.0;
	double F_L_J_A_1_4_12 = 0.0;
	double F_L_J_A_1_4_6 = 0.0;
	double F_L_J_1_4_k[3];
	double len;
	double A_i_j,B_i_j;

	FILE *out;

	for(alpha=0; alpha<3; ++alpha)
	{
		F_L_J_1_4_k[alpha] = 0.0;
	}

	ii = prot.L_J_parm.atomtype[num_a_prot]-1;
	jj = prot.L_J_parm.atomtype[NUM_A_PROT]-1;

	index = prot.L_J_parm.atomtypeIndex[ii][jj]-1;

//	A_i_j = prot.L_J_parm.A[index]/6.022142e23/*mol*/*4.18e3/*J/kcal*/*1.0e-150;
//	B_i_j = prot.L_J_parm.B[index]/6.022142e23/*mol*/*4.18e3/*J/kcal*/*1.0e-60;

	A_i_j = prot.L_J_parm.A[index]/6.022142e23/*mol*////2.3889e-4/*Kcal/J*//**1.0e-150;*/
					*4.184000*1000.0
				   /1.660539e-27/1.0e-20*1.0e-24;
	B_i_j = prot.L_J_parm.B[index]/6.022142e23/*mol*////2.3889e-4/*Kcal/J*//**1.0e-60;*/
					*4.184000*1000.0
				   /1.660539e-27/1.0e-20*1.0e-24;

	len = len_q[num_a_prot][NUM_A_PROT]/*A*//**1.0e-10*/;

	F_L_J_A_1_4_12/*m^-2*/ = /*2.0**//*24.0*/12.0*A_i_j/pow/*er*/(len,/*13*/14);

	F_L_J_A_1_4_6/*m^-2*/ = 6.0*B_i_j/pow/*er*/(len,/*7*/8);

	F_L_J_A_1_4/*J*m^-2*/ = /*24.0**/ /*(*/ F_L_J_A_1_4_12 
									 /*-*/- F_L_J_A_1_4_6/*)*/ /*m^-2*/
										;

//	out = fopen("FLJA.out", "a");
//	fprintf(out,"%e \n",F_L_J_A);
//	fclose(out);

	for(alpha=0;alpha<3;++alpha)
	{
		F_L_J_1_4_k[alpha]/*N=J/m*/
		             = F_L_J_A_1_4/*J*m^-2*/
					   *q[num_a_prot][NUM_A_PROT][alpha]/*A*/
					   /**1.0e-10*//*m/A*/;

		clust[nNumClut].f_c.f_1_4_L_J[i_c][alpha]/*N=J/m*/
		            /*-*/-= F_L_J_1_4_k[alpha]/2.0/*N=J/m*/;
	}

//	int alpha,ii,jj,index;
//
//	double F_L_J_A;
//	double F_L_J_A_12, F_L_J_A_6;
//	double F_L_J_k[3];
//	double gammma_len;
//	double len;
//	double A_i_j,B_i_j;
//
//	for(alpha=0; alpha<3; ++alpha)
//	{
//		F_L_J_k[alpha] = 0.0;
//	}
//
//	ii = prot.L_J_parm.atomtype[num_a_prot]-1;
//	jj = prot.L_J_parm.atomtype[NUM_A_PROT]-1;
//
//	index = prot.L_J_parm.atomtypeIndex[ii][jj]-1;
//
////	A_i_j = prot.L_J_parm.A[index]/6.022142e23/*mol*/*4.18e3/*J/kcal*/*1.0e-150;
////	B_i_j = prot.L_J_parm.B[index]/6.022142e23/*mol*/*4.18e3/*J/kcal*/*1.0e-60;
//
//	A_i_j = prot.L_J_parm.A[index]/6.022142e23/*mol*//2.3889e-4/*Kcal/J*//**1.0e-150;*/
//				   /1.660539e-27/1.0e-20*1.0e-30;
//	B_i_j = prot.L_J_parm.B[index]/6.022142e23/*mol*//2.3889e-4/*Kcal/J*//**1.0e-60;*/
//				   /1.660539e-27/1.0e-20*1.0e-30;
//
////	gammma_len = ganmma_i_j/*A*//len_q[num_a_prot][NUM_A_PROT]/*A*/;
//
////	F_L_J_A_12/*m^-2*/ = 2*A_i_j
////	                    /power( len_q[num_a_prot][NUM_A_PROT]/*A*/
////	                           *1.0e-10,14);
//
////	F_L_J_A_6/*m^-2*/ = B_i_j
////	                    /power( len_q[num_a_prot][NUM_A_PROT]/*A*/
////	                           *1.0e-10,8);
//
//	len = len_q[num_a_prot][NUM_A_PROT]/*A*//**1.0e-10*/;
//
//	F_L_J_A_12/*m^-2*/ = /*2.0**/24.0*A_i_j/power(len,14);
//
//	F_L_J_A_6/*m^-2*/ = 6.0*B_i_j/power(len,8);
//
//	F_L_J_A/*J*m^-2*/ = /*24.0*(*/F_L_J_A_12
//								 -F_L_J_A_6/*)*//*m^-2*/
//								;
//
////	for(alpha=0;alpha<3;++alpha)
////	{
////		F_L_J_k[alpha]/*N=J/m*/
////		             = F_L_J_A/*J*m^-2*/
////					   *q[num_a_prot][NUM_A_PROT][alpha]/*A*/
////					   /**1.0e-10*//*m/A*/;
////
////		clust[nNumClut].f_c.f_1_4_L_J[i_c][alpha]/*N=J/m*/
////		            /*-*/+= 1.0/2.0*F_L_J_k[alpha]/*N=J/m*/;
////	}
}


