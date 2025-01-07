#include <stdio.h>
#include <math.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "MD.h"
#include "force.h"
#include "physics.h"

double FPIeata_o;

double Calc_stat_elect_pote_1_4(int nNumClut, int nNumClut2, int i_c, int ii_c);

// 静電相互作用の計算を行う関数
double Calc_ele_sta_PotentialandForce(void)
{
	int i,j,ii,alpha;

	int nNumClut, nNumClut2;

	FPIeata_o = 4*PI*eata_o;/*F^-1*m=m^3*kg*s^-4*A^-2*/

	int num_non_bond_clust;
	double P_e_s[MAXA];
	double dihed_ang[10];

	num_a_prot = 0;
	NUM_A_PROT = 0;

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

	num_a_prot = 0;
	for(nNumClut=0; nNumClut<prot.DOF; ++nNumClut)
	{
		for(i=0; i<clust[nNumClut].num_atom_clust; ++i)
		{
			NUM_A_PROT = 0;
			for(nNumClut2=0; nNumClut2<prot.DOF; ++nNumClut2)
			{
				for(ii=0; ii<clust[nNumClut2].num_atom_clust; ++ii)
				{
//					if (num_a_prot < NUM_A_PROT)
//					{
						if(which_calc_non_bonding_pote(nNumClut, i, NUM_A_PROT)==1)
						{
							potential_pro.p_elesta[num_a_prot]/*J*/ 
							  // 静電相互作用ポテンシャルエネルギーの計算
							  += Calc_stat_elect_pote(nNumClut, nNumClut2, i, ii)
												 /**1.0e-3*/;
//							potential_pro.p_elestat += Calc_stat_elect_pote(nNumClut, nNumClut2, i, ii);
							// 静電相互作用力の計算を行う
							Calc_stat_elect_f(nNumClut, nNumClut2, i, ii);
												 /*N=J/m*/
//							Calc_stat_elect_f(nNumClut2, nNumClut, ii, i);
												 /*N=J/m*/
						}
						else if(which_calc_1_4_non_bonding_pote(nNumClut,
						                                         i, NUM_A_PROT)==1)
						{
							potential_pro.p_1_4_elesta[num_a_prot]/*J*/ 
							  // 静電相互作用ポテンシャルエネルギーの計算
							  += 1.0/1.2*Calc_stat_elect_pote_1_4(nNumClut, nNumClut2, i, ii)
												 /**1.0e-3*/;
							// 静電相互作用力の計算を行う
							Calc_stat_elect_f_1_4(nNumClut, nNumClut2, i, ii);
												 /*N=J/m*/
//							Calc_stat_elect_f_1_4(nNumClut2, nNumClut, ii, i);
//												 /*N=J/m*/
						}
//					}
					++NUM_A_PROT;
				}
			}
			++num_a_prot;
		}
	}

	for (i=0;i<prot.num_atom;++i)
	{
		potential_pro.p_elestat += potential_pro.p_elesta[i]/2.0;
	}

	for (i=0;i<prot.num_atom;++i)
	{
		potential_pro.p_1_4_elestat += potential_pro.p_1_4_elesta[i]/2.0;
	}
//		potential_pro.p_elestat = potential_pro.p_elestat/2.0;
}

// 静電相互作用ポテンシャルエネルギーの計算を行う関数
//double Calc_stat_elect_pote(int nNumClut, int nNumClut2, int i_c, int ii_c)
//{
//	FILE *out;
//	double elect_pote = 0.0;
//
//	if(len_q[num_a_prot][NUM_A_PROT]>0)
//	{
//		elect_pote /*J*/ = clust[nNumClut].f_p_clust.e_f[i_c]/*e*/
//		                  *1.602e-19/*C=s*A*/
//				 	      *clust[nNumClut2].f_p_clust.e_f[ii_c]/*e*/
//				          *1.602e-19/*C=s*A*/
//				          /FPIeata_o/*m^3*kg*s^-4*A^-2*/
//					      /(len_q[num_a_prot][NUM_A_PROT]/*A*/*1.0e-10/*m/A*/);
//	}
//	else
//	{
//		elect_pote = 0.0; /*J*/
////		printf("error: length of %d th atom and %d th atom is not reasnoable\n", 
////									                      num_a_prot, NUM_A_PROT);
//	}
//
//	elect_pote = elect_pote*2.3889e-4/*Kcal/J*/
//						   *6.022142e23/*/mol*/;
////	out = fopen("testtest.out","a");
////	fprintf(out, "%d %d %lf %lf \n", num_a_prot+1, NUM_A_PROT+1, elect_pote, 1.0/(len_q[num_a_prot][NUM_A_PROT]*len_q[num_a_prot][NUM_A_PROT]));
////	fclose(out);
//
//	return elect_pote; /*J*/
//
//}
//
//// 静電相互作用ポテンシャルエネルギーの計算を行う関数
//double Calc_stat_elect_pote_1_4(int nNumClut, int nNumClut2, int i_c, int ii_c)
//{
//	double elect_pote = 0.0;
//
//	if(len_q[num_a_prot][NUM_A_PROT]>0)
//	{
//		elect_pote /*J*/ = clust[nNumClut].f_p_clust.e_f[i_c]/*e*/
//		                  *1.602e-19/*C=s*A*/
//				 	      *clust[nNumClut2].f_p_clust.e_f[ii_c]/*e*/
//				          *1.602e-19/*C=s*A*/
//				          /FPIeata_o/*m^3*kg*s^-4*A^-2*/
//					      /(len_q[num_a_prot][NUM_A_PROT]/*A*/*1.0e-10/*m/A*/);
//	}
//	else
//	{
//		elect_pote = 0.0; /*J*/
////		printf("error: length of %d th atom and %d th atom is not reasnoable\n", 
////									                      num_a_prot, NUM_A_PROT);
//	}
//
//	elect_pote = elect_pote*2.3889e-4/*Kcal/J*/
//						   *6.022142e23/*/mol*/;
//
//	return elect_pote; /*J*/
//
//}
//
//// 静電相互作用力の計算を行う関数
//void Calc_stat_elect_f(int nNumClut, int nNumClut2, int i_c, int ii_c)
//{
//	int alpha;
//
//	double Th_len_q;
//
//	FILE *out;
//
//	Th_len_q = power(len_q[num_a_prot][NUM_A_PROT]*1.0e-10,3);
//
//	out = fopen("FelestatA.out", "a");
//	fprintf(out,"%e \n",( clust[nNumClut].f_p_clust.e_f[i_c]/*e*/*1.602e-19
//				      *clust[nNumClut2].f_p_clust.e_f[ii_c]/*e*/*1.602e-19 )
//				     /FPIeata_o/*m^3*kg*s^-4*A^-2*/
//				     /Th_len_q/*A^3*/);
//	fclose(out);
//
//	for(alpha=0; alpha<3; ++alpha)
//	{
//		clust[nNumClut].f_c.f_elesta[i_c][alpha]/*N*/
//		+= /*-(7_9)*/( clust[nNumClut].f_p_clust.e_f[i_c]/*e*/*1.602e-19
//				      *clust[nNumClut2].f_p_clust.e_f[ii_c]/*e*/*1.602e-19 )
//				     /FPIeata_o/*m^3*kg*s^-4*A^-2*/
//				     /Th_len_q/*A^3*/
//				     *(q[num_a_prot][NUM_A_PROT][alpha]/*A*/*1.0e-10/*m/A*/)/**10e-1*/;
//	}
//}
//
//// 静電相互作用力の計算を行う関数
//void Calc_stat_elect_f_1_4(int nNumClut, int nNumClut2, int i_c, int ii_c)
//{
//	int alpha;
//
//	double Th_len_q;
//
//	Th_len_q = power(len_q[num_a_prot][NUM_A_PROT]*1.0e-10,3);
//
//	for(alpha=0; alpha<3; ++alpha)
//	{
//		clust[nNumClut].f_c.f_1_4_elesta[i_c][alpha]/*N*/
//		+= /*-(7_9)*/( clust[nNumClut].f_p_clust.e_f[i_c]/*e*/*1.602e-19
//				       *clust[nNumClut2].f_p_clust.e_f[ii_c]/*e*/*1.602e-19 )
//				       /FPIeata_o/*m^3*kg*s^-4*A^-2*/
//				       /Th_len_q/*A^3*/
//				       *(q[num_a_prot][NUM_A_PROT][alpha]/*A*/*1.0e-10/*m/A*/)/**10e-1*/
//					   /1.2;
//	}
//}

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
