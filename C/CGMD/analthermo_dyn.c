

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "physics.h"
#include "MD.h"
//#include "UmsSan.h"

// 現在の全エネルギー
double Energy_total;

// 現在の全運動エネルギー NVE
double Energy_kinetic;
double Energy_kinetic2;
double Energy_kinetic3;
double Energy_kinetic_o4;
double Energy_kinetic5;
double Energy_kinetic6;
double Energy_kinetic_o7[20];
//double Energy_kinetic8;
double Energy_kinetic8_dummy[MAXDOF][6];
double Energy_kinetic8_dummy2[MAXDOF];
double Energy_kinetictest;
double Energy_total_vertial;
double kinetic_energy_vertial;
double potential_energy_vertial;

double p_test;
// 現在の全運動エネルギー NVT
double Energy_kinetic_o;

double T_Kelvin_Now2;

double Energy_kinetic_tra;
double Energy_kinetic_rot;

double velo2[100][3];
double velo3[100][3];

void calc_velo(void);
//void calc_velo2(void);
void calc_velo_free_term(double vel_Term[3]);
void outprod(double v1[3],double v2[3], double v1x2[3]);

// 温度、運動エネルギー、
// ポテンシャルエネルギー、全エネルギーの計算とファイル出力を行う関数

double analthermo_dyn_properties(FILE *output, double vel_Term[3])
{
        int i,nNumClut;
	double kinetic_ene_t;
	FILE *outputfile;
	int degF;

	pre_dyn();
	/****************************/
        /* printf("a:55\n");	    */
        /****************************/
	CalcTotal_Energy(vel_Term);
	/****************************/
        /* printf("a:57\n");	    */
        /****************************/
	CalcT(vel_Term);
	/****************************/
        /* printf("a:59\n");	    */
        /****************************/

	Energy_total =  potential_pro.p_total + Energy_kinetic8;

	if (MODE == NVT)
	{
//		q_NVT = tau_NVT*deltat*tau_NVT*deltat*(prot.DOF-1)*k_B*T_Kelvin*2.3889e-4/*Kcal/J*/*6.022142e23/*/mol*/;
	  kinetic_energy_vertial = 0.5*q_NVT/4.18407/100.0
	    /**************************************************************************************************************/
            /* *1.660539e-27*1.0e-20/1.0e-24*2.3889e-4*6.022142e23/\**dot_s_NVT*dot_s_NVT/\*\//\**xi_NVT*xi_NVT*\/	  */
            /**************************************************************************************************************/
	    *(dot_s_NVT/s_NVT)*(dot_s_NVT/s_NVT);

	  if (TermMoveMode2==12)
	    degF=prot.DOF-1+6;
	  else
	    degF=prot.DOF-1;

          potential_energy_vertial = /*(prot.DOF-1)*//*deg*/DOFOFPROT*k_B_kcm/*1.38065e-23*/*T_Kelvin*log(s_NVT)
	  //	  potential_energy_vertial = Energy_kinetic8/4.18407/100.0*log(s_NVT)
	    /*********************************/
            /* *2.3889e-4*6.022142e23	     */
            /*********************************/
	    ;
	  //	  potential_energy_vertial = (prot.DOF-1/*+1*/)/**2.0*/*1.38065e-23*T_Kelvin*2.3889e-4/*Kcal/J*/*6.022142e23/*/mol*/*log(s_NVT);
//		s_NVT += xi_NVT*deltat*s_NVT*s_NVT;
//		s_NVT = s_NVT + xi_NVT*deltat*s_NVT;
//		dot_s_NVT = xi_NVT*s_NVT*s_NVT;
		Energy_total_vertial/*kcal/mol*/ =  potential_pro.p_total + Energy_kinetic8/*6*/
		                                   +potential_energy_vertial + kinetic_energy_vertial
		                                   ;
	}

	if (restflag==ON) {
	  potential_pro.p_restt = 0.0;
	  for (nNumClut=0;nNumClut<prot.DOF;++nNumClut) {
	    potential_pro.p_restt += potential_pro.p_rest[nNumClut];
	  }
	  Energy_total_vertial += potential_pro.p_restt;
	  Energy_total += potential_pro.p_restt;
	}


	if ((nNumStep % out_put_steps_thomo) == 0)
	{

	  if((outputfile = fopen(/*"thermo_dyn_properties.out"*/OutfilTHMO,"a")) == NULL) {
		  printf("error: can't open thermo_dyn_properties.out file");
		  exit(1);
	  	}


		// 温度、運動エネルギー、
		// ポテンシャルエネルギー、全エネルギーのファイル出力
		fprintf(outputfile,"/***********************************************/\n");
		fprintf(outputfile,"steps            = %d  \n",nNumStep);
		fprintf(outputfile,"total time       = %10.3lf ps  \n"
		                            ,deltat*(double)nNumStep/**1.0e12*/);
		fprintf(outputfile,"T_kelvin         = %e K  \n"
		                                      ,T_Kelvin_Now);
//		fprintf(output,"T_kelvin2         = %e K  \n"
//		                                      ,T_Kelvin_Now2);
		fprintf(outputfile,"toal_energy      = %e kcal/mol  \n"
		                                      ,Energy_total);
		if (MODE == NVT)
			fprintf(outputfile,"toal_vertial_energy      = %e kcal/mol  \n"
			                                      ,Energy_total_vertial);
		else 
			fprintf(outputfile,"toal_vertial_energy      = %e kcal/mol  \n"
			                                      ,0.0);
//		fprintf(output,"kinetic_energy1   = %e kcal/mol  \n"
//		                                    ,Energy_kinetic);
//		fprintf(output,"kinetic_energy2   = %e kcal/mol  \n"
//		                                    ,Energy_kinetic2);
//		fprintf(output,"kinetic_energy3   = %e kcal/mol  \n"
//		                                    ,Energy_kinetic3);
//		fprintf(output,"kinetic_energy4   = %e kcal/mol  \n"
//		                                    ,Energy_kinetic_o4);
//		fprintf(output,"kinetic_energy5   = %e kcal/mol  \n"
//		                                    ,Energy_kinetic5);
		fprintf(outputfile,"kinetic_energy6   = %e kcal/mol  \n"
		                                    ,Energy_kinetic6);
//		fprintf(output,"kinetic_energytest   = %e kcal/mol  \n"
//		                                    ,Energy_kinetictest);
//		fprintf(output,"p_test   = %e kcal/mol  \n"
//		                                    ,p_test);
		fprintf(outputfile,"kinetic_energy8   = %e kcal/mol  \n"
		                                    ,Energy_kinetic8);
		fprintf(outputfile,"kinetic_energy_vertial   = %e kcal/mol  \n"
		                                    ,kinetic_energy_vertial);
		fprintf(outputfile,"potential_energy_real = %e kcal/mol  \n"
		                                    ,potential_pro.p_total);
		fprintf(outputfile,"potential_energy_vertial   = %e kcal/mol  \n"
		                                    ,potential_energy_vertial);
		fprintf(outputfile,"s_NVT1   = %e   \n"
		                                    ,s_NVT);
		fprintf(outputfile,"s_NVT2   = %e   \n"
		                                    ,dot_s_NVT);
		fprintf(outputfile,"s_NVT3   = %e   \n"
		                                    ,acc_s_NVT);
		fprintf(outputfile,"xi   = %e   \n"
		                                    ,xi_NVT);
		fprintf(outputfile,"xi_dummy   = %e   \n"
		                                    ,xi_dummy);
///		for (i=0;i<prot.DOF;++i)
//		{
//			fprintf(output,"kinetic_energy7   = %d %e kcal/mol  \n"
//		                                    ,i,Energy_kinetic_o7[i]);
//		}
		fprintf(outputfile,"dihedral_energy  = %e kcal/mol  \n"
		                            ,potential_pro.p_dihedt);
		fprintf(outputfile,"elect_energy     = %e kcal/mol  \n"
		                           ,potential_pro.p_elestat);
		fprintf(outputfile,"VDW_energy       = %e kcal/mol  \n"
		                              ,potential_pro.p_L_Jt);
		fprintf(outputfile,"1_4_elect_energy = %e kcal/mol  \n"
		                           ,potential_pro.p_1_4_elestat);
		fprintf(outputfile,"1_4_VDW_energy   = %e kcal/mol  \n"
		                              ,potential_pro.p_1_4_L_Jt);

		if (restflag==ON) {
		  fprintf(outputfile,"dihedral_rest_energy  = %e kcal/mol  \n"
			  ,potential_pro.p_restt);
		}
		fprintf(outputfile,"\n");
 
		fclose(outputfile);
	}
	/*****************************/
        /* printf("a:178\n");	     */
        /*****************************/

//	potential_pro.p_total = 0.0;
//	potential_pro.p_dihedt = 0.0;
//	potential_pro.p_L_Jt = 0.0;
//	potential_pro.p_elestat = 0.0;
//	potential_pro.p_1_4_L_Jt = 0.0;
//	potential_pro.p_1_4_elestat = 0.0;

	return T_Kelvin_Now; // K

}

// 温度、運動エネルギー、
// ポテンシャルエネルギー、全エネルギーの計算を行う関数
void CalcTotal_Energy(double vel_Term[3])
{
	int i;

	int nNumClut;
	int nNumClut2;

	FILE *out;
	FILE *ouu;

	potential_pro.p_total = 0.0;
	potential_pro.p_dihedt = 0.0;
	potential_pro.p_L_Jt = 0.0;
	potential_pro.p_elestat = 0.0;
	potential_pro.p_1_4_L_Jt = 0.0;
	potential_pro.p_1_4_elestat = 0.0;

	double potential_atom[MAXA];
	double dot_s_NVT;

	Energy_total = 0.0;

//	if ((out=fopen("energy.out","a")) == NULL)
//	{
//		exit(1);
//	}

	for(i=0; i<prot.num_atom; ++i)
	{
		potential_atom[i]/*J*/ = 0.0;
	}

	for(i=0; i<prot.num_atom; ++i)
	{
		potential_atom[i]/*J*/    
		             =   potential_pro.p_L_J[i]/*J*/+ potential_pro.p_elesta[i]
		               + potential_pro.p_1_4_L_J[i] + potential_pro.p_1_4_elesta[i];
//		fprintf(out, "%d, %lf\n",i, potential_atom[i]*1.0e20);
	}
	/*****************************/
        /* printf("a:232\n");	     */
        /*****************************/
//	fprintf(out, "\n");
//	fclose(out);

//	if ((ouu=fopen("dihene.out","a")) == NULL)
//	{
//		exit(1);
//	}

//	for(i=0; i<prot.DOF-1; ++i)
//	{
//		fprintf(ouu, "%d, %lf\n",i, potential_pro.p_dihedc[i]*1.0e20);
//	}
//	fprintf(ouu, "\n");
//	fclose(ouu);

	// 非結合性相互作用相互作用ポテンシャルエネルギー
//	for(i=0; i<prot.num_atom; ++i)
//	{
//		potential_pro.p_L_Jt/*J*/    
//		             += potential_pro.p_L_J[i]/*J*/;
//		potential_pro.p_elestat/*J*/
//		          += potential_pro.p_elesta[i]/*J*/;
//	}

//	potential_pro.p_L_Jt/*J*/    
//		         = potential_pro.p_L_Jt*0.5/*J*/;
//	potential_pro.p_elestat/*J*/
//		         = potential_pro.p_elestat*0.5/*J*/;

//	for(i=0; i<prot.num_atom; ++i)
//	{
//		potential_pro.p_1_4_L_Jt/*J*/    
//		             += potential_pro.p_1_4_L_J[i]/*J*/;
//		potential_pro.p_1_4_elestat/*J*/
//		          += potential_pro.p_1_4_elesta[i]/*J*/;
//	}

//	potential_pro.p_1_4_L_Jt/*J*/    
//		         = potential_pro.p_1_4_L_Jt*0.5/*J*/;
//	potential_pro.p_1_4_elestat/*J*/
//		         = potential_pro.p_1_4_elestat*0.5/*J*/;

	// 2 面角相互作用ポテンシャルエネルギー
//	for(nNumClut=0; nNumClut<prot.DOF-1; ++nNumClut)
//	{
//		potential_pro.p_dihedt/*J*/
//		   += potential_pro.p_dihedc[nNumClut]/*J*/;
//	}

	// 温度、運動エネルギーの計算
	CalcT(vel_Term);
	/*****************************/
        /* printf("a:284\n");	     */
        /*****************************/

	// ポテンシャルエネルギー 単位系の変換 (J -> kcal/mol)
//	potential_pro.p_L_Jt    /*kcal/mol*/ = potential_pro.p_L_Jt/*J*/
//	                                       *2.3889e-4/*Kcal/J*/
//	                                       *6.022142e23/*/mol*/;

//	potential_pro.p_elestat /*kcal/mol*/ = potential_pro.p_elestat/*J*/
//	                                       *2.3889e-4/*Kcal/J*/
//						                   *6.022142e23/*/mol*/;

//	potential_pro.p_1_4_L_Jt /*kcal/mol*/ = potential_pro.p_1_4_L_Jt/*J*/
//	                                       *2.3889e-4/*Kcal/J*/
//	                                       *6.022142e23/*/mol*/;

//	potential_pro.p_1_4_elestat /*kcal/mol*/ = potential_pro.p_1_4_elestat/*J*/
//	                                       *2.3889e-4/*Kcal/J*/
//						                   *6.022142e23/*/mol*/;

//	potential_pro.p_dihedt  /*kcal/mol*/ = potential_pro.p_dihedt/*J*/
//	                                       *2.3889e-4/*Kcal/J*/
//						                   *6.022142e23/*//*mol*/;

//	p_test = 0.0;
	for (i=0;i<prot.num_atom;++i)
	{
		potential_pro.p_elestat += potential_pro.p_elesta[i]/2.0;
//		p_test += potential_pro.p_elesta[i]/2.0;
	}

	for (i=0;i<clust[0].num_atom_clust;++i)
	{
//		potential_pro.p_elestat -= potential_pro.p_elesta[i]/2.0;
//		p_test -= potential_pro.p_elesta[i]/2.0;
	}

	for (i=0;i<prot.num_atom;++i)
	{
		potential_pro.p_1_4_elestat += potential_pro.p_1_4_elesta[i]/2.0;
	}

	for (i=0;i<clust[0].num_atom_clust;++i)
	{
//		potential_pro.p_1_4_elestat -= potential_pro.p_1_4_elesta[i]/2.0;
	}

///////////////////////////////////////////////////////////////////////////////

	for (i=0;i<prot.num_atom;++i)
	{
		potential_pro.p_L_Jt += potential_pro.p_L_J[i]/2.0;
	}

	for (i=0;i<clust[0].num_atom_clust;++i)
	{
//		potential_pro.p_L_Jt -= potential_pro.p_L_J[i]/2.0;
	}

	for (i=0;i<prot.num_atom;++i)
	{
		potential_pro.p_1_4_L_Jt += potential_pro.p_1_4_L_J[i]/2.0;
	}

	for (i=0;i<clust[0].num_atom_clust;++i)
	{
//		potential_pro.p_L_Jt -= potential_pro.p_L_J[i]/2.0;
	}

	for (i=/*1*/0;i<prot.DOF;++i)
	{
		potential_pro.p_dihedt += potential_pro.p_dihedc[i];
	}

//	OldEnergy = potential_pro.p_total;

	// 全ポテンシャルエネルギー
	potential_pro.p_total/*kcal/mol*/
	  =   potential_pro.p_L_Jt    /*kcal/mol*/
	  + potential_pro.p_elestat /*kcal/mol*/
	  + potential_pro.p_dihedt  /*kcal/mol*/
	  + potential_pro.p_1_4_L_Jt    /*kcal/mol*/
	  + potential_pro.p_1_4_elestat /*kcal/mol*/
	  //	                                     + potential_US
	  ;
//	DeltaEnergy = potential_pro.p_total - OldEnergy;
	// 運動エネルギー 単位系の変換 (J -> kcal/mol)
//	Energy_kinetic/*kcal/mol*/ = Energy_kinetic/*J*/
//	                            *2.3889e-4/*Kcal/J*/
//	                            *6.022142e23/*/mol*/;

	// 全エネルギー
	Energy_total/*kcal/mol*/ =   potential_pro.p_total/*kcal/mol*/ + Energy_kinetic8/*kcal/mol*/;

//	if (MODE == NVT)
//	{
//		q_NVT = tau_NVT*deltat*tau_NVT*deltat*(prot.DOF-1)*k_B*T_Kelvin*2.3889e-4/*Kcal/J*/*6.022142e23/*/mol*/;
//		kinetic_energy_vertial = 0.5*q_NVT/**dot_s_NVT*dot_s_NVT/*//**xi_NVT*xi_NVT*/
//		                                        *(dot_s_NVT/s_NVT)*(dot_s_NVT/s_NVT);
//		potential_energy_vertial = (prot.DOF-1/*+1*/)/**2.0*/*k_B*T_Kelvin
//                                   *2.3889e-4/*Kcal/J*/
//	                               *6.022142e23/*/mol*/*log(s_NVT);
////		s_NVT += xi_NVT*deltat*s_NVT*s_NVT;
////		s_NVT = s_NVT + xi_NVT*deltat*s_NVT;
////		dot_s_NVT = xi_NVT*s_NVT*s_NVT;
//		Energy_total_vertial/*kcal/mol*/ =  potential_pro.p_total + Energy_kinetic8
//		                                   + 0.5*q_NVT/**dot_s_NVT*dot_s_NVT/*//**xi_NVT*xi_NVT*/
//		                                        *(dot_s_NVT/s_NVT)*(dot_s_NVT/s_NVT)
//		                                   /*+*/+ (prot.DOF-1/*+1*/)/**2.0*/*k_B*T_Kelvin
//                                            *2.3889e-4/*Kcal/J*/
//	                                         *6.022142e23/*/mol*/*log(s_NVT)
//		                                   ;
//	}
	/*****************************/
        /* printf("a:396\n");	     */
        /*****************************/

}

// 温度、運動エネルギーの計算を行う関数
void CalcT(double vel_Term[3])
{
	int nNumClut;

	int nNumClut2;

	int i,alpha,alpha2;

	int nNumAtomOrigClut;

	int nNumAtom;

	int nNumAtomClut;

	double sumMass;

	double Inertia = 0.0;

	double 	Energy_kinetic_clust[MAXDOF];

	int nNumA;

	int nNumBod;

	double dtheta;

	double velo[3];
	double angvelo[3];

	double kinetic_ene_t;

	FILE *out;
	FILE *out2;

	int DegOfFed;
	int degF;

	Energy_kinetic_o = 0.0;
	Energy_kinetic_o2 = 0.0;
	Energy_kinetic_o4 = 0.0;
	Energy_kinetic_tra = 0.0;
	Energy_kinetic_rot = 0.0;
	T_Kelvin_Now = 0.0;

//	for(nNumClut=0; nNumClut<prot.DOF/*-1*/; ++nNumClut)
//	{
//		nNumA = clust[nNumClut].terminal_atom_a[0] - 1;
//		for (i=0; i<clust[nNumClut].num_atom_clust; ++i)
//		{
//			for (alpha=0; alpha<3; ++alpha)
//			{
//				Energy_kinetic_o2
//				+= 0.50 * (prot.velo[nNumA+i][alpha])
//						* (prot.velo[nNumA+i][alpha])
//						* (clust[nNumClut].mass_clust[i]
//						*1.660539e-27);
//			}
//		}
//	}
//
////	for (alpha=0; alpha<3; ++alpha)
////	{
////		Energy_kinetic_o2
////			-= 0.50 * (prot.veloCOM[alpha])
////					* (prot.veloCOM[alpha])
////					* prot.sumMass
////					*1.660539e-27;
////	}
//
//	Energy_kinetic2 = Energy_kinetic_o2;
//	Energy_kinetic2 = Energy_kinetic2*2.3889e-4/*Kcal/J*/
//	                                 *6.022142e23/*/mol*/;
//
////	for(nNumClut=0; nNumClut<prot.DOF/*-1*/; ++nNumClut)
////	{
////		Energy_kinetic_o4
////				+= 0.50 * (clust[nNumClut].ddihedang[0])*1.0e15
////						* (clust[nNumClut].ddihedang[0])*1.0e15
////						* (clust[nNumClut].Inertia_clust_total)
////						*1.660539e-27*1.0e-20;
////	}
//
//	Energy_kinetic_o4 = Energy_kinetic_o4*2.3889e-4/*Kcal/J*/
//	                                     *6.022142e23/*/mol*/;
//
//	T_Kelvin_Now2 = Energy_kinetic_o2;
//
//	/////////////////////////////////////////////////////////////////////////////
//  /////////////////////////////////////////////////////////////////////////////
//
////	for(nNumClut=0; nNumClut<prot.DOF-1; ++nNumClut)
////	{
////		nNumA = clust[nNumClut].terminal_atom_a[0] - 1;
////		for (i=0; i<clust[nNumClut].num_atom_clust; ++i)
////		{
////			for (alpha=0; alpha<3; ++alpha)
////			{
////				Energy_kinetic_o
////				+= 0.50 * (prot.velo[nNumA+i][alpha]
////				           /**1.0e-10/1.0e-15*/)
////						* (prot.velo[nNumA+i][alpha]
////						   /**1.0e-10/1.0e-15*/)
////						* (clust[nNumClut].mass_clust[i]
////						   *1.660539e-27);
////			}
////		}
////	}
//
//	// 運動エネルギーを計算
//  /////////////////////////////////////////////////////////////////////////////
//  /////////////////////////////////////////////////////////////////////////////
////	for (nNumClut=0; nNumClut<prot.DOF-1; ++nNumClut)
////	{
////		nNumAtomOrigClut = clust[nNumClut].origin_atom_a;
////
//// 		sumMass = 0.0;
////		for (i=0;i<clust[nNumClut].num_atom_clust;++i)
////		{
////			sumMass += clust[nNumClut].mass_clust[i];
////		}
////
////		// 並進運動エネルギー
////		for (alpha=0; alpha<3; ++alpha)
////		{
////			Energy_kinetic_tra
////			 += 0.5*(prot.velo[nNumAtomOrigClut][alpha]*1.0e-10)
////			       *(prot.velo[nNumAtomOrigClut][alpha]*1.0e-10)
////			       *sumMass*1.660539e-27;
////		}
////
////		// 回転運動エネルギー
////		Energy_kinetic_rot
////		  += 0.5*clust[nNumClut].ddihedang[0]
////				*clust[nNumClut].ddihedang[0]
////			    *clust[nNumClut].Inertia_clust[2][2];
////	}
////
////	Energy_kinetic_o = Energy_kinetic_tra + Energy_kinetic_rot;
//
//	for (nNumClut=0; nNumClut<prot.DOF/*-1*/; ++nNumClut)
//	{
//		nNumAtomOrigClut = clust[nNumClut].origin_atom_a;
//
// 		sumMass = 0.0;
//		for (i=0;i<clust[nNumClut].num_atom_clust;++i)
//		{
//			sumMass += clust[nNumClut].mass_clust[i];
//		}
//
//		for (alpha=0; alpha<3; ++alpha)
//		{
//			velo[alpha] = 0.0;
//			angvelo[alpha] = 0.0;
//		}
//
//
//		for (alpha=0; alpha<3; ++alpha)
//		{
//			for (alpha2=0; alpha2<3; ++alpha2)
//			{
//				velo[alpha] += clust[nNumClut].trans_A_to_CN[0][alpha2][alpha]
//							  *clust[nNumClut].sp_velo[alpha2+3];
//				angvelo[alpha] += clust[nNumClut].trans_A_to_CN[0][alpha2][alpha]
//							  *clust[nNumClut].sp_velo[alpha2];
//			}
//		}
//		
//		// 並進運動エネルギー
//		for (alpha=0; alpha<3; ++alpha)
//		{
//			Energy_kinetic_tra
//			 += 0.5*(velo[alpha])
//			       *(velo[alpha])
//			       *sumMass*1.660539e-27;
//		}
//
//		// 回転運動エネルギー
//		for (alpha=0; alpha<3; ++alpha)
//		{
//			for (alpha2=0; alpha2<3; ++alpha2)
//			{
//				Energy_kinetic_rot
//				  += 0.5*angvelo[alpha]
//						*clust[nNumClut].Inertia_clust[alpha][alpha2]
//					    *angvelo[alpha2];
//			}
//		}
//	}
//
//	Energy_kinetic_o = Energy_kinetic_tra + Energy_kinetic_rot;
//
////	for (nNumClut=0; nNumClut<prot.DOF-1; ++nNumClut)
////	{
////		nNumAtomOrigClut = clust[nNumClut].origin_atom_a;
////
//// 		sumMass = 0.0;
////		for (i=0;i<clust[nNumClut].num_atom_clust;++i)
////		{
////			sumMass += clust[nNumClut].mass_clust[i];
////		}
////
////		// 並進運動エネルギー
////		for (alpha=4; alpha<6; ++alpha)
////		{
////			Energy_kinetic_tra
////			 += 0.5*(clust[nNumClut].sp_velo[alpha])
////			       *(clust[nNumClut].sp_velo[alpha])
////			       *sumMass*1.660539e-27;
////		}
////
////		// 回転運動エネルギー
////		Energy_kinetic_rot
////		  += 0.5*clust[nNumClut].ddihedang[0]
////				*clust[nNumClut].ddihedang[0]
////			    *clust[nNumClut].Inertia_clust[2][2];
////	}
////
////	Energy_kinetic_o = Energy_kinetic_tra + Energy_kinetic_rot;
//
//  /////////////////////////////////////////////////////////////////////////////
//  /////////////////////////////////////////////////////////////////////////////
//
//	// 運動エネルギーを計算
// /////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////////////////////////
///*	for (nNumClut=0; nNumClut<prot.DOF-1; ++nNumClut)
//	{
//		// 回転運動エネルギー
//		Energy_kinetic_o
//		  += 0.5
//		        *clust[nNumClut].sp_velo[0]
//				*clust[nNumClut].sp_velo[0]
//			    *clust[nNumClut].Inertia_clust[0][0]
//			   +
//			    clust[nNumClut].sp_velo[0]
//				*clust[nNumClut].sp_velo[1]
//			    *clust[nNumClut].Inertia_clust[0][1]
//			   +
//			    clust[nNumClut].sp_velo[0]
//				*clust[nNumClut].sp_velo[2]
//			    *clust[nNumClut].Inertia_clust[0][2]
//			   +
//			    clust[nNumClut].sp_velo[1]
//				*clust[nNumClut].sp_velo[1]
//			    *clust[nNumClut].Inertia_clust[1][1]
//			   +
//			    clust[nNumClut].sp_velo[1]
//				*clust[nNumClut].sp_velo[0]
//			    *clust[nNumClut].Inertia_clust[1][0]
//			   +
//			    clust[nNumClut].sp_velo[1]
//				*clust[nNumClut].sp_velo[2]
//			    *clust[nNumClut].Inertia_clust[1][2]
//			   +
//			    clust[nNumClut].sp_velo[2]
//				*clust[nNumClut].sp_velo[2]
//			    *clust[nNumClut].Inertia_clust[2][2]
//			   +
//			    clust[nNumClut].sp_velo[2]
//				*clust[nNumClut].sp_velo[0]
//			    *clust[nNumClut].Inertia_clust[2][0]
//			   +
//			    clust[nNumClut].sp_velo[2]
//				*clust[nNumClut].sp_velo[1]
//			    *clust[nNumClut].Inertia_clust[2][1];
//	}*/
//	///////////////////////////////////////////////////////////////////////
//  ///////////////////////////////////////////////////////////////////////
//
//	// 運動エネルギーを計算
//	///////////////////////////////////////////////////////////////////////
//  /////////////////////////////////////////////////////////////////////
///*	for (nNumClut=0; nNumClut<prot.DOF-1; ++nNumClut)
//	{
//		if (clust[nNumClut].terminal == TERMINAL)
//		{
//			// 回転運動エネルギー
//			Energy_kinetic_clust[nNumClut]
//		  		= 0.5*clust[nNumClut].ddihedang[0]
//				*clust[nNumClut].ddihedang[0]
//			    *clust[nNumClut].Inertia_clust[2][2];
//
//		}
//		else if (clust[nNumClut+1].terminal == TERMINAL)
//		{
//		for (nNumClut2=nNumClut+1; nNumClut2<prot.DOF-1; ++nNumClut2)
//		{
//			Inertia += clust[nNumClut2].Inertia_clust[2][2];
//		}
//
//		// 回転運動エネルギー
//		Energy_kinetic_clust[nNumClut]
//	  		= 0.5*clust[nNumClut].ddihedang[0]
//			*clust[nNumClut].ddihedang[0]
//		    *Inertia;
//
//		Inertia = 0.0;
//	}
//	else
//	{
//		for (nNumClut2=nNumClut+1; nNumClut2<prot.DOF-1; ++nNumClut2)
//		{
//			Inertia += clust[nNumClut2].Inertia_clust[2][2];
//		}
//
//		// 回転運動エネルギー
//		Energy_kinetic_clust[nNumClut]
//	  		= 0.5*clust[nNumClut].ddihedang[0]
//			*clust[nNumClut].ddihedang[0]
//		    *Inertia;
//
//		Inertia = 0.0;
//	}
//
//}
//
//Energy_kinetic_rot = 0.0;
//
//for (nNumClut=1;nNumClut<prot.DOF-1;++nNumClut)
//{
//	Energy_kinetic_rot += Energy_kinetic_clust[nNumClut];
//}
//
//Energy_kinetic_o = Energy_kinetic_rot;*/
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
//for (nNumClut=0; nNumClut<prot.DOF/*-1*/; ++nNumClut)
//{
//	for (nNumClut2=nNumClut; nNumClut2<prot.DOF/*-1*/; ++nNumClut2)
//	{
//		// 回転運動エネルギー
//		Energy_kinetic_clust[nNumClut]
//	  		= 0.5*clust[nNumClut].ddihedang[0]
//			*clust[nNumClut].ddihedang[0]
//		    *clust[nNumClut2].Inertia_clust[2][2];
//	}
//}
//
//Energy_kinetic_rot = 0.0;
//
//	for (nNumClut=0;nNumClut<prot.DOF/*-1*/;++nNumClut)
//	{
//		Energy_kinetic_rot += Energy_kinetic_clust[nNumClut];
//	}
//
//	Energy_kinetic3 = Energy_kinetic_rot*2.3889e-4/*Kcal/J*/
//	                                    *6.022142e23/*/mol*/;
//
//	T_Kelvin_Now/*K*/
//	      = Energy_kinetic6/*_o*//*J*/
//			/(k_B/*J/K*/
//			*(prot.num_atom*3-clust[0].num_atom_clust*3));
//
//	T_Kelvin_Now2/*K*/
//	      = Energy_kinetic_o2/*J*/
//			/(k_B/*J/K*/
//			*(prot.num_atom*3-clust[0].num_atom_clust*3));
//
//	// NVTの場合、速度スケーリングを行う
//	if ((MODE == NVT /*|| T_Kelvin_Now > T_Kelvin + 300.0*/) /*&& nNumStep > 2*/)
//	{
//		Energy_kinetic/*J*/
//		 = velocity_scaling(Energy_kinetic_o/*J*/);
//	}
//	else
//	{
//		Energy_kinetic/*J*/
//		  = Energy_kinetic_o/*J*/;
//	}
//
//	// 温度の計算を行う
////	T_Kelvin_Now/*K*/
////	      = Energy_kinetic/*J*/
////			/(k_B/*J/K*/
////			*(prot.num_atom*3-clust[0].num_atom_clust*3));
//
//	// 温度の計算を行う
////	T_Kelvin_Now/*K*/
////	      = Energy_kinetic/*J*/
//			/(k_B/*J/K*/
//			*(prot.DOF-1)*3);
//
//	// 温度の計算を行う
//	T_Kelvin_Now/*K*/
//	      = Energy_kinetic/*J*/
//			/(k_B/*J/K*/
//			*(prot.num_atom*3-clust[0].num_atom_clust*3));
//
//	T_Kelvin_Now = 0.5*Energy_kinetic6/(prot.DOF*k_B);
//

//	calc_velo();
	if ((nNumStep % out_put_steps)==0 && veloutflag == ON) {
          if (TermMoveMode2==12)
	    calc_velo_free_term(vel_Term);
	  else
	    calc_velo2();///*030410*/
	  /*****************************/
          /* printf("a:793\n");	       */
          /*****************************/
	  nNumAtom = 0;
	  Energy_kinetic6 = 0.0;
	  for (nNumClut=0;nNumClut<prot.DOF;++nNumClut) {
	    for (nNumAtomClut=0;nNumAtomClut<clust[nNumClut].num_atom_clust;++nNumAtomClut) {
	      for (alpha=0;alpha<3;++alpha) {
		Energy_kinetic6 += 0.5*clust[nNumClut].mass_clust[nNumAtomClut]*prot.velo[nNumAtom][alpha]*prot.velo[nNumAtom][alpha];
	      }
	      ++nNumAtom;
	    }
	  }
	  Energy_kinetictest = Energy_kinetic6;
	  Energy_kinetic6 = Energy_kinetic6/(4.18407*100.0);
	}
	/*****************************/
        /* printf("a:807\n");	     */
        /*****************************/

//
//	nNumAtom = 0;
//	Energy_kinetic5 = 0.0;
//	for (nNumClut=0;nNumClut<prot.DOF;++nNumClut)
//	{
////		for (nNumAtomClut=0;nNumAtomClut<clust[nNumClut].num_atom_clust;++nNumAtomClut)
////		for (nNumClutJ=0;nNumClutJ<prot.DOF;++nNumClutJ)
////		{
//			for (alpha=0;alpha<3;++alpha)
//			{
//				Energy_kinetic5 += 0.5*clust[nNumClut].sum_mass/*.mass_clust[nNumAtomClut]*/*1.660539e-27
//								      *velo2[/*nNumAtom*/nNumClut][alpha]*1.0e-10*1.0e12
//								      *velo2[/*nNumAtom*/nNumClut][alpha]*1.0e-10*1.0e12;
//			}
////			++nNumAtom;
////		}
//	}
//	Energy_kinetic5 = Energy_kinetic5*2.3889e-4/*Kcal/J*/
//	                                 *6.022142e23/*/mol*/;
//
////	if ((out = fopen("velo.txt","a")) == NULL)
////	{
////		exit(1);
////	}
//
////	for (nNumAtom=0;nNumAtom<prot.num_atom;++nNumAtom)
////	{
////		for (alpha=0;alpha<3;++alpha)
////		{
////			fprintf(out,"%lf ",velo2[nNumAtom][alpha]);
////		}
////		fprintf(out,"\n ");
////	}
////	fclose(out);
//


//	T_Kelvin_Now = 0.5*Energy_kinetic6/(prot.DOF*k_B);

//	if ((out2 = fopen("velo2.txt","a")) == NULL)
//	{
//		exit(1);
//	}
//
//	for (nNumAtom=0;nNumAtom<prot.num_atom;++nNumAtom)
//	{
//		for (alpha=0;alpha<3;++alpha)
//		{
//			fprintf(out2,"%lf ",velo3[nNumAtom][alpha]);
//		}
//		fprintf(out,"\n ");
//	}
//	fclose(out);

//	Energy_kinetic_o7 = 0.0;
	for(nNumClut=0; nNumClut<prot.DOF/*-1*/; ++nNumClut)
	{
		Energy_kinetic_o7[nNumClut] = 0.0;
		nNumA = clust[nNumClut].terminal_atom_a[0] - 1;
		for (i=0; i<clust[nNumClut].num_atom_clust; ++i)
		{
			for (alpha=0; alpha<3; ++alpha)
			{
			  //030410				Energy_kinetic_o7[nNumClut]
			  //+= 0.50 * (velo3[nNumA+i][alpha]*1.0e2)
			  //		* (velo3[nNumA+i][alpha]*1.0e2)
			  //		* (clust[nNumClut].mass_clust[i]
			  //		*1.660539e-27)
			  //		/4.184000/1000.0
			  //*6.022142e23/*/mol*/;
			}
		}
	}

	Energy_kinetic8 = 0.0;
	for(nNumClut=/*1*/0; nNumClut<prot.DOF/*-1*/; ++nNumClut){
	  for (alpha=0; alpha<6; ++alpha){
	    Energy_kinetic8_dummy[nNumClut][alpha] = 0.0;
	  }
	  Energy_kinetic8_dummy2[nNumClut] = 0.0;
	}

	for(nNumClut=/*1*/0; nNumClut<prot.DOF/*-1*/; ++nNumClut){
	  for (alpha=0; alpha<6; ++alpha){
	    for (alpha2=0; alpha2<6; ++alpha2){
	      Energy_kinetic8_dummy[nNumClut][alpha]
		+=       clust[nNumClut].sp_velo[alpha2]
		* clust[nNumClut].InertiaMatrix[alpha][alpha2];
	    }
	  }

	  for (alpha=0; alpha<6; ++alpha){
	    Energy_kinetic8_dummy2[nNumClut]
	      +=   Energy_kinetic8_dummy[nNumClut][alpha]
	      * clust[nNumClut].sp_velo[alpha];
	  }
	  Energy_kinetic8 += Energy_kinetic8_dummy2[nNumClut];
	}

	/***********************************************************/
        if (TermMoveMode2 == 3) {
	  Energy_kinetic8 = 0.0;
	  for (alpha=0; alpha<3; ++alpha)	{
	    Energy_kinetic8_dummy[0][alpha] = 0.0;
	  }
	  Energy_kinetic8_dummy2[0] = 0.0;
	  for (alpha=0; alpha<3; ++alpha) {
	    for (alpha2=0; alpha2<3; ++alpha2) {
	      Energy_kinetic8_dummy[0][alpha]
		+=       clust[0].sp_velo[alpha2]
		* clust[0].InertiaMatrix[alpha][alpha2];
	    }
	  }
	
	  for (alpha=0; alpha<3; ++alpha){
	    Energy_kinetic8_dummy2[0]
	    +=   Energy_kinetic8_dummy[0][alpha]
	      * clust[0].sp_velo[alpha];
	  }
	  Energy_kinetic8 += Energy_kinetic8_dummy2[0];
	/************/
        }
        /************/
        /***********************************************************/
	//	DegOfFed=(prot.num_atom-clust[0].num_atom_clust-1)*3;
 	if (TermMoveMode2==12)
 	  degF=prot.DOF-1+6;
 	else
 	  degF=prot.DOF-1;
	T_Kelvin_Now = /*2.0(3III10)**/Energy_kinetic8
	                                        /4.18407/100.0
	                                        /**************************/
                                                /* *1.660539e-27	  */
						/* *1.0e-10*1.0e12	  */
						/* *1.0e-10*1.0e12	  */
                                                /**************************/
	                                        // /((/*prot.DOF-1*/degF)*k_B_kcm);
	                                        /((DOFOFPROT)*k_B_kcm);
	Energy_kinetic8 = 0.5*Energy_kinetic8/(4.18407*100.0);
 	/*****************************/
        /* printf("a:922\n");	     */
        /*****************************/
	  //						*1.660539e-27
	  //					*1.0e-10*1.0e12
	  //					*1.0e-10*1.0e12
	  //					*2.3889e-4
//						/4.184000/1000.0
	  //					*6.022142e23/*/mol*/;

	//	kinetic_ene_t =     0.5*1.0*1.0*1.0*clust[1].ddihedang[0]*clust[1].ddihedang[0]*1.660539e-27*1.0e-20/1.0e-24*2.3889e-4*6.022142e23/*kcal/mol*/
	//	               +0.5*1.0*1.0*1.0*(clust[2].ddihedang[0]+clust[1].ddihedang[0])*(clust[2].ddihedang[0]+clust[1].ddihedang[0])*1.660539e-27*1.0e-20/1.0e-24*2.3889e-4*6.022142e23
	//                   +0.5*1.0*1.0*1.0*clust[1].ddihedang[0]*clust[1].ddihedang[0]*1.660539e-27*1.0e-20/1.0e-24*2.3889e-4*6.022142e23
	//                   -1.0*1.0*1.0*cos(clust[2].dihedang[0])*(clust[1].ddihedang[0]*clust[2].ddihedang[0]+clust[1].ddihedang[0]*clust[1].ddihedang[0])*1.660539e-27*1.0e-20/1.0e-24*2.3889e-4*6.022142e23/*kcal/mol*/
	//                    ;

////	Energy_kinetic_o7 = Energy_kinetic_o7*2.3889e-4/*Kcal/J*/
//	                                     *6.022142e23/*/mol*/;

}


// 現在の構造でのタンパク質のポテンシャルエネルギー、
// タンパク質に及ぼす力の計算を行う関数
void calc_velo(void)
{
	int i,ii,i_c,alpha,alpha2,alpha3,j,k;

	int nNumClut;
	int nNumClutI;
	int nNumClutJ;
	int nNumAtom;

	int nNumOrigc;
	int nNumAbsoc;

	int origin_of_this_branch = 0;

	int nNumClutOfParent;
	int nNmAtomOfCN_1_A;

	int num;

	double q_c[MAXNAC][3];

	double qq[3], qq2[3];

	double RotatnNumtonNumMiOn[3][3];

	double mat[4][4];
	double dummy[4][4];
	double dummy_xoord[4];
	double dummy_xoord2[4];
	double dotTransQOne[10][100][4];
	double dotTransQTwo[10][100][4];
	double dotTransQ[10][100][4];
	double dotTransQ_mat[10][40];


	for (nNumClut=1;nNumClut<=prot.DOF-1;++nNumClut)
	{
//		for (nNumClutJ=1;nNumClutJ<=prot.DOF-1;++nNumClutJ)
		for (nNumAtom=1;nNumAtom<=prot.num_atom;++nNumAtom)
		{
			for (i=0;i<4;++i)
			{
				dotTransQ[nNumClut][nNumAtom][i] = 0.0;
			}
		}
	}

	for (nNumClutI=0;nNumClutI<=prot.DOF-1;++nNumClutI)
	{
		for (nNumClutJ=1;nNumClutJ<=nNumClutI;++nNumClutJ)
		{
//			for (nNumAtom=/*1*/clust[nNumClutJ].origin_xoord_a;
//				 nNumAtom</*=*/clust[nNumClutJ].num_xoord_a;
//				 ++nNumAtom)
//			{
				for (i=0;i<4;++i)
				{
					for (j=0;j<4;++j)
					{
						mat[i][j] = 0.0;
					}
				}
				calc_dot_Pseduo_TransMatrix(nNumClutJ,nNumClutJ,mat);

				for (i=0;i<3;++i)
				{
					dummy_xoord[i] = clust[nNumClutI].qCOM[i];
//					dummy_xoord[i] = clust[nNumClutI].xoord_clust/*[0]*/[nNumAtom][i];
				}
				dummy_xoord[3] = 1.0;

				if (nNumClutI != nNumClutJ)
				{
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
						dotTransQ[nNumClutI][/*nNumAtom*/nNumClutJ][i]
								+= mat[i][j]
								*dummy_xoord[j];
					}
				}
			}
//		}
	}

//		for (nNumAtom=0;nNumAtom<=prot.num_atom;++nNumAtom)
		for(nNumClut=0;nNumClut<=prot.DOF-1;++nNumClut)
		{
			for (alpha=0;alpha<3;++alpha)
			{
				velo2[/*nNumAtom*/nNumClut][alpha] = 0.0;
			}
		}

		num = clust[0].num_atom_clust;
		for(nNumClutI=1;nNumClutI<=prot.DOF-1;++nNumClutI)
//		for (nNumAtom=1;nNumAtom<=prot.num_atom-num;++nNumAtom)
		{
//			for(nNumClut=1;nNumClut<=prot.DOF-1;++nNumClut)
			for(nNumClutJ=1;nNumClutJ<=prot.DOF-1;++nNumClutJ)
			{
				for (alpha=0;alpha<3;++alpha)
				{
					velo2[/*nNumAtom+num*/nNumClutI][alpha] += dotTransQ[nNumClutI][/*nNumAtom*/nNumClutJ][alpha]
									             *clust[nNumClutJ].ddihedang[0];
				}
			}
		}

//		for(nNumClut=1;nNumClut<=prot.DOF-1;++nNumClut)
//		{
//			for (nNumAtom=1;nNumAtom<=prot.num_atom;++nNumAtom)
//			{
//				for (alpha=0;alpha<3;++alpha)
//				{
//					velo2[nNumAtom+num][alpha] += dotTransQ[nNumClut][nNumAtom][alpha]
//									             *clust[nNumClut].ddihedang[0];
//				}
//			}
//		}
}

// 現在の構造でのタンパク質のポテンシャルエネルギー、
// タンパク質に及ぼす力の計算を行う関数
void calc_velo2(void) {
  int i,ii,i_c,alpha,alpha2,alpha3,j,k,n;

  int nNumClut,nNumClutI,nNumClutJ,nNumCltCoo,nNumClutdummy;
  int nNumAtom,nNumAtomtotal,nNumAtomClut;
  int nNumOrigc, nNumAbsoc;

  double mat[4][4];
  double dummy[4][4];
  double dummy_xoord[4];
  double dummy_xoord2[4];
  double *dotTransQ;
  double *velo;
  /**************************/
  /* double velo_d[100][3]; */
  /**************************/

  /******************************/
  /* printf("%p \n",dotTransQ); */
  /******************************/
  if((dotTransQ=(double *)malloc(sizeof(double)*prot.DOF*prot.num_atom*4))==NULL) {
    printf("error: malloc\n");
    exit(1);
  }
  if((velo=(double *)malloc(sizeof(double)*prot.num_atom*4))==NULL) {
      printf("error: malloc\n");
      exit(1);
  }

  for (nNumClutI=0;nNumClutI<prot.DOF;++nNumClutI) {
    //    set_pseduo_trans_matrix(nNumClutI);
  }
  /***********************/
  /* printf("a:1114\n"); */
  /***********************/

  for (nNumClut=0;nNumClut<prot.DOF;++nNumClut) {
    for (nNumAtom=0;nNumAtom<prot.num_atom;++nNumAtom)	{
      for (i=0;i<4;++i) {
	/****************************************************/
        /* printf("%d %d %d\n",nNumClut,nNumAtom,i);	    */
	/* printf("%p \n",dotTransQ);			    */
	/* printf("%p \n",velo);			    */
        /****************************************************/
	dotTransQ[nNumClut*prot.num_atom*4+nNumAtom*4+i] = 0.0;
      }
    }
  }

  for (nNumClutI=0;nNumClutI<prot.DOF;++nNumClutI) {
    nNumAtom=0;
    for (nNumClutJ=0;nNumClutJ<prot.DOF;++nNumClutJ) {
      if (/*nNumClutI != 0 && */nNumClutI<=nNumClutJ && clust[nNumClutI].join <= clust[nNumClutJ].join) {
	for (i=0;i<3;++i) {
	  for (j=0;j<3;++j) {
	    mat[i][j] = 0.0;
	  }
	}
	calc_dot_Pseduo_TransMatrix(nNumClutI,nNumClutJ,mat);
	/******************************/
        /* printf("a:1135\n");	      */
        /******************************/
      }
      else {
	for (i=0;i<3;++i) {
	  for (j=0;j<3;++j) {
	    mat[i][j] = 0.0;
	  }
	}
      }
      /***********************/
      /* printf("a:1143\n"); */
      /***********************/
      for (nNumAtomClut=0;nNumAtomClut<clust[nNumClutJ].num_atom_clust;++nNumAtomClut) {
	if (/*nNumClutI != 0 &&*/ nNumClutI<=nNumClutJ && clust[nNumClutI].join <= clust[nNumClutJ].join) {
          /***********************************************************************/
          /* if (nNumAtomClut==0) {						 */
	  /*   printf("nNumClutI=%d nNumClutJ=%d\n",nNumClutI,nNumClutJ);	 */
	  /* }									 */
          /***********************************************************************/

	  for (i=0;i<3;++i)	{
	    dummy_xoord[i] = clust[nNumClutJ].xoord_clust/*[0]*/[nNumAtomClut][i];
	  }
	  dummy_xoord[3] = 1.0;
	  
	  for (i=0;i<4;++i) {
	    for (j=0;j<4;++j) {
	      dotTransQ[nNumClutI*prot.num_atom*4+nNumAtom*4+i]+= mat[i][j]*dummy_xoord[j];
	    }
	  }
	}
	else {
	  for (i=0;i<4;++i) {
	    for (j=0;j<4;++j) {
	      dotTransQ[nNumClutI*prot.num_atom*4+nNumAtom*4+i]=0.0;
	    }
	  }
	}
	++nNumAtom;
      }
    }
  }
  /***********************/
  /* printf("a:1174\n"); */
  /***********************/

  for (nNumAtom=0;nNumAtom<prot.num_atom;++nNumAtom){
    for (alpha=0;alpha<4;++alpha){
      velo[nNumAtom*4+alpha] = 0.0;
    }
  }

  for(nNumClutJ=0;nNumClutJ<prot.DOF;++nNumClutJ) {
    for (nNumAtom=0;nNumAtom<prot.num_atom;++nNumAtom) {
      for (alpha=0;alpha<3;++alpha) {
	velo[nNumAtom*4+alpha] += dotTransQ[nNumClutJ*prot.num_atom*4+nNumAtom*4+alpha]*clust[nNumClutJ].ddihedang[0];
      }
    }
  }
  /***********************/
  /* printf("a:1189\n"); */
  /***********************/

  
  for (nNumAtom=0;nNumAtom<prot.num_atom;++nNumAtom){
    for (alpha=0;alpha<3;++alpha){
      prot.velo[nNumAtom][alpha]=velo[nNumAtom*4+alpha];
    }
  }
  /***********************/
  /* printf("a:1196\n"); */
  /***********************/

  free(dotTransQ);
  free(velo);
  /*****************************************************************************************************/
  /* for (nNumAtom=0;nNumAtom<prot.num_atom;++nNumAtom){					       */
  /*   for (alpha=0;alpha<3;++alpha){								       */
  /*     velo_d[nNumAtom][alpha]=(prot.coord[nNumAtom][alpha]-prot.old_coord[nNumAtom][alpha])/deltat; */
  /*   }											       */
  /* }												       */
  /*****************************************************************************************************/
}

void calc_velo_free_term(double vel_Term[3]) {
  int i,ii,i_c,alpha,alpha2,alpha3,j,k,n;

  int nNumClut,nNumClutI,nNumClutJ,nNumCltCoo,nNumClutdummy;
  int nNumAtom,nNumAtomtotal,nNumAtomClut;
  int nNumOrigc, nNumAbsoc;

  double mat[4][4];
  double dummy[4][4];
  double dummy_xoord[4];
  double dummy_xoord2[4];
  double *dotTransQ;
  double *velo,*velo_temp;

  int nNumAtomOrig;
  double omeg_t[3],omeg_t2[3];
  double velo_t[3],velo_t2[3];
  double avelo[3],avelo2[3];
  double vect[3],vect2[3];

  if((dotTransQ=(double *)malloc(sizeof(double)*prot.DOF*prot.num_atom*4))==NULL) {
    printf("error: malloc\n");
    exit(1);
  }
  if((velo=(double *)malloc(sizeof(double)*prot.num_atom*4))==NULL) {
      printf("error: malloc\n");
      exit(1);
  }

  for (nNumClutI=0;nNumClutI<prot.DOF;++nNumClutI) {
    //    set_pseduo_trans_matrix(nNumClutI);
  }

  for (nNumClut=0;nNumClut<prot.DOF;++nNumClut) {
    for (nNumAtom=0;nNumAtom<prot.num_atom;++nNumAtom)	{
      for (i=0;i<4;++i) {
	dotTransQ[nNumClut*prot.num_atom*4+nNumAtom*4+i] = 0.0;
      }
    }
  }

  for (nNumClutI=0;nNumClutI<prot.DOF;++nNumClutI) {
    nNumAtom=0;
    for (nNumClutJ=0;nNumClutJ<prot.DOF;++nNumClutJ) {
      if (nNumClutI<=nNumClutJ && clust[nNumClutI].join <= clust[nNumClutJ].join) {
	for (i=0;i<3;++i) {
	  for (j=0;j<3;++j) {
	    mat[i][j] = 0.0;
	  }
	}
	calc_dot_Pseduo_TransMatrix(nNumClutI,nNumClutJ,mat);
      }
      else {
	for (i=0;i<3;++i) {
	  for (j=0;j<3;++j) {
	    mat[i][j] = 0.0;
	  }
	}
      }
      for (nNumAtomClut=0;nNumAtomClut<clust[nNumClutJ].num_atom_clust;++nNumAtomClut) {
	if (nNumClutI<=nNumClutJ && clust[nNumClutI].join <= clust[nNumClutJ].join) {
	  for (i=0;i<3;++i)	{
	    dummy_xoord[i] = clust[nNumClutJ].xoord_clust/*[0]*/[nNumAtomClut][i];
	  }
	  dummy_xoord[3] = 1.0;

	  for (i=0;i<4;++i) {
	    for (j=0;j<4;++j) {
	      dotTransQ[nNumClutI*prot.num_atom*4+nNumAtom*4+i]+= mat[i][j]*dummy_xoord[j];
	    }
	  }
	}
	else {
	  for (i=0;i<4;++i) {
	    for (j=0;j<4;++j) {
	      dotTransQ[nNumClutI*prot.num_atom*4+nNumAtom*4+i]=0.0;
	    }
	  }
	}
	++nNumAtom;
      }
    }
  }
  for (nNumAtom=0;nNumAtom<prot.num_atom;++nNumAtom){
    for (alpha=0;alpha<4;++alpha){
      velo[nNumAtom*4+alpha] = 0.0;
    }
  }

  for(nNumClutJ=0;nNumClutJ<prot.DOF;++nNumClutJ) {
    for (nNumAtom=0;nNumAtom<prot.num_atom;++nNumAtom) {
      for (alpha=0;alpha<3;++alpha) {
	velo[nNumAtom*4+alpha] += dotTransQ[nNumClutJ*prot.num_atom*4+nNumAtom*4+alpha]*clust[nNumClutJ].ddihedang[0];
      }
    }
  }

  for (nNumAtom=0;nNumAtom<prot.num_atom;++nNumAtom){
    for (alpha=0;alpha<3;++alpha){
      prot.velo[nNumAtom][alpha]=velo[nNumAtom*4+alpha];
    }
  }

  nNumAtomOrig=clust[0].origin_atom_a-1;
  for (alpha=0;alpha<3;++alpha) {
    omeg_t[alpha]=vel_Term[alpha];
    velo_t[alpha]=vel_Term[alpha+3];
  }
  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      //      Rot_Term[i][j]=clust[0].trans_A_to_CN[0][i][j];
      ;
  for (i=0;i<3;++i)
    velo_t2[i]=0.0;
  for (i=0;i<3;++i)
    omeg_t2[i]=0.0;
  for (i=0;i<3;++i)
    vect2[i]=0.0;
  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      velo_t2[i]+=Rot_Term[j][i]*velo_t[j];
  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      omeg_t2[i]+=Rot_Term[i][j]*omeg_t[j];

  velo_temp=(double *)calloc(sizeof(double),prot.num_atom*3);

  for (nNumAtom=0;nNumAtom<prot.num_atom;++nNumAtom){
    for (alpha=0;alpha<3;++alpha){
      vect[alpha]=prot.coord[nNumAtom][alpha]-prot.coord[nNumAtomOrig][alpha];
    }
    for (i=0;i<3;++i) {
      vect2[i]=0.0;
    }
    for (i=0;i<3;++i) {
      for (j=0;j<3;++j) {
	//	vect2[i]+=Rot_Term[i][j]*vect[j];
	vect2[i]+=clust[0].trans_A_to_CN[0][i][j]*vect[j];
      }
    }
    for (alpha=0;alpha<3;++alpha){
      vect[alpha]=clust[0].xoord_clust/*[0]*/[nNumAtom][alpha];
    }
    //    outprod(omeg_t2,vect,avelo2);
    outprod(omeg_t ,vect2,avelo);
    for (alpha=0;alpha<3;++alpha){
      //      prot.velo[nNumAtom][alpha]+=velo_t2[alpha]+avelo2[alpha];
      prot.velo[nNumAtom][alpha]+=velo_t[alpha]+avelo[alpha];
      //      velo_temp[nNumAtom*3+alpha]=prot.velo[nNumAtom][alpha]+velo_t[alpha]+avelo[alpha];
    }
  }

  /********************************************************************************/
  /* for (nNumAtom=0;nNumAtom<prot.num_atom;++nNumAtom) {			  */
  /*   for (i=0;i<3;++i) {							  */
  /*     velo_temp[nNumAtom*3+i]=prot.velo[nNumAtom][i];			  */
  /*     prot.velo[nNumAtom][i]=0.0;						  */
  /*   }									  */
  /* }										  */
  /* for (nNumAtom=0;nNumAtom<prot.num_atom;++nNumAtom)				  */
  /*   for (i=0;i<3;++i)							  */
  /*     for (j=0;j<3;++j)							  */
  /* 	prot.velo[nNumAtom][i]+=Rot_Term[j][i]*velo_temp[nNumAtom*3+j];		  */
  /********************************************************************************/

  free(dotTransQ);
  free(velo);

}

void outprod(double v1[3],double v2[3], double v1x2[3]) {
  int alpha,alpha2;
  double v1x[3][3];

  v1x[0][0]= 0.0;
  v1x[0][1]=-v1[2];
  v1x[0][2]= v1[1];
  v1x[1][0]= v1[2];
  v1x[1][1]= 0.0;
  v1x[1][2]=-v1[0];
  v1x[2][0]=-v1[1];
  v1x[2][1]= v1[0];
  v1x[2][2]= 0.0;

  for (alpha=0;alpha<3;++alpha)
    v1x2[alpha]=0.0;

  for (alpha=0;alpha<3;++alpha)
    for (alpha2=0;alpha2<3;++alpha2)
      v1x2[alpha]+=v1x[alpha][alpha2]*v2[alpha2];
}


//// 剛体の慣性モーメントの設定
//void set_Inertia_clust2(int nNumClut)
//{
//	int i,j;
//	int nNumAtomLoca=0,nNumParent;
//
//	clust[nNumClut].Inertia_clust[0][0]=0.0;
//	clust[nNumClut].Inertia_clust[1][1]=0.0;
//	clust[nNumClut].Inertia_clust[2][2]=0.0;
//	clust[nNumClut].Inertia_clust[0][1]=0.0;
//	clust[nNumClut].Inertia_clust[0][2]=0.0;
//	clust[nNumClut].Inertia_clust[1][2]=0.0;
//
//	nNumParent = clust[nNumClut].nNumClutOfParent-1;
//
//	nNumAtomLoca =  clust[nNumParent].num_atom_clust
//				   -clust[nNumParent].terminal_atom_a[0]
//				   +clust[nNumParent].origin_atom_a;
//
//	for (i=1;i<nNumClut-nNumParent;++i)
//	{
//		nNumAtomLoca += clust[nNumParent+i].num_atom_clust;
//	}
//
//	if (clust[nNumClut].num_atom_clust > 1)
//	{
//	for(i=0; i<=clust[nNumClut].num_atom_clust; ++i)
//	{
//		clust[nNumClut].Inertia_clust[0][0]/*kg*m^2*/ += 
//					clust[nNumClut].mass_clust[/*nNumAtomLoca+*/i]/*u*//**1.660539e-27*/
//					*  (    clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][1]/*A*//**1.0e-10*/
//					       *clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][1]/*A*//**1.0e-10*/
//					     +  clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][2]/*A*//**1.0e-10*/
//					       *clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][2]/*A*//**1.0e-10*/);
//
//		clust[nNumClut].Inertia_clust[1][1]/*kg*m^2*/ += 
//					clust[nNumClut].mass_clust[/*nNumAtomLoca+*/i]/*u*//**1.660539e-27*/
//					*   (   clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][0]/*A*//**1.0e-10*/
//					       *clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][0]/*A*//**1.0e-10*/
//					      + clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][2]/*A*//**1.0e-10*/
//					       *clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][2]/*A*//**1.0e-10*/);
//
//		clust[nNumClut].Inertia_clust[2][2]/*kg*m^2*/ += 
//					clust[nNumClut].mass_clust[/*nNumAtomLoca+*/i]/*u*//**1.660539e-27*/
//					*(      clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][0]/*A*//**1.0e-10*/
//					       *clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][0]/*A*//**1.0e-10*/
//				         +  clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][1]/*A*//**1.0e-10*/
//					       *clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][1]/*A*//**1.0e-10*/ );
//
////		if (clust[nNumClut].Inertia_clust[2][2]<1.0e-70)
////		{
////			clust[nNumClut].Inertia_clust[2][2]=1.0e-47;
////		}
//		clust[nNumClut].Inertia_clust[0][1]/*kg*m^2*/ -= 
//						clust[nNumClut].mass_clust[/*nNumAtomLoca+*/i]/*u*//**1.660539e-27*/
//					      * clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][0]/*A*//**1.0e-10*/
//					      * clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][1]/*A*//**1.0e-10*/;
//
//		clust[nNumClut].Inertia_clust[0][2]/*kg*m^2*/ -= 
//						clust[nNumClut].mass_clust[/*nNumAtomLoca+*/i]/*u*//**1.660539e-27*/
//					      * clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][0]/*A*//**1.0e-10*/ 
//					      * clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][2]/*A*//**1.0e-10*/;
//
//		clust[nNumClut].Inertia_clust[1][2]/*kg*m^2*/ -= 
//					        clust[nNumClut].mass_clust[/*nNumAtomLoca+*/i]/*u*//**1.660539e-27*/
//					      * clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][1]/*A*//**1.0e-10*/
//					      * clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][2]/*A*//**1.0e-10*/;
//	}
//
//	clust[nNumClut].Inertia_clust[1][0]/*kg*m^2*/=clust[nNumClut].Inertia_clust[0][1]/*kg*m^2*/;
//	clust[nNumClut].Inertia_clust[2][0]/*kg*m^2*/=clust[nNumClut].Inertia_clust[0][2]/*kg*m^2*/;
//	clust[nNumClut].Inertia_clust[2][1]/*kg*m^2*/=clust[nNumClut].Inertia_clust[1][2]/*kg*m^2*/;
//	}
//	else
//	{
//		clust[nNumClut].Inertia_clust[0][0]/*kg*m^2*/ += 
//					clust[nNumClut].mass_clust[/*nNumAtomLoca+*/0]/*u*//**1.660539e-27*/
//					*  (  0.1*0.1/*A*//**1.0e-10*/
//					     + 0.1*0.1/*A*//**1.0e-10*/);
//
//		clust[nNumClut].Inertia_clust[1][1]/*kg*m^2*/ += 
//					clust[nNumClut].mass_clust[/*nNumAtomLoca+*/0]/*u*//**1.660539e-27*/
//					*  (  0.1*0.1/*A*//**1.0e-10*/
//					     + 0.1*0.1/*A*//**1.0e-10*/);
//
//		clust[nNumClut].Inertia_clust[2][2]/*kg*m^2*/ += 
//					clust[nNumClut].mass_clust[/*nNumAtomLoca+*/0]/*u*//**1.660539e-27*/
//					*  (  0.1*0.1/*A*//**1.0e-10*/
//					     + 0.1*0.1/*A*//**1.0e-10*/);
//
//		clust[nNumClut].Inertia_clust[0][1]/*kg*m^2*/ -= 
//						clust[nNumClut].mass_clust[/*nNumAtomLoca+*/0]/*u*//**1.660539e-27*/
//					*  (  0.1*0.1/*A*//**1.0e-10*/);
//
//		clust[nNumClut].Inertia_clust[1][0]/*kg*m^2*/ -= 
//						clust[nNumClut].mass_clust[/*nNumAtomLoca+*/0]/*u*//**1.660539e-27*/
//					*  (  0.1*0.1/*A*//**1.0e-10*/);
//
//		clust[nNumClut].Inertia_clust[0][2]/*kg*m^2*/ -= 
//						clust[nNumClut].mass_clust[/*nNumAtomLoca+*/0]/*u*//**1.660539e-27*/
//					*  (  0.1*0.1/*A*//**1.0e-10*/);
//
//		clust[nNumClut].Inertia_clust[1][2]/*kg*m^2*/ -= 
//					        clust[nNumClut].mass_clust[/*nNumAtomLoca+*/0]/*u*//**1.660539e-27*/
//					*  (  0.1*0.1/*A*//**1.0e-10*/);
//
//		clust[nNumClut].Inertia_clust[2][1]/*kg*m^2*/ -= 
//					        clust[nNumClut].mass_clust[/*nNumAtomLoca+*/0]/*u*//**1.660539e-27*/
//					*  (  0.1*0.1/*A*//**1.0e-10*/);
//	}
//
//}

