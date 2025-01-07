#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gener.h"
#include "ABA.h" // #include "ABA_multimer.h"
#include "physics.h"
#include "MD.h"
#include "force.h"

// NVTの場合、速度スケーリングを行う関数
double velocity_scaling(double Energy_kinetic_o)
{
	int i, alpha, nNumClut, nNumClut2, nNumBod;

	int nNumAtomOrigClut;

	int nNumA = 0;

	double sumMass;

	double Inertia;

	double Energy_kinetic_scaled = 0.0;

	double Energy_kinetic_tra = 0.0;

	double Energy_kinetic_rot = 0.0;

	double Energy_kinetic_clust[MAXDOF];

	// スケーリング係数
	double z;
	double beta;

	for (nNumClut=0; nNumClut<prot.DOF; ++nNumClut)
	{
		nNumAtomOrigClut = clust[nNumClut].terminal_atom_a[0];
		z = (3*prot.num_atom-1)
		      *k_B_J*T_Kelvin
		      /Energy_kinetic_o;
//		if (z<0)
//		{
//			z = -1.0*z;
//		}
		if (z == 0.0)
		{
			beta = 0.0;
		}
		else
		{
			beta = sqrt(z);
		}

		// 角速度のスケーリングを行う
		clust[nNumClut].ddihedang[0]
			  = clust[nNumClut].ddihedang[0]*beta;

		for (alpha=4; alpha<6; ++alpha)
		{
		 	clust[nNumClut].sp_velo[alpha]
			  = clust[nNumClut].sp_velo[alpha]*beta;
		}

		// 速度のスケーリングを行う
		for (alpha=0; alpha<3; ++alpha)
		{
			prot.velo[nNumAtomOrigClut][alpha]
			      = prot.velo[nNumAtomOrigClut][alpha]*beta;
		}
	}

	// スケーリング後の運動エネルギーを計算
//	for (nNumClut=0; nNumClut<prot.DOF-1; ++nNumClut)
//	{
//		nNumAtomOrigClut = clust[nNumClut].terminal_atom_a[0];

// 		sumMass = 0.0;
//		for (i=0;i<clust[nNumClut].num_atom_clust;++i)
//		{
//			sumMass += clust[nNumClut].mass_clust[i];
//		}

		// スケーリング後の並進運動エネルギー
//		for (alpha=0; alpha<3; ++alpha)
//		{
//			Energy_kinetic_tra
//			 += 0.5*(prot.velo[nNumAtomOrigClut][alpha]*1.0e-10)
//			       *([nNumAtomOrigClut][alpha]*1.0e-10)
//			       *sumMass*1.660539e-27;
//		}

		// スケーリング後の回転運動エネルギー
//		Energy_kinetic_rot
//		  += 0.5*clust[nNumClut].ddihedang[0]
//		        *clust[nNumClut].ddihedang[0]
//			    *clust[nNumClut].Inertia_clust[2][2];
//	}

	// スケーリング後の全運動エネルギー
//	Energy_kinetic_scaled
//	 = Energy_kinetic_tra + Energy_kinetic_rot;

	// 運動エネルギーを計算
	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
/*	for (nNumClut=0; nNumClut<prot.DOF-1; ++nNumClut)
	{
		if (clust[nNumClut].terminal == TERMINAL)
		{
			// 回転運動エネルギー
			Energy_kinetic_clust[nNumClut]
		  		= 0.5*clust[nNumClut].ddihedang[0]
				*clust[nNumClut].ddihedang[0]
			    *clust[nNumClut].Inertia_clust[2][2];

		}
		else if (clust[nNumClut+1].terminal == TERMINAL)
		{
			for (nNumClut2=nNumClut+1; nNumClut2<prot.DOF-1; ++nNumClut2)
			{
				Inertia += clust[nNumClut2].Inertia_clust[2][2];
			}

			// 回転運動エネルギー
			Energy_kinetic_clust[nNumClut]
		  		= 0.5*clust[nNumClut].ddihedang[0]
				*clust[nNumClut].ddihedang[0]
			    *Inertia;

			Inertia = 0.0;
 		}
		else
		{
			for (nNumClut2=nNumClut+1; nNumClut2<prot.DOF-1; ++nNumClut2)
			{
				Inertia += clust[nNumClut2].Inertia_clust[2][2];
			}

			// 回転運動エネルギー
			Energy_kinetic_clust[nNumClut]
		  		= 0.5*clust[nNumClut].ddihedang[0]
				*clust[nNumClut].ddihedang[0]
			    *Inertia;

			Inertia = 0.0;
 		}

	}

	Energy_kinetic_rot = 0.0;

	for (nNumClut=1;nNumClut<prot.DOF-1;++nNumClut)
	{
		Energy_kinetic_rot += Energy_kinetic_clust[nNumClut];
	}*/
	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////


	// 運動エネルギーを計算
	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	Energy_kinetic_o2 = 0.0;
	for(nNumClut=0; nNumClut<prot.DOF-1; ++nNumClut)
	{
		nNumA = clust[nNumClut].terminal_atom_a[0] - 1;
		for (i=0; i<clust[nNumClut].num_atom_clust; ++i)
		{
			for (alpha=0; alpha<3; ++alpha)
			{
				Energy_kinetic_o2
				+= 0.50 * (prot.velo[nNumA+i][alpha])
						* (prot.velo[nNumA+i][alpha])
						* (clust[nNumClut].mass_clust[i]
						   *1.660539e-27);
			}
		}
	}

	Energy_kinetic_o = 0.0;
	for (nNumClut=0; nNumClut<prot.DOF-1; ++nNumClut)
	{
		nNumAtomOrigClut = clust[nNumClut].origin_atom_a;

 		sumMass = 0.0;
		for (i=0;i<clust[nNumClut].num_atom_clust;++i)
		{
			sumMass += clust[nNumClut].mass_clust[i];
		}

		// 並進運動エネルギー
		for (alpha=4; alpha<6; ++alpha)
		{
			Energy_kinetic_tra
			 += 0.5*(clust[nNumClut].sp_velo[alpha]/**1.0e-10*/)
			       *(clust[nNumClut].sp_velo[alpha]/**1.0e-10*/)
			       *sumMass*1.660539e-27;
		}

		// 回転運動エネルギー
		Energy_kinetic_rot
		  += 0.5*clust[nNumClut].ddihedang[0]
				*clust[nNumClut].ddihedang[0]
			    *clust[nNumClut].Inertia_clust[2][2];
	}

	Energy_kinetic_o = Energy_kinetic_tra + Energy_kinetic_rot;

	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////

//	Energy_kinetic_scaled
//	 = Energy_kinetic_rot;

	Energy_kinetic_scaled
	 = Energy_kinetic_o;

	return Energy_kinetic_scaled;
}
