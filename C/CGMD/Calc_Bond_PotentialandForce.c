#include <stdio.h>
#include <math.h>
#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "MD.h"
#include "force.h"
#include "physics.h"

// Œ‹‡’·‚ÌLk‹y‚ÑŒ‹‡Šp‚Ì•ÏŠp‚ÌƒGƒlƒ‹ƒM[
void Calc_Bond_PotentialandForce(void)
{
	int alpha;
	int nNumClut;

	potential_pro.bond_length = 0.0;
	potential_pro.bond_angle = 0.0;

	for (nNumClut=0;nNumClut<prot.DOF;++nNumClut)
	{
		sub_Calc_Bond_Potential(nNumClut);
		sub_Calc_Bond_Force(nNumClut);
	}

	potential_pro.bond_length += clust[nNumClust].f_c.bond;
	potential_pro.bond_angle += clust[nNumClust].f_c.angle;
}

// Œ‹‡’·‚ÌLk‚ÌƒGƒlƒ‹ƒM[
void sub_Calc_Bond_length_Potential(int nNumClut)
{
	int alpha;
	int nNumBond_Clust;
	double delta_length;

	clust[nNumClust].f_c.bond = 0.0;

	for (nNumBond_Clust=0;nNumBond_Clust<clust[nNumClust].numBond;
                                                 ++nNumBond_Clust)
	{
		// Œ‹‡’··‚ÌŒvZ
		delta_length =  bond_length[nNumBond_Clust]
		              - clust[nNumClust].f_p_clust.bond_length_refrence;

		// Œ‹‡’·‚ÌƒGƒlƒ‹ƒM[‚ÌŒvZ
		clust[nNumClust].f_c.bond += clust[nNumClust].f_p_clust.k_bond_2
		                            *(delta_length)*(delta_length);
	}
}

// Œ‹‡Šp‚Ì•ÏŠp‚ÌƒGƒlƒ‹ƒM[
void sub_Calc_Bond_angle_Force(int nNumClut)
{
	int alpha;
	int nNumBond_Clust;
	double delta_angle;

	clust[nNumClust].f_c.angle = 0.0;

	for (nNumBond_Clust=0;nNumBond_Clust<clust[nNumClust].numBond;
                                                 ++nNumBond_Clust)
	{
		// Šp“x·‚ÌŒvZ
		delta_angle =  bond_angle[nNumBond_Clust]
		              - clust[nNumClust].f_p_clust.bond_angle_refrence;

		// •ÏŠp‚ÌƒGƒlƒ‹ƒM[‚ÌŒvZ
		clust[nNumClust].f_c.angle += clust[nNumClust].f_p_clust.k_angle_2
		                            *(delta_angle)*(delta_angle);
	}
}
