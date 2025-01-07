#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gener.h"
#include "ABA.h" // #include "ABA_multimer.h"
#include "physics.h"
#include "MD.h"
#include "force.h"
#include "BD.h"

// –€CŒW”‚Ìİ’è‚ğs‚¤ŠÖ”
void set_friction_tensor(void)
{
	int nNumClut;

	for (nNumClut = 0; nNumClut < prot.DOF; ++nNumClut)
	{
		// –€CŒW”‚Ìİ’è‚ğs‚¤
		sub_set_friction_tensor(nNumClut);
	}
}

// –€CŒW”‚Ìİ’è‚Ì•â•‚ğs‚¤ŠÖ”
void sub_set_friction_tensor(int nNumClut)
{
	double radius;
	int alpha;
	int nNumAtom;

	// ”¼Œa‚Ìæ“¾‚ğs‚¤
//	radius = pick_clust_radius(nNumClut);
	radius = 1.0;

	for (nNumAtom=0;nNumAtom<clust[nNumClut].num_atom_clust;++nNumClut)
	{
		// –€CŒW”‚Ìİ’è‚Ì•â•‚ğs‚¤(•Ài)
	  //		clust[nNumClut].friction_tensor_tra/*[nNumAtom]*//*kgA^2/s*/
	  //             = 6*PI*visco_wat
	  //			    /clust[nNumClut].InertiaMatrix[3][3]/*kgm/s*/
	  //			    *radius/*A*/*1.0e-10/*A^2/m^2*/;

		// –€CŒW”‚Ìİ’è‚Ì•â•‚ğs‚¤(‰ñ“])
		//clust[nNumClut].friction_tensor_rot/*[nNumAtom]*//*kgA^2/s*/
	  //           = 8*PI*visco_wat
	  //			    /clust[nNumClut].InertiaMatrix[3][3]/*kgm/s*/
	  //                *radius*radius*radius*1.0e-30/*A^2/m^2*//*A*/;
	}
}

// „‘Ì‚Ì‹…‹ß—‚ğs‚¤ŠÖ”
double pick_clust_radius(int nNumClut)
{
	double radius;

	if (clust[nNumClut].num_atom_clust < 2)
	{
		radius = sub_pick_clust_radius_1(nNumClut);
	}
	else if (clust[nNumClut].num_atom_clust < 3)
	{
		radius = sub_pick_clust_radius_2(nNumClut);
	}
	else
	{
		radius = sub_pick_clust_radius_3(nNumClut);
	}

	return radius;
}

// „‘Ì‚Ì‹…‹ß—‚Ì•â•‚ğs‚¤ŠÖ”_1
double sub_pick_clust_radius_1(int nNumClut)
{
	return 1.0/*A*/;
}

// „‘Ì‚Ì‹…‹ß—‚Ì•â•‚ğs‚¤ŠÖ”_2
double sub_pick_clust_radius_2(int nNumClut)
{
	int alpha;
	double radius = 0.0;

	for (alpha=0;alpha<3;++alpha)
	{
		radius += (clust[nNumClut].xoord_clust/*[0]*/[0][alpha]
	              -clust[nNumClut].xoord_clust/*[0]*/[1][alpha])
				 *(clust[nNumClut].xoord_clust/*[0]*/[0][alpha]
			      -clust[nNumClut].xoord_clust/*[0]*/[1][alpha]);
	}

	radius = sqrt(radius);

	return radius/2.0/*A*/;
}

// „‘Ì‚Ì‹…‹ß—‚Ì•â•‚ğs‚¤ŠÖ”_3
double sub_pick_clust_radius_3(int nNumClut)
{
/*	int i;
	int alpha;
	double originofcube[3];
	double radius;
	double S;

	for (i=0;i<clust[nNumClut].num_atom_clust;++i)
	{
		for (alpha=0;alpha<3;++alpha)
		{
			originofcube[alpha]
			 += clust[nNumClut].xoord_clust[0][i][alpha];
		}
	}

	for (alpha=0; alpha<3; ++alpha)
	{
		originofcube[alpha]
		 += originofcube[alpha]
		   /clust[nNumClut].num_atom_clust;
	}

	for (i=0;i<clust[nNumClut].num_atom_clust;++i)
	{
		S = 0.0;
		for (alpha=0;alpha<3;++alpha)
		{
			S += (  clust[nNumClut].xoord_clust[0][i][alpha]
				  - originofcube[alpha])
				 *(  clust[nNumClut].xoord_clust[0][i][alpha]
				  - originofcube[alpha]);
		}

		if (TargetS < S)
		{
			TargetS = S;
		}
	}

	radius = TargetS;

	return radius;*/

	int i;
	int n;
	int alpha;

	double S;
	double TargetS = 0.0;
	double TargetS_new = 0.0;
	double radius;
	double DeltaS;
	double originofcube[3];
	double originofcube_new[3];

	for (alpha=0;alpha<3;++alpha)
	{
		originofcube[alpha] = 0.0;
	}

	// ‹…‚Ì’†S‚ğ’Tõ
	for (n=0;n<1000;++n)
	{
		for (alpha=0;alpha<3;++alpha)
		{
			originofcube_new[alpha] = 0.0;
		}

		// ‹…‚Ì”¼Œa‚ÌŒó•â
		for (i=0;i<clust[nNumClut].num_atom_clust;++i)
		{
			S = 0.0;
			for (alpha=0;alpha<3;++alpha)
			{
				S += (  clust[nNumClut].xoord_clust/*[0]*/[i][alpha]
				      - originofcube[alpha])
				    *(  clust[nNumClut].xoord_clust/*[0]*/[i][alpha]
				      - originofcube[alpha]);
			}

			if (TargetS < S)
			{
				TargetS = S;
			}
		}

		for (alpha = 0; alpha < 3; ++alpha)
		{
			originofcube_new[alpha] = originofcube[alpha];
		}

		// —”‚Å“®‚©‚·
		mc_move(originofcube_new, 0.5);

		// ‹…‚Ì”¼Œa‚ÌŒó•â_2
		for (i=0;i<clust[nNumClut].num_atom_clust;++i)
		{
			S = 0.0;
			for (alpha=0;alpha<3;++alpha)
			{
				S += (  clust[nNumClut].xoord_clust/*[0]*/[i][alpha]
				      - originofcube_new[alpha])
				    *(  clust[nNumClut].xoord_clust/*[0]*/[i][alpha]
				      - originofcube_new[alpha]);
			}

			if (TargetS_new < S)
			{
				TargetS_new = S;
			}
		}

		DeltaS = TargetS_new - TargetS;

		// V‚µ‚¢”¼Œa‚ğÌ—p‚©
		if (DeltaS < 0.0)
		{
			originofcube[alpha] = originofcube_new[alpha];
		}
		else if (DeltaS > -0.1 && DeltaS < 0.0)
		{
			break;
		}
	}

	radius = TargetS;

	return radius/*A*/;
}

// ’†S‚Ì’Tõ‚ğs‚¤ŠÖ”
void mc_move(double origin[3], double step_limit)
{
	int alpha;

	for (alpha = 0;alpha < 3;++alpha)
	{
		origin[alpha] += randum_delta_q(step_limit);
	}
}

// ”÷¬‚È•ÏˆÊ‚Ì”­¶‚ğs‚¤ŠÖ”
double randum_delta_q(double step_limit)
{
	double u;
	double delta_q;

	u = (10.0/(RAND_MAX + 10.0)) * rand();
	delta_q = u*step_limit;

	return delta_q;
}
