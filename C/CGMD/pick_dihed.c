#include <stdio.h>
#include <math.h>
#include "gener.h"
#include "ABA.h" // #include "ABA_multimer.h"
#include "physics.h"
#include "MD.h"

void sub_pick_dihed(int nNumClutOrigin, int nNumClutTarget);
double sub_pick_dihed2(int nNumClutOrigin, int nNumClutTarget);



// “ñ–ÊŠp‚Ìæ“¾‚ğs‚¤ŠÖ”
void pick_dihed(int nNumClut)
{
	int origin_of_this_branch;
	int nNumClutOfParent;

	if (clust[nNumClut].num_branch > 1)
	{
		origin_of_this_branch = nNumClut;
	}


	nNumClutOfParent = clust[nNumClut].nNumClutOfParent-1;


	// I’[‚Ì„‘Ì‚Ìê‡
//	if(clust[nNumClut].terminal == TERMINAL)
//	{
//		sub_pick_dihed(nNumClut, origin_of_this_branch);
//	}
	// ’Êí‚Ì„‘Ì‚Ìê‡
//	else
//	{
		sub_pick_dihed(nNumClut, nNumClutOfParent);
//	}
}

/*********************************************************************************/

// “ñ–ÊŠp‚Ìæ“¾‚Ì•â•‚ğs‚¤ŠÖ”
void sub_pick_dihed(int nNumClutOrigin, int nNumClutTarget)
{
	int alpha;

	int nNumAtomOrigin;
	int nNumAtomTarget;
	int nNumAtomTargettwo;

	double CN_A[3];
	double HN_A[3];
	double CN_1_A[3];
	double HN_1_A[3];

	double cc1[3];
	double cc2[3];

	double d1=0.0;
	double d2=0.0;
	double d4=0.0;
	double cs=0.0;

	double det;

	double theta;

	double sintheta=0.0;
	double sinphi=0.0;

	// Œ´q”Ô†‚Ìæ“¾
	nNumAtomOrigin = clust[nNumClutOrigin].origin_atom_a-1;
	nNumAtomTarget = clust[nNumClutTarget].terminal_atom_a[0]-1;
//	nNumAtomTarget = clust[nNumClutTarget].origin_atom_a-1;

	if (nNumAtomTarget+1 == nNumAtomOrigin)
	{
		nNumAtomTargettwo = nNumAtomTarget - 1;
	}
	else
	{
		nNumAtomTargettwo = nNumAtomTarget + 1;
	}

	// Œ´qÀ•W‚Ìæ“¾
	for(alpha=0;alpha<3;++alpha)
	{
		CN_A[alpha]=prot.coord[nNumAtomOrigin][alpha];
		HN_A[alpha]=prot.coord[nNumAtomOrigin + 1][alpha];
		CN_1_A[alpha]=prot.coord[nNumAtomTarget][alpha];
		HN_1_A[alpha]=prot.coord[nNumAtomTargettwo][alpha];
	}

	// Œ‹‡’·‚Ìæ“¾_1
	for(alpha=0;alpha<3;++alpha)
	{
		d1 += (CN_1_A[alpha]-CN_A[alpha])*(CN_1_A[alpha]-CN_A[alpha]);
		d2 += (HN_A[alpha]-CN_A[alpha])*(HN_A[alpha]-CN_A[alpha]);
		d4 += (HN_1_A[alpha]-CN_1_A[alpha])*(HN_1_A[alpha]-CN_1_A[alpha]);
	}

	// Œ‹‡’·‚Ìæ“¾_2
	d1 = sqrt(d1);
	d2 = sqrt(d2);
	d4 = sqrt(d4);

	// Šp(CN_1_A-CN_A-HN_A)‚ÆŠp(CN_1_A-HN_A-CN_A)‚Ìæ“¾_1
	for(alpha=0;alpha<3;++alpha)
	{
		sintheta += (CN_1_A[alpha]-CN_A[alpha])*(HN_A[alpha]-CN_A[alpha]);
		sinphi += (CN_A[alpha]-CN_1_A[alpha])*(HN_1_A[alpha]-CN_1_A[alpha]);
	}

	// Šp(CN_1_A-CN_A-HN_A)‚ÆŠp(CN_1_A-HN_A-CN_A)‚Ìæ“¾_2
	sintheta = sintheta/(d1*d2);
	sintheta = 1.0 - sintheta*sintheta;
	if (sintheta != 0)
		sintheta = sqrt(sintheta);
	else 
		sintheta = 0.0;

	sinphi = sinphi/(d1*d4);
	sinphi = 1.0 - sinphi*sinphi;
	if (sinphi != 0)
		sinphi = sqrt(sinphi);
	else 
		sinphi = 0.0;

	// (CN_1_A-CN_A)x(HN_A-CN_A)
	cc1[0]=((CN_1_A[1]-CN_A[1])*(HN_A[2]-CN_A[2])
	       -(CN_1_A[2]-CN_A[2])*(HN_A[1]-CN_A[1]))
	       /(d1*d2*sintheta);
	cc1[1]=((CN_1_A[2]-CN_A[2])*(HN_A[0]-CN_A[0])
	       -(CN_1_A[0]-CN_A[0])*(HN_A[2]-CN_A[2]))
	       /(d1*d2*sintheta);
	cc1[2]=((CN_1_A[0]-CN_A[0])*(HN_A[1]-CN_A[1])
	       -(CN_1_A[1]-CN_A[1])*(HN_A[0]-CN_A[0]))
	       /(d1*d2*sintheta);

	// (CN_A-CN_1_A)x(HN_1_A-CN_1_A)
	cc2[0]=((CN_A[1]-CN_1_A[1])*(HN_1_A[2]-CN_1_A[2])
	       -(CN_A[2]-CN_1_A[2])*(HN_1_A[1]-CN_1_A[1]))
	       /(d1*d4*sinphi);
	cc2[1]=((CN_A[2]-CN_1_A[2])*(HN_1_A[0]-CN_1_A[0])
	       -(CN_A[0]-CN_1_A[0])*(HN_1_A[2]-CN_1_A[2]))
	       /(d1*d4*sinphi);
	cc2[2]=((CN_A[0]-CN_1_A[0])*(HN_1_A[1]-CN_1_A[1])
	       -(CN_A[1]-CN_1_A[1])*(HN_1_A[0]-CN_1_A[0]))
	       /(d1*d4*sinphi);

	for(alpha=0;alpha<3;++alpha)
	{
		cs += cc1[alpha]*cc2[alpha];
	}

	if (cs < -1.0)
	{
		cs = -1.0;
	}
	else if (cs > 1.0)
	{
		cs = 1.0;
	}

//			det = cc1[0]*(HN_A[0]-HN_1_A[0])
//			     +cc1[1]*(HN_A[1]-HN_1_A[1])
// 		     	 +cc1[2]*(HN_A[2]-HN_1_A[2]);
//
//			if (det < 0)
//			{
//				theta/*rad*/ = PI/*rad*/ - acos(cs)/*rad*/;
//			}
//			else
//			{
//				theta/*rad*/ = PI/*rad*/ + acos(cs)/*rad*/;
//			}

			det = cc1[0]*(HN_A[0]-CN_A[0])
			     +cc1[1]*(HN_A[1]-CN_A[1])
 		     	 +cc1[2]*(HN_A[2]-CN_A[2]);

//			if (det < 0)
//			{
//				theta/*rad*/ = PI/*rad*/ + acos(cs)/*rad*/;
//			}
//			else
//			{
				theta/*rad*/ = PI/*rad*/ - acos(cs)/*rad*/;
//			}

	// “ñ–ÊŠp‚Ì‘ã“ü
	clust[nNumClutOrigin].dihedang[0] = theta;
/*********************************************************************************/
}

// “ñ–ÊŠp‚Ìæ“¾‚ğs‚¤ŠÖ”
void pick_dihed2(int nNumClut)
{
	int origin_of_this_branch;
	int nNumClutOfParent;
	double theta;
	FILE *outtest2;

	if (clust[nNumClut].num_branch > 1)
	{
		origin_of_this_branch = nNumClut;
	}


	nNumClutOfParent = clust[nNumClut].nNumClutOfParent-1;


	// I’[‚Ì„‘Ì‚Ìê‡
//	if(clust[nNumClut].terminal == TERMINAL)
//	{
//		sub_pick_dihed(nNumClut, origin_of_this_branch);
//	}
	// ’Êí‚Ì„‘Ì‚Ìê‡
//	else
//	{
		theta = sub_pick_dihed2(nNumClut, nNumClutOfParent);
//	}
//		if ((outtest2=fopen("dihedang2.out","a")) == NULL)
//		{
//			printf("in\n");
//			exit(1);
//		}
//		for(nNumClut=1; nNumClut<prot.DOF; ++nNumClut)
//		{
//			fprintf(outtest2, "%d %e %e %e \n",nNumClut,clust[nNumClut].dihedang[0],clust[nNumClut].ddihedang[0],clust[nNumClut].correct_dihedang[2]*2.0/deltat/deltat);
//		}
//		fclose(outtest2);
}

/*********************************************************************************/

// “ñ–ÊŠp‚Ìæ“¾‚Ì•â•‚ğs‚¤ŠÖ”
double sub_pick_dihed2(int nNumClutOrigin, int nNumClutTarget)
{
	int alpha,i;

	int nNumAtomOrigin;
	int nNumAtomTarget;
	int nNumAtomTargettwo;

	double CN_A[3];
	double HN_A[3];
	double CN_1_A[3];
	double HN_1_A[3];
	double vect[3];

	double cc1[3];
	double cc2[3];

	double d1=0.0;
	double d2=0.0;
	double d4=0.0;
	double cs=0.0;

	double det[3];

	double theta;

	double sintheta=0.0;
	double sinphi=0.0;

	// Œ´q”Ô†‚Ìæ“¾
	nNumAtomOrigin = clust[nNumClutOrigin].origin_atom_a-1;
	nNumAtomTarget = clust[nNumClutTarget].terminal_atom_a[0]-1;
//	nNumAtomTarget = clust[nNumClutTarget].origin_atom_a-1;

	if (nNumAtomTarget+1 == nNumAtomOrigin)
	{
		nNumAtomTargettwo = nNumAtomTarget - 1;
	}
	else
	{
		nNumAtomTargettwo = nNumAtomTarget + 1;
	}

	// Œ´qÀ•W‚Ìæ“¾
	for(alpha=0;alpha<3;++alpha)
	{
		CN_A[alpha]=prot.coord[nNumAtomOrigin][alpha];
		HN_A[alpha]=prot.coord[nNumAtomOrigin + 1][alpha];
		CN_1_A[alpha]=prot.coord[nNumAtomTarget][alpha];
		HN_1_A[alpha]=prot.coord[nNumAtomTargettwo][alpha];
	}

	// Œ‹‡’·‚Ìæ“¾_1
	for(alpha=0;alpha<3;++alpha)
	{
		d1 += (CN_1_A[alpha]-CN_A[alpha])*(CN_1_A[alpha]-CN_A[alpha]);
		d2 += (HN_A[alpha]-CN_A[alpha])*(HN_A[alpha]-CN_A[alpha]);
		d4 += (HN_1_A[alpha]-CN_1_A[alpha])*(HN_1_A[alpha]-CN_1_A[alpha]);
	}

	// Œ‹‡’·‚Ìæ“¾_2
	d1 = sqrt(d1);
	d2 = sqrt(d2);
	d4 = sqrt(d4);

	// Šp(CN_1_A-CN_A-HN_A)‚ÆŠp(CN_1_A-HN_A-CN_A)‚Ìæ“¾_1
	for(alpha=0;alpha<3;++alpha)
	{
		sintheta += (CN_1_A[alpha]-CN_A[alpha])*(HN_A[alpha]-CN_A[alpha]);
		sinphi += (CN_A[alpha]-CN_1_A[alpha])*(HN_1_A[alpha]-CN_1_A[alpha]);
	}

	// Šp(CN_1_A-CN_A-HN_A)‚ÆŠp(CN_1_A-HN_A-CN_A)‚Ìæ“¾_2
	sintheta = sintheta/(d1*d2);
	sintheta = 1.0 - sintheta*sintheta;
	if (sintheta != 0)
		sintheta = sqrt(sintheta);
	else 
		sintheta = 0.0;

	sinphi = sinphi/(d1*d4);
	sinphi = 1.0 - sinphi*sinphi;
	if (sinphi != 0)
		sinphi = sqrt(sinphi);
	else 
		sinphi = 0.0;

	// (CN_1_A-CN_A)x(HN_A-CN_A)
	cc1[0]=((CN_1_A[1]-CN_A[1])*(HN_A[2]-CN_A[2])
	       -(CN_1_A[2]-CN_A[2])*(HN_A[1]-CN_A[1]))
	       /(d1*d2*sintheta);
	cc1[1]=((CN_1_A[2]-CN_A[2])*(HN_A[0]-CN_A[0])
	       -(CN_1_A[0]-CN_A[0])*(HN_A[2]-CN_A[2]))
	       /(d1*d2*sintheta);
	cc1[2]=((CN_1_A[0]-CN_A[0])*(HN_A[1]-CN_A[1])
	       -(CN_1_A[1]-CN_A[1])*(HN_A[0]-CN_A[0]))
	       /(d1*d2*sintheta);

	// (CN_A-CN_1_A)x(HN_1_A-CN_1_A)
	cc2[0]=((CN_A[1]-CN_1_A[1])*(HN_1_A[2]-CN_1_A[2])
	       -(CN_A[2]-CN_1_A[2])*(HN_1_A[1]-CN_1_A[1]))
	       /(d1*d4*sinphi);
	cc2[1]=((CN_A[2]-CN_1_A[2])*(HN_1_A[0]-CN_1_A[0])
	       -(CN_A[0]-CN_1_A[0])*(HN_1_A[2]-CN_1_A[2]))
	       /(d1*d4*sinphi);
	cc2[2]=((CN_A[0]-CN_1_A[0])*(HN_1_A[1]-CN_1_A[1])
	       -(CN_A[1]-CN_1_A[1])*(HN_1_A[0]-CN_1_A[0]))
	       /(d1*d4*sinphi);

	for(alpha=0;alpha<3;++alpha)
	{
		cs += cc1[alpha]*cc2[alpha];
	}

	if (cs < -1.0)
	{
		cs = -1.0;
	}
	else if (cs > 1.0)
	{
		cs = 1.0;
	}

//			det = cc1[0]*(HN_A[0]-HN_1_A[0])
//			     +cc1[1]*(HN_A[1]-HN_1_A[1])
// 		     	 +cc1[2]*(HN_A[2]-HN_1_A[2]);
//
//			if (det < 0)
//			{
//				theta/*rad*/ = PI/*rad*/ - acos(cs)/*rad*/;
//			}
//			else
//			{
//				theta/*rad*/ = PI/*rad*/ + acos(cs)/*rad*/;
//			}

//			det = cc1[0]*(HN_A[0]-CN_A[0])
//			     +cc1[1]*(HN_A[1]-CN_A[1])
// 		     	 +cc1[2]*(HN_A[2]-CN_A[2]);

//			if (det < 0)
//			{
//				theta/*rad*/ = PI/*rad*/ + acos(cs)/*rad*/;
//			}
//			else
//			{
//				theta/*rad*/ = PI/*rad*/ - acos(cs)/*rad*/;
//			}


			for (i=0;i<3;++i)
			{
				vect[i] = CN_A[i] + cc2/*1*/[i];
			}

//			if (det[0]*vect[0]+det[1]*vect[1]+det[2]*vect[2]-det[3] < 0.0)
//			{
//				if (det[0]*HN_A[0]+det[1]*HN_A[1]+det[2]*HN_A[2]-det[3] > /*<*/ 0.0)
//				{
//					theta/*rad*/ = PI/*rad*/ + acos(cs)/*rad*/;
//				}
//				else
//				{
//					theta/*rad*/ = PI/*rad*/ - acos(cs)/*rad*/;
//				}
//			}
//			else
//			{
//				if (det[0]*HN_A[0]+det[1]*HN_A[1]+det[2]*HN_A[2]-det[3] < /*>*/ 0.0)
//				{
//					theta/*rad*/ = PI/*rad*/ + acos(cs)/*rad*/;
//				}
//				else
//				{
					theta/*rad*/ = PI/*rad*/ - acos(cs)/*rad*/;
//				}
//			}

	// “ñ–ÊŠp‚Ì‘ã“ü
	return theta;
/*********************************************************************************/
}
