#include <stdio.h>
#include <math.h>

#include "Vis_MD.h"

void TransAbsoToLocal(int nNumAtom1, int nNumAtom2, int nSumNumAtom)
{
	int nNumAtom;
	int alpha;

	int kkk, kkkk, lll, llll, nn;

	double CN_A[3];
	double HN_A[3];
	double CN_1_A[3];

	double Mat[/*MAXNATOM*/10][3];

	double ii[3];
	double jj[3];
	double kk[3];

	double jj_x;
	double jj_y;
	double jj_z;

	double d1=0.0;
	double d2=0.0;
	double sn=0.0;
	double cs=0.0;

	for(alpha=0; alpha<3; ++alpha)
	{
		CN_A[alpha]   = prot.atm[nNumAtom1].coord[alpha];
		HN_A[alpha]   = prot.atm[nNumAtom1+1].coord[alpha];
		CN_1_A[alpha] = prot.atm[nNumAtom2].coord[alpha];
	}

	for(alpha=0; alpha<3; ++alpha)
	{
		kk[alpha] = CN_1_A[alpha]-CN_A[alpha];
		d1 += (CN_1_A[alpha]-CN_A[alpha])*(CN_1_A[alpha]-CN_A[alpha]);
		d2 += (HN_A[alpha]-CN_A[alpha])*(HN_A[alpha]-CN_A[alpha]);
		cs += (CN_1_A[alpha]-CN_A[alpha])*(HN_A[alpha]-CN_A[alpha]);
	}

	d1 = sqrt(d1);
	d2 = sqrt(d2);
	cs = -cs/(d1*d2);
	sn = 1.0-cs*cs;
	sn = sqrt(sn);

	for(alpha=0; alpha<3; ++alpha)
	{
		kk[alpha]=kk[alpha]/d1;
	}

	jj[0]=((CN_1_A[1]-CN_A[1])*(HN_A[2]-CN_A[2])-(CN_1_A[2]-CN_A[2])*(HN_A[1]-CN_A[1]))/(d1*d2*sn);
	jj[1]=((CN_1_A[2]-CN_A[2])*(HN_A[0]-CN_A[0])-(CN_1_A[0]-CN_A[0])*(HN_A[2]-CN_A[2]))/(d1*d2*sn);
	jj[2]=((CN_1_A[0]-CN_A[0])*(HN_A[1]-CN_A[1])-(CN_1_A[1]-CN_A[1])*(HN_A[0]-CN_A[0]))/(d1*d2*sn);

	ii[0]=(jj[1]*(CN_1_A[2]-CN_A[2])-jj[2]*(CN_1_A[1]-CN_A[1]))/d1;
	ii[1]=(jj[2]*(CN_1_A[0]-CN_A[0])-jj[0]*(CN_1_A[2]-CN_A[2]))/d1;
	ii[2]=(jj[0]*(CN_1_A[1]-CN_A[1])-jj[1]*(CN_1_A[0]-CN_A[0]))/d1;

	for(nNumAtom=0; nNumAtom < nSumNumAtom; ++nNumAtom)
	{
		for(alpha=0;alpha<3;++alpha)
		{
			Mat[nNumAtom][alpha] =   prot.atm[nNumAtom1+nNumAtom].coord[alpha]
			                       - CN_A[alpha];
		}
	}

	for(alpha=0;alpha<3;++alpha)
	{
		MatTransAbsoToLocal[0][alpha]=ii[alpha];
		MatTransAbsoToLocal[1][alpha]=jj[alpha];
		MatTransAbsoToLocal[2][alpha]=kk[alpha];
	}

	for(nNumAtom=0; nNumAtom < nSumNumAtom; ++nNumAtom)
	{
		xoord[nNumAtom][0] =   MatTransAbsoToLocal[0][0]*Mat[nNumAtom][0]
				    		     + MatTransAbsoToLocal[0][1]*Mat[nNumAtom][1]
				    		     + MatTransAbsoToLocal[0][2]*Mat[nNumAtom][2];

		xoord[nNumAtom][1] =   MatTransAbsoToLocal[1][0]*Mat[nNumAtom][0]
						         + MatTransAbsoToLocal[1][1]*Mat[nNumAtom][1]
						         + MatTransAbsoToLocal[1][2]*Mat[nNumAtom][2];

		xoord[nNumAtom][2] =   MatTransAbsoToLocal[2][0]*Mat[nNumAtom][0]
						  		 + MatTransAbsoToLocal[2][1]*Mat[nNumAtom][1]
						  		 + MatTransAbsoToLocal[2][2]*Mat[nNumAtom][2];
	}
}
