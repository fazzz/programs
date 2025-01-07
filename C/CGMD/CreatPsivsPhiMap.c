#include <stdio.h>
#include <math.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "physics.h"
#include "MD.h"

#define PSI 0
#define PHI 1

double sub_CreatPsivsPhiMap(int nNumAtom_N, int nNumAtom_C, int flag);

void CreatPsivsPhiMap(FILE *outputMap)
{
	int i, nNumAtom_N, nNumAtom_C;

	double phi=0.0, psi=0.0;

	for ( i=0; i<prot.nNumPeptide; ++i )
	{
		nNumAtom_N = nNumAtomPeptide_N[i];
		nNumAtom_C = nNumAtomPeptide_C[i];

		psi = sub_CreatPsivsPhiMap(nNumAtom_N, nNumAtom_C, PSI);
		phi = sub_CreatPsivsPhiMap(nNumAtom_N, nNumAtom_C, PHI);

		if (psi > PI)
		{
			psi -= 2*PI;
		}

		if (phi > PI)
		{
			phi -= 2*PI;
		}

		fprintf(outputMap, "%6.3lf ", psi*180/PI);
		fprintf(outputMap, "%6.3lf", phi*180/PI);
		fprintf(outputMap, "\n");
	}
	fprintf(outputMap, "\n\n");
}

double sub_CreatPsivsPhiMap(int nNumAtom_N, int nNumAtom_C, int flag)
{
	int i,j,kkk,l,lll,m,n,kkkk,llll;

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
	double Theta;

	double sintheta=0.0;
	double sinphi=0.0;

	double anglepsi = 0.0;

	if (flag == 0)
	{
		kkk = nNumAtom_N;
		lll = nNumAtom_N + 2;

		kkkk = nNumAtom_N + 1;
		llll = nNumAtom_C;
	}
	else
	{
		kkk = nNumAtom_N + 2;
		lll = nNumAtom_C;

		kkkk = nNumAtom_N;
		llll = nNumAtom_C + 1;
	}

	for(i=0;i<3;++i)
	{
		CN_A[i]=prot.coord[kkk][i];
		HN_A[i]=prot.coord[kkkk][i];
		CN_1_A[i]=prot.coord[lll][i];
		HN_1_A[i]=prot.coord[llll][i];
	}

	for(i=0;i<3;++i)
	{
		d1 += (CN_1_A[i]-CN_A[i])*(CN_1_A[i]-CN_A[i]);
		d2 += (HN_A[i]-CN_A[i])*(HN_A[i]-CN_A[i]);
		d4 += (HN_1_A[i]-CN_1_A[i])*(HN_1_A[i]-CN_1_A[i]);
	}

	d1 = sqrt(d1);
	d2 = sqrt(d2);
	d4 = sqrt(d4);

	for(i=0;i<3;++i)
	{
		sintheta += (CN_1_A[i]-CN_A[i])*(HN_A[i]-CN_A[i]);
		sinphi += (CN_A[i]-CN_1_A[i])*(HN_1_A[i]-CN_1_A[i]);
	}

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


	cc1[0]=((CN_1_A[1]-CN_A[1])*(HN_A[2]-CN_A[2])-(CN_1_A[2]-CN_A[2])*(HN_A[1]-CN_A[1]))/(d1*d2*sintheta);
	cc1[1]=((CN_1_A[2]-CN_A[2])*(HN_A[0]-CN_A[0])-(CN_1_A[0]-CN_A[0])*(HN_A[2]-CN_A[2]))/(d1*d2*sintheta);
	cc1[2]=((CN_1_A[0]-CN_A[0])*(HN_A[1]-CN_A[1])-(CN_1_A[1]-CN_A[1])*(HN_A[0]-CN_A[0]))/(d1*d2*sintheta);

	cc2[0]=((CN_A[1]-CN_1_A[1])*(HN_1_A[2]-CN_1_A[2])-(CN_A[2]-CN_1_A[2])*(HN_1_A[1]-CN_1_A[1]))/(d1*d4*sinphi);
	cc2[1]=((CN_A[2]-CN_1_A[2])*(HN_1_A[0]-CN_1_A[0])-(CN_A[0]-CN_1_A[0])*(HN_1_A[2]-CN_1_A[2]))/(d1*d4*sinphi);
	cc2[2]=((CN_A[0]-CN_1_A[0])*(HN_1_A[1]-CN_1_A[1])-(CN_A[1]-CN_1_A[1])*(HN_1_A[0]-CN_1_A[0]))/(d1*d4*sinphi);


	for(i=0;i<3;++i)
	{
		cs += cc1[i]*cc2[i];
	}

	det = cc1[0]*(HN_A[0]-HN_1_A[0])+cc1[1]*(HN_A[1]-HN_1_A[1])+cc1[2]*(HN_A[2]-HN_1_A[2]);

	if (det < 0)
	{
		theta = PI - acos(cs);
	}
	else
	{
		theta = PI + acos(cs);
	}

	anglepsi = theta;

	return anglepsi;
/*********************************************************************************/
}
