#include <stdio.h>
#include <math.h>

#include "Vis_MD.h"

void TransLocalToAbso (int nNumAtom1, int nNumAtom2,
                       int nSumNumAtom, double delta_dihed)
{
	int i,alpha;

	int n=0, angle;

	int nNumBod;

	int nNumClut;

	int nNumAtomLoca;

	int nNumAtomAbsoOrig, nNumAtomAbsoTerm;

	double sn2;
	double cs2;

	double Coord[/*MAXNATOM*/10][3];
	double Coord2[/*MAXNATOM*/10][3];
	double O_A[3];

	sn2 = sin(delta_dihed);
	cs2 = cos(delta_dihed);

	nNumAtomAbsoOrig = prot.clt[nNumAtom1].origin-1;
	nNumAtomAbsoTerm = prot.clt[nNumAtom2].term[0]-1;

	for(alpha=0;alpha<3;++alpha)
	{
		O_A[alpha] = prot.atm[nNumAtomAbsoOrig].coord[alpha];
	}

	for(i = 0; i < nSumNumAtom ; ++i)
	{
	  //		Coord[i][0] =   cs2*prot.clt[nNumAtom1].xoord_clust[nNumBod][nNumAtomLoca+i][0]
	  //		              + sn2*prot.clt[nNumAtom1].xoord_clust[nNumBod][nNumAtomLoca+i][1];
	  //
	  //	Coord[i][1] = - sn2*prot.clt[nNumAtom1].xoord_clust[nNumBod][nNumAtomLoca+i][0]
	  //	              + cs2*prot.clt[nNumAtom1].xoord_clust[nNumBod][nNumAtomLoca+i][1];
	  //
	  //	Coord[i][2] = prot.clt[nNumAtom1].xoord_clust[nNumBod][nNumAtomLoca+i][2];
	}

	for(i = 0; i < nSumNumAtom ; ++i)
	{
		Coord2[i][0] =   MatTransAbsoToLocal[0][0]*Coord[i][0]
		               + MatTransAbsoToLocal[1][0]*Coord[i][1]
		               + MatTransAbsoToLocal[2][0]*Coord[i][2];

		Coord2[i][1] =   MatTransAbsoToLocal[0][1]*Coord[i][0]
		               + MatTransAbsoToLocal[1][1]*Coord[i][1]
		               + MatTransAbsoToLocal[2][1]*Coord[i][2];

		Coord2[i][2] =   MatTransAbsoToLocal[0][2]*Coord[i][0]
		               + MatTransAbsoToLocal[1][2]*Coord[i][1]
		               + MatTransAbsoToLocal[2][2]*Coord[i][2];
	}

	for(i = 0; i < nSumNumAtom ; ++i)
	{
		for(alpha=0 ;alpha<3; ++alpha)
		{
			prot.atm[nNumAtomAbsoTerm+i].coord[alpha]
			 = Coord2[i][alpha] + O_A[alpha];
		}
	}

	for(nNumClut=0; nNumClut<prot.nNumClut; ++nNumClut)
	{
/*		TransAbsoToLocal(nNumClut);*/
	}
}
