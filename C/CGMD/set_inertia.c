#include <stdio.h>
#include <math.h>

#include "gener.h"
#include "ABA.h" // #include "ABA_multimer.h"

void PsedouInertiaMatrix(int nNumClut);
//int setJoin(int nNumClut, int joinflag);

// 慣性モーメント、慣性行列の設定を行う関数
void set(void)
{
	int alpha, alpha2;
	int nNumClut;
	int num_atom_clust_total;
	int joinflag;

	ident[0][0]=1.0;
	ident[0][1]=0.0;
	ident[0][2]=0.0;
	ident[1][0]=0.0;
	ident[1][1]=1.0;
	ident[1][2]=0.0;
	ident[2][0]=0.0;
	ident[2][1]=0.0;
	ident[2][2]=1.0;
	ident[0][3]=0.0;
	ident[0][4]=0.0;
	ident[0][5]=0.0;
	ident[1][3]=0.0;
	ident[1][4]=0.0;
	ident[1][5]=0.0;
	ident[2][3]=0.0;
	ident[2][4]=0.0;
	ident[2][5]=0.0;
	ident[3][0]=0.0;
	ident[3][1]=0.0;
	ident[3][2]=0.0;
	ident[3][3]=1.0;
	ident[3][4]=0.0;
	ident[3][5]=0.0;
	ident[4][0]=0.0;
	ident[4][1]=0.0;
	ident[4][2]=0.0;
	ident[4][3]=0.0;
	ident[4][4]=1.0;
	ident[4][5]=0.0;
	ident[5][0]=0.0;
	ident[5][1]=0.0;
	ident[5][2]=0.0;
	ident[5][3]=0.0;
	ident[5][4]=0.0;
	ident[5][5]=1.0;


	// 慣性モーメントの初期化を行う
	for(nNumClut=0; nNumClut < prot.DOF ; ++nNumClut)
	{
		for (alpha=0;alpha<3;++alpha)
		{
			for (alpha2=0;alpha2<3;++alpha2)
			{
				clust[nNumClut].Inertia_clust[alpha][alpha2] = 0.0;
			}
		}
	}

	// 慣性モーメント、慣性行列の設定を行う
	//	joinflag=0;
	for(nNumClut=0; nNumClut<prot.DOF; ++nNumClut)
	{
		num_atom_clust_total = clust[nNumClut].num_xoord_a;
		set_Inertia_clust(nNumClut, num_atom_clust_total);
		InertiaMatrix(nNumClut);
		PsedouInertiaMatrix(nNumClut);
		//		joinflag=setJoin(nNumClut,joinflag);
	}

//	clust[0].InertiaMatrix[0][0]=0.0;
//	clust[0].InertiaMatrix[0][1]=0.0;
//	clust[0].InertiaMatrix[0][2]=0.0;
//	clust[0].InertiaMatrix[0][3]=0.0;
//	clust[0].InertiaMatrix[0][4]=0.0;
//	clust[0].InertiaMatrix[0][5]=0.0;
//	clust[0].InertiaMatrix[1][0]=0.0;
//	clust[0].InertiaMatrix[1][1]=0.0;
//	clust[0].InertiaMatrix[1][2]=0.0;
//	clust[0].InertiaMatrix[1][3]=0.0;
//	clust[0].InertiaMatrix[1][4]=0.0;
//	clust[0].InertiaMatrix[1][5]=0.0;
//	clust[0].InertiaMatrix[2][0]=0.0;
//	clust[0].InertiaMatrix[2][1]=0.0;
//	clust[0].InertiaMatrix[2][2]=0.0;
//	clust[0].InertiaMatrix[2][3]=0.0;
//	clust[0].InertiaMatrix[2][4]=0.0;
//	clust[0].InertiaMatrix[2][5]=0.0;
//	clust[0].InertiaMatrix[3][0]=0.0;
//	clust[0].InertiaMatrix[3][1]=0.0;
//	clust[0].InertiaMatrix[3][2]=0.0;
//	clust[0].InertiaMatrix[3][3]=0.0;
//	clust[0].InertiaMatrix[3][4]=0.0;
//	clust[0].InertiaMatrix[3][5]=0.0;
//	clust[0].InertiaMatrix[4][0]=0.0;
//	clust[0].InertiaMatrix[4][1]=0.0;
//	clust[0].InertiaMatrix[4][2]=0.0;
//	clust[0].InertiaMatrix[4][3]=0.0;
//	clust[0].InertiaMatrix[4][4]=0.0;
//	clust[0].InertiaMatrix[4][5]=0.0;
//	clust[0].InertiaMatrix[5][0]=0.0;
//	clust[0].InertiaMatrix[5][1]=0.0;
//	clust[0].InertiaMatrix[5][2]=0.0;
//	clust[0].InertiaMatrix[5][3]=0.0;
//	clust[0].InertiaMatrix[5][4]=0.0;
//	clust[0].InertiaMatrix[5][5]=0.0;

}

// 剛体の慣性モーメントの設定
void set_Inertia_clust(int nNumClut,int  num_atom_clust_total)
{
	int i,j,k;
	int nNumAtomLoca=0,nNumParent;

	clust[nNumClut].Inertia_clust[0][0]=0.0;
	clust[nNumClut].Inertia_clust[1][1]=0.0;
	clust[nNumClut].Inertia_clust[2][2]=0.0;
	clust[nNumClut].Inertia_clust[0][1]=0.0;
	clust[nNumClut].Inertia_clust[0][2]=0.0;
	clust[nNumClut].Inertia_clust[1][2]=0.0;

	nNumParent = clust[nNumClut].nNumClutOfParent-1;

	//	nNumAtomLoca =  clust[nNumParent].num_atom_clust
	//		   -clust[nNumParent].terminal_atom_a[0]
	//		   +clust[nNumParent].origin_atom_a;
	//
	//	for (i=1;i<nNumClut-nNumParent;++i)
	//{
	//	nNumAtomLoca += clust[nNumParent+i].num_atom_clust;
	//}
	nNumAtomLoca = 0;

	if (clust[nNumClut].num_atom_clust > 1)
	{
	for(i=0; i<=clust[nNumClut].num_atom_clust; ++i)
	{
		clust[nNumClut].Inertia_clust[0][0]/*kg*m^2*/ += 
					clust[nNumClut].mass_clust[/*nNumAtomLoca+*/i]/*u*//**1.660539e-27*/
					*  (    clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][1]/*A*//**1.0e-10*/
					       *clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][1]/*A*//**1.0e-10*/
					     +  clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][2]/*A*//**1.0e-10*/
					       *clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][2]/*A*//**1.0e-10*/);

		clust[nNumClut].Inertia_clust[1][1]/*kg*m^2*/ += 
					clust[nNumClut].mass_clust[/*nNumAtomLoca+*/i]/*u*//**1.660539e-27*/
					*   (   clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][0]/*A*//**1.0e-10*/
					       *clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][0]/*A*//**1.0e-10*/
					      + clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][2]/*A*//**1.0e-10*/
					       *clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][2]/*A*//**1.0e-10*/);

		clust[nNumClut].Inertia_clust[2][2]/*kg*m^2*/ += 
					clust[nNumClut].mass_clust[/*nNumAtomLoca+*/i]/*u*//**1.660539e-27*/
					*(      clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][0]/*A*//**1.0e-10*/
					       *clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][0]/*A*//**1.0e-10*/
				         +  clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][1]/*A*//**1.0e-10*/
					       *clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][1]/*A*//**1.0e-10*/ );

//		if (clust[nNumClut].Inertia_clust[2][2]<1.0e-70)
//		{
//			clust[nNumClut].Inertia_clust[2][2]=1.0e-47;
//		}
		clust[nNumClut].Inertia_clust[0][1]/*kg*m^2*/ -= 
						clust[nNumClut].mass_clust[/*nNumAtomLoca+*/i]/*u*//**1.660539e-27*/
					      * clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][0]/*A*//**1.0e-10*/
					      * clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][1]/*A*//**1.0e-10*/;

		clust[nNumClut].Inertia_clust[0][2]/*kg*m^2*/ -= 
						clust[nNumClut].mass_clust[/*nNumAtomLoca+*/i]/*u*//**1.660539e-27*/
					      * clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][0]/*A*//**1.0e-10*/ 
					      * clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][2]/*A*//**1.0e-10*/;

		clust[nNumClut].Inertia_clust[1][2]/*kg*m^2*/ -= 
					        clust[nNumClut].mass_clust[/*nNumAtomLoca+*/i]/*u*//**1.660539e-27*/
					      * clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][1]/*A*//**1.0e-10*/
					      * clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][2]/*A*//**1.0e-10*/;
	}

	clust[nNumClut].Inertia_clust[1][0]/*kg*m^2*/=clust[nNumClut].Inertia_clust[0][1]/*kg*m^2*/;
	clust[nNumClut].Inertia_clust[2][0]/*kg*m^2*/=clust[nNumClut].Inertia_clust[0][2]/*kg*m^2*/;
	clust[nNumClut].Inertia_clust[2][1]/*kg*m^2*/=clust[nNumClut].Inertia_clust[1][2]/*kg*m^2*/;

	j = 0;
	k = 0;
	for(i=0; i<num_atom_clust_total; ++i)
	{
		if (k < clust[nNumClut+j].num_atom_clust)
		{
			k = 0;
			++j;
		}
		clust[nNumClut].Inertia_clust_total/*kg*m^2*/ += 
					clust[nNumClut+j].mass_clust[/*nNumAtomLoca+*/k]/*u*//**1.660539e-27*/
					*(      clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][0]/*A*//**1.0e-10*/
					       *clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][0]/*A*//**1.0e-10*/
				         +  clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][1]/*A*//**1.0e-10*/
					       *clust[nNumClut].xoord_clust/*[0]*/[nNumAtomLoca+i][1]/*A*//**1.0e-10*/ );
		++k;
	}

	}
	else
	{
	  clust[nNumClut].Inertia_clust[0][0] = 
	    clust[nNumClut].mass_clust[0]
	    *(1.0*1.0+1.0*1.0);
	  
	  clust[nNumClut].Inertia_clust[1][1] = 
	    clust[nNumClut].mass_clust[0]
	    *(1.0*1.0+1.0*1.0);

	  clust[nNumClut].Inertia_clust[2][2] = 
	    clust[nNumClut].mass_clust[0]
	    *(1.0*1.0+1.0*1.0);
	}
}

// 剛体内原子の質量和の設定
double Sum_Mass(int nNumClut)
{
	int i;
	double sum_mass_clut=0.0;

	for(i=0; i<clust[nNumClut].num_atom_clust; ++i)
	{
		sum_mass_clut/*u*/+=clust[nNumClut].mass_clust[i]/*u*/;
	}

	clust[nNumClut].sum_mass = sum_mass_clut;

	return sum_mass_clut/*u*/;
}

// 剛体の慣性行列の設定
void InertiaMatrix(int nNumClut)
{
	int nNumAtom;
	int alpha, alpha2;
	int nNumAtomLoca=0,nNumParent;
	int nNumAtomLoca2;

	for(alpha=0;alpha<3;++alpha)
	{
		for(alpha2=0;alpha2<3;++alpha2)
		{
			clust[nNumClut].InertiaMatrix[alpha][alpha2]/*kg*m^2*/ 
						= /*2.0**/clust[nNumClut].Inertia_clust[alpha][alpha2]/*kg*m^2*/;
//						= 2.0*clust[nNumClut].Inertia_clust[alpha][alpha2]/*kg*m^2*/;
		}
	}

//	for(alpha=0;alpha<3;++alpha)
//	{
//		for(alpha2=0;alpha2<3;++alpha2)
//		{
//			clust[nNumClut].InertiaMatrix[alpha][alpha2]/*kg*m^2*/ 
//						= 0.0;
//		}
//	}

	clust[nNumClut].InertiaMatrix[3][3]/*kg*/=Sum_Mass(nNumClut)/*u*//**1.660539e-27*/;
	clust[nNumClut].InertiaMatrix[3][4]=0.0;
	clust[nNumClut].InertiaMatrix[3][5]=0.0;
	clust[nNumClut].InertiaMatrix[4][3]=0.0;
	clust[nNumClut].InertiaMatrix[4][4]/*kg*/=clust[nNumClut].InertiaMatrix[3][3]/*kg*/;
	clust[nNumClut].InertiaMatrix[4][5]=0.0;
	clust[nNumClut].InertiaMatrix[5][3]=0.0;
	clust[nNumClut].InertiaMatrix[5][4]=0.0;
	clust[nNumClut].InertiaMatrix[5][5]/*kg*/=clust[nNumClut].InertiaMatrix[3][3]/*kg*/;

//	for(i=0; i<clust[nNumClut].num_atom_clust; ++i)
//	{
//		for(j=0;j<3;++j)
//		{
//			mq[j]/*kg*m*/ +=  clust[nNumClut].mass_clust[i]/*u*/*1.660539e-27
//					  *(clust[nNumClut].xoord_clust[i][j]
//						-clust[nNumClut].xoord_clust/*[0]*/[j])/*A*/
//					  *1.0e-10;
//		}
//	}

	nNumParent = clust[nNumClut].nNumClutOfParent-1;

       	nNumAtomLoca =  0;/*04_10*///clust[nNumParent].num_atom_clust
				  // -clust[nNumParent].terminal_atom_a[0]
				  // +clust[nNumParent].origin_atom_a;


	//	nNumAtomLoca2 =  clust[nNumClut].origin_atom_a-1;


	for (alpha=0;alpha<3;++alpha)
	{
			mq[alpha] = 0.0;
			clust[nNumClut].qCOM[alpha] = 0.0;
	}

	for(nNumAtom=0; nNumAtom<clust[nNumClut].num_atom_clust; ++nNumAtom)
	{
		mq[0]/*kg*m*/ +=  clust[nNumClut].mass_clust[nNumAtom]/*u*//**1.660539e-27*/
					     *( clust[nNumClut].xoord_clust/*[0]*/[nNumAtom+nNumAtomLoca][0]
						/*-clust[nNumClut].xoord_clust[0][nNumAtomLoca][0]*/)/*A*/
					      /**1.0e-10*/;
//		mq[0]/*kg*m*/ +=  clust[nNumClut].mass_clust[nNumAtom]/*u*//**1.660539e-27*/
//					     *( prot.coord[nNumAtom+nNumAtomLoca2][0]
//						   -prot.coord[nNumAtomLoca2][0])/*A*/
//					      /**1.0e-10*/;
//		clust[nNumClut].qCOM[0] += mq[0]/Sum_Mass(nNumClut);
	}
		clust[nNumClut].qCOM[0] /*+*/= mq[0]/Sum_Mass(nNumClut);

	for(nNumAtom=0; nNumAtom<clust[nNumClut].num_atom_clust; ++nNumAtom)
	{
		mq[1]/*kg*m*/ +=  clust[nNumClut].mass_clust[nNumAtom]/*u*//**1.660539e-27*/
					     *( clust[nNumClut].xoord_clust/*[0]*/[nNumAtom+nNumAtomLoca][1]
						/*-clust[nNumClut].xoord_clust[0][nNumAtomLoca][1]*/)/*A*/
					      /**1.0e-10*/;
//		mq[1]/*kg*m*/ +=  clust[nNumClut].mass_clust[nNumAtom]/*u*//**1.660539e-27*/
//					     *( prot.coord[nNumAtom+nNumAtomLoca2][1]
//						   -prot.coord[nNumAtomLoca2][1])/*A*/
//					      /**1.0e-10*/;
//		clust[nNumClut].qCOM[1] += mq[1]/Sum_Mass(nNumClut);
	}
		clust[nNumClut].qCOM[1] /*+*/= mq[1]/Sum_Mass(nNumClut);

	for(nNumAtom=0; nNumAtom<clust[nNumClut].num_atom_clust; ++nNumAtom)
	{
		mq[2]/*kg*m*/ +=  clust[nNumClut].mass_clust[nNumAtom]/*u*//**1.660539e-27*/
					     *( clust[nNumClut].xoord_clust/*[0]*/[nNumAtom+nNumAtomLoca][2]
						/*-clust[nNumClut].xoord_clust[0][nNumAtomLoca][2]*/)/*A*/
					      /**1.0e-10*/;
//		mq[2]/*kg*m*/ +=  clust[nNumClut].mass_clust[nNumAtom]/*u*//**1.660539e-27*/
//					     *( prot.coord[nNumAtom+nNumAtomLoca2][2]
//						   -prot.coord[nNumAtomLoca2][2])/*A*/
//					      /**1.0e-10*/;
//		clust[nNumClut].qCOM[2] += mq[2]/Sum_Mass(nNumClut);
	}
		clust[nNumClut].qCOM[2] /*+*/= mq[2]/Sum_Mass(nNumClut);

	if (clust[nNumClut].num_atom_clust == 1)
	{
		for (alpha=0;alpha<3;++alpha)
		{
			mq[alpha] = 0.1*Sum_Mass(nNumClut);
			clust[nNumClut].qCOM[alpha] = mq[alpha]/Sum_Mass(nNumClut);
		}
	}

	clust[nNumClut].InertiaMatrix[0][3]/*kg*m*/=0.0;
	clust[nNumClut].InertiaMatrix[0][4]/*kg*m*/=-mq[2]/*kg*m*/;
	clust[nNumClut].InertiaMatrix[0][5]/*kg*m*/=mq[1]/*kg*m*/;
	clust[nNumClut].InertiaMatrix[1][3]/*kg*m*/=mq[2]/*kg*m*/;
	clust[nNumClut].InertiaMatrix[1][4]/*kg*m*/=0.0;
	clust[nNumClut].InertiaMatrix[1][5]/*kg*m*/=-mq[0]/*kg*m*/;
	clust[nNumClut].InertiaMatrix[2][3]/*kg*m*/=-mq[1]/*kg*m*/;
	clust[nNumClut].InertiaMatrix[2][4]/*kg*m*/=mq[0]/*kg*m*/;
	clust[nNumClut].InertiaMatrix[2][5]/*kg*m*/=0.0;

	clust[nNumClut].InertiaMatrix[3][0]/*kg*m*/=0.0;
	clust[nNumClut].InertiaMatrix[3][1]/*kg*m*/=mq[2]/*kg*m*/;
	clust[nNumClut].InertiaMatrix[3][2]/*kg*m*/=-mq[1]/*kg*m*/;
	clust[nNumClut].InertiaMatrix[4][0]/*kg*m*/=-mq[2]/*kg*m*/;
	clust[nNumClut].InertiaMatrix[4][1]/*kg*m*/=0.0;
	clust[nNumClut].InertiaMatrix[4][2]/*kg*m*/=mq[0]/*kg*m*/;
	clust[nNumClut].InertiaMatrix[5][0]/*kg*m*/=mq[1]/*kg*m*/;
	clust[nNumClut].InertiaMatrix[5][1]/*kg*m*/=-mq[0]/*kg*m*/;
	clust[nNumClut].InertiaMatrix[5][2]/*kg*m*/=0.0;

//	clust[nNumClut].InertiaMatrix[0][3]/*kg*m*/=0.0;
//	clust[nNumClut].InertiaMatrix[0][4]/*kg*m*/=mq[2]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[0][5]/*kg*m*/=-mq[1]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[1][3]/*kg*m*/=-mq[2]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[1][4]/*kg*m*/=0.0;
//	clust[nNumClut].InertiaMatrix[1][5]/*kg*m*/=mq[0]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[2][3]/*kg*m*/=mq[1]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[2][4]/*kg*m*/=-mq[0]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[2][5]/*kg*m*/=0.0;
//
//	clust[nNumClut].InertiaMatrix[3][0]/*kg*m*/=0.0;
//	clust[nNumClut].InertiaMatrix[3][1]/*kg*m*/=-mq[2]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[3][2]/*kg*m*/=mq[1]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[4][0]/*kg*m*/=mq[2]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[4][1]/*kg*m*/=0.0;
//	clust[nNumClut].InertiaMatrix[4][2]/*kg*m*/=-mq[0]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[5][0]/*kg*m*/=-mq[1]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[5][1]/*kg*m*/=mq[0]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[5][2]/*kg*m*/=0.0;

//	clust[nNumClut].InertiaMatrix[0][3]/*kg*m*/=0.0;
//	clust[nNumClut].InertiaMatrix[0][4]/*kg*m*/=-mq[2]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[0][5]/*kg*m*/=mq[1]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[1][3]/*kg*m*/=mq[2]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[1][4]/*kg*m*/=0.0;
//	clust[nNumClut].InertiaMatrix[1][5]/*kg*m*/=-mq[0]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[2][3]/*kg*m*/=-mq[1]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[2][4]/*kg*m*/=mq[0]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[2][5]/*kg*m*/=0.0;
//
//	clust[nNumClut].InertiaMatrix[3][0]/*kg*m*/=0.0;
//	clust[nNumClut].InertiaMatrix[3][1]/*kg*m*/=-mq[2]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[3][2]/*kg*m*/=mq[1]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[4][0]/*kg*m*/=mq[2]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[4][1]/*kg*m*/=0.0;
//	clust[nNumClut].InertiaMatrix[4][2]/*kg*m*/=-mq[0]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[5][0]/*kg*m*/=-mq[1]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[5][1]/*kg*m*/=mq[0]/*kg*m*/;
//	clust[nNumClut].InertiaMatrix[5][2]/*kg*m*/=0.0;

//	clust[nNumClut].InertiaMatrix[0][0]/*kg*m^2*/ 
//					=( clust[nNumClut].qCOM[1]*clust[nNumClut].qCOM[1]/*kg*m^2*/
//					  +clust[nNumClut].qCOM[2]*clust[nNumClut].qCOM[2])*Sum_Mass(nNumClut)/*kg*m^2*/;
//	clust[nNumClut].InertiaMatrix[1][1]/*kg*m^2*/ 
//					=( clust[nNumClut].qCOM[0]*clust[nNumClut].qCOM[0]/*kg*m^2*/
//					  +clust[nNumClut].qCOM[2]*clust[nNumClut].qCOM[2])*Sum_Mass(nNumClut)/*kg*m^2*/;
//	clust[nNumClut].InertiaMatrix[2][2]/*kg*m^2*/ 
//					=( clust[nNumClut].qCOM[0]*clust[nNumClut].qCOM[0]/*kg*m^2*/
//					  +clust[nNumClut].qCOM[1]*clust[nNumClut].qCOM[1])*Sum_Mass(nNumClut)/*kg*m^2*/;
//	clust[nNumClut].InertiaMatrix[0][1]/*kg*m^2*/ 
//					=  -clust[nNumClut].qCOM[0]*clust[nNumClut].qCOM[1]*Sum_Mass(nNumClut)/*kg*m^2*/;
//	clust[nNumClut].InertiaMatrix[1][0]/*kg*m^2*/ 
//					=  -clust[nNumClut].qCOM[1]*clust[nNumClut].qCOM[0]*Sum_Mass(nNumClut)/*kg*m^2*/;
//	clust[nNumClut].InertiaMatrix[0][2]/*kg*m^2*/ 
//					=  -clust[nNumClut].qCOM[0]*clust[nNumClut].qCOM[2]*Sum_Mass(nNumClut)/*kg*m^2*/;
//	clust[nNumClut].InertiaMatrix[2][0]/*kg*m^2*/ 
//					=  -clust[nNumClut].qCOM[2]*clust[nNumClut].qCOM[0]*Sum_Mass(nNumClut)/*kg*m^2*/;
//	clust[nNumClut].InertiaMatrix[1][2]/*kg*m^2*/ 
//					=  -clust[nNumClut].qCOM[1]*clust[nNumClut].qCOM[2]*Sum_Mass(nNumClut)/*kg*m^2*/;
//	clust[nNumClut].InertiaMatrix[2][1]/*kg*m^2*/ 
//					=  -clust[nNumClut].qCOM[2]*clust[nNumClut].qCOM[1]*Sum_Mass(nNumClut)/*kg*m^2*/;

//	clust[nNumClut].InertiaMatrix[0][0]= 2*clust[nNumClut].Inertia_clust[0][0];
//	clust[nNumClut].InertiaMatrix[1][1]= 2*clust[nNumClut].Inertia_clust[1][1];
//	clust[nNumClut].InertiaMatrix[2][2]= 2*clust[nNumClut].Inertia_clust[2][2];
//	clust[nNumClut].InertiaMatrix[0][1]= -2*clust[nNumClut].Inertia_clust[0][1];
//	clust[nNumClut].InertiaMatrix[1][0]= -2*clust[nNumClut].Inertia_clust[1][0];
//	clust[nNumClut].InertiaMatrix[0][2]= -2*clust[nNumClut].Inertia_clust[0][2];
//	clust[nNumClut].InertiaMatrix[2][0]= -2*clust[nNumClut].Inertia_clust[2][0];
//	clust[nNumClut].InertiaMatrix[1][2]= -2*clust[nNumClut].Inertia_clust[1][2];
//	clust[nNumClut].InertiaMatrix[2][1]= -2*clust[nNumClut].Inertia_clust[1][2];

//	clust[nNumClut].InertiaMatrix[0][1]= -clust[nNumClut].Inertia_clust[0][1];
//	clust[nNumClut].InertiaMatrix[1][0]= -clust[nNumClut].Inertia_clust[1][0];
//	clust[nNumClut].InertiaMatrix[0][2]= -clust[nNumClut].Inertia_clust[0][2];
//	clust[nNumClut].InertiaMatrix[2][0]= -clust[nNumClut].Inertia_clust[2][0];
//	clust[nNumClut].InertiaMatrix[1][2]= -clust[nNumClut].Inertia_clust[1][2];
//	clust[nNumClut].InertiaMatrix[2][1]= -clust[nNumClut].Inertia_clust[2][1];

/*	if (nNumClut==5){
	clust[nNumClut].InertiaMatrix[0][0]=	5.4510446423439847;
	clust[nNumClut].InertiaMatrix[0][1]=	-7.6269643086891174e-21; 
	clust[nNumClut].InertiaMatrix[0][2]=	9.2540328312366107; 
	clust[nNumClut].InertiaMatrix[0][3]=	0;
	clust[nNumClut].InertiaMatrix[0][4]=    4.8885238689512001; 
	clust[nNumClut].InertiaMatrix[0][5]=	12.595775803491005; 
	clust[nNumClut].InertiaMatrix[1][0]=   -7.6269643086891174e-21;
	clust[nNumClut].InertiaMatrix[1][1]=   21.161267041986736;
	clust[nNumClut].InertiaMatrix[1][2]=    4.4926275592943815e-21;
	clust[nNumClut].InertiaMatrix[1][3]=    -4.8885238689512001;
	clust[nNumClut].InertiaMatrix[1][4]=    0;
	clust[nNumClut].InertiaMatrix[1][5]=    -7.1168233407018455; 
	clust[nNumClut].InertiaMatrix[2][0]=    9.2540328312366107;
	clust[nNumClut].InertiaMatrix[2][1]=    4.4926275592943815e-21;
        clust[nNumClut].InertiaMatrix[2][2]=   15.710222399642751;
	clust[nNumClut].InertiaMatrix[2][3]=    -12.595775803491005;
	clust[nNumClut].InertiaMatrix[2][4]=         7.1168233407018455;
	clust[nNumClut].InertiaMatrix[2][5]=         0;
	clust[nNumClut].InertiaMatrix[3][0]=         0 ;
	clust[nNumClut].InertiaMatrix[3][1]=         -4.8885238689512001;
	clust[nNumClut].InertiaMatrix[3][2]=         -12.595775803491005;
	clust[nNumClut].InertiaMatrix[3][3]=         16.02;
	clust[nNumClut].InertiaMatrix[3][4]=        0;
	clust[nNumClut].InertiaMatrix[3][5]=        0;
	clust[nNumClut].InertiaMatrix[4][0]=       4.8885238689512001;
	clust[nNumClut].InertiaMatrix[4][1]=        0;
	clust[nNumClut].InertiaMatrix[4][2]=         7.1168233407018455;
	clust[nNumClut].InertiaMatrix[4][3]=       0;
	clust[nNumClut].InertiaMatrix[4][4]=	 16.02;
	clust[nNumClut].InertiaMatrix[4][5]=       0 ;
	clust[nNumClut].InertiaMatrix[5][0]=         12.595775803491005;
	clust[nNumClut].InertiaMatrix[5][1]= -7.1168233407018455;
	clust[nNumClut].InertiaMatrix[5][2]= 0;
	clust[nNumClut].InertiaMatrix[5][3]= 0;
	clust[nNumClut].InertiaMatrix[5][4]= 0;
	clust[nNumClut].InertiaMatrix[5][5]= 16.02;
	clust[nNumClut].qCOM[0] = 0.44424615110498411;
	clust[nNumClut].qCOM[1] = 0.7862531712541202;
	clust[nNumClut].qCOM[2] = -0.30515130268109864;
	}*/

}

// 剛体の擬似慣性行列の設定
void PsedouInertiaMatrix(int nNumClut)
{
	int i;

	clust[nNumClut].PsedoInertia[0][0] =0.5*( -clust[nNumClut].InertiaMatrix[0][0]
										      +clust[nNumClut].InertiaMatrix[1][1]
										      +clust[nNumClut].InertiaMatrix[2][2] );
	clust[nNumClut].PsedoInertia[1][1] =0.5*(  clust[nNumClut].InertiaMatrix[0][0]
										      -clust[nNumClut].InertiaMatrix[1][1]
										      +clust[nNumClut].InertiaMatrix[2][2] );
	clust[nNumClut].PsedoInertia[2][2] =0.5*(  clust[nNumClut].InertiaMatrix[0][0]
										      +clust[nNumClut].InertiaMatrix[1][1]
										      -clust[nNumClut].InertiaMatrix[2][2] );

	clust[nNumClut].PsedoInertia[3][3] =clust[nNumClut].sum_mass;

	clust[nNumClut].PsedoInertia[0][1] = clust[nNumClut].InertiaMatrix[0][1];
	clust[nNumClut].PsedoInertia[0][2] = clust[nNumClut].InertiaMatrix[0][2];
	clust[nNumClut].PsedoInertia[1][2] = clust[nNumClut].InertiaMatrix[1][2];
	clust[nNumClut].PsedoInertia[1][0] = clust[nNumClut].InertiaMatrix[0][1];
	clust[nNumClut].PsedoInertia[2][0] = clust[nNumClut].InertiaMatrix[0][2];
	clust[nNumClut].PsedoInertia[2][1] = clust[nNumClut].InertiaMatrix[2][1];

	clust[nNumClut].PsedoInertia[0][3] = clust[nNumClut].InertiaMatrix[4][2];
	clust[nNumClut].PsedoInertia[1][3] = clust[nNumClut].InertiaMatrix[5][0];
	clust[nNumClut].PsedoInertia[2][3] = clust[nNumClut].InertiaMatrix[3][1];
	clust[nNumClut].PsedoInertia[3][0] = clust[nNumClut].InertiaMatrix[4][2];
	clust[nNumClut].PsedoInertia[3][1] = clust[nNumClut].InertiaMatrix[5][0];
	clust[nNumClut].PsedoInertia[3][2] = clust[nNumClut].InertiaMatrix[3][1];
}

int setJoin(int nNumClut, int joinflag) {

  clust[nNumClut].join = joinflag;
  
  if (clust[nNumClut].num_branch > 1) {
    joinflag +=1;
  }
  if (clust[nNumClut].terminal == 0) {
    joinflag -=1;
  }

  return joinflag;
}
