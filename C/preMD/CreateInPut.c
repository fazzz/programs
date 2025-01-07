#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Vis_MD.h"

#define ConstScaleOfCharge 18.2223

void CreatePDBFile(FILE *input);
int CheckDihedPair(int atoms1[4], int atoms2[4]);
void ascatomnum(int atomnum[20],int num);
void CreateIsouFile(FILE* input);
void CreatemolFile(FILE *input);

// �C���v�b�g�̍쐬���s���֐�
void CreateInPut(int TotalAtomType)
{
	FILE *inputfile;

	// ���W�C���v�b�g�̍쐬
	///////////////////////////////////////////////////////////////////////////
	/***********************************************************************************/
        /* if ((inputfile = fopen(/\*InpfilCOORD/\**\/"crd.in"/\**\/,"w")) == NULL)	   */
	/* {										   */
	/* 	printf("error crd.in cannot open\n");					   */
	/* 	exit(1);								   */
	/* }										   */
	/* 										   */
	/* CreateCoordFile(inputfile);							   */
	/* 										   */
	/* fclose(inputfile);								   */
        /***********************************************************************************/
	///////////////////////////////////////////////////////////////////////////

	// ���W�C���v�b�g�̍쐬
	///////////////////////////////////////////////////////////////////////////
	/**************************************************************/
        /* if ((inputfile = fopen("protein.pdb","w")) == NULL)	      */
	/* {							      */
	/* 	printf("error crd.pdb cannot open\n");		      */
	/* 	exit(1);					      */
	/* }							      */
	/* 							      */
	/* CreatePDBFile(inputfile);				      */
	/* 							      */
	/* fclose(inputfile);					      */
        /**************************************************************/
	///////////////////////////////////////////////////////////////////////////

	// �c����C���v�b�g�̍쐬
	///////////////////////////////////////////////////////////////////////////
	/************************************************************************/
        /* if ((inputfile = fopen(/\*InpfilSEQ*\/"seq.in","w")) == NULL)        */
	/* {								        */
	/* 	printf("error seq.in cannot open\n");			        */
	/* 	exit(1);						        */
	/* }								        */
	/* 								        */
	/* CreateSeqFile(inputfile);					        */
	/* 								        */
	/* fclose(inputfile);						        */
        /************************************************************************/
	///////////////////////////////////////////////////////////////////////////

	// ���̏��C���v�b�g�̍쐬
	///////////////////////////////////////////////////////////////////////////
	//	InpfilCLUST="clust.in";
	if ((inputfile = fopen(InpfilCLUST,"w")) == NULL)
	{
	  		printf("error clust.in cannot open\n");
			exit(1);
	}

	CreateClustFile(inputfile);

	fclose(inputfile);
	///////////////////////////////////////////////////////////////////////////

	// �͏���C���v�b�g�̍쐬
	///////////////////////////////////////////////////////////////////////////
	/*********************************************************************************/
        /* if ((inputfile = fopen(/\*InpfilTOP/\**\/"top.in"/\**\/,"w")) == NULL)	 */
	/* {										 */
	/* 	printf("error top.in cannot open\n");					 */
	/* 	exit(1);								 */
	/* }										 */
	/* 										 */
	/* CreateTopFile(inputfile, TotalAtomType);					 */
	/* 										 */
	/* fclose(inputfile);								 */
        /*********************************************************************************/
	///////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////
	/***********************************************************/
        /* if ((inputfile = fopen("isou.dat","w")) == NULL)	   */
	/* {							   */
	/* 	printf("error isou.dat cannot open\n");		   */
	/* 	exit(1);					   */
	/* }							   */
	/* 							   */
	/* CreateIsouFile(inputfile);				   */
	/* 							   */
	/* fclose(inputfile);					   */
        /***********************************************************/
	///////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////
	/**********************************************************/
        /* if ((inputfile = fopen("mol.dat","w")) == NULL)	  */
	/* {							  */
	/* 	printf("error mol.dat cannot open\n");		  */
	/* 	exit(1);					  */
	/* }							  */
	/* 							  */
	/* CreatemolFile(inputfile);				  */
	/* 							  */
	/* fclose(inputfile);					  */
        /**********************************************************/
	///////////////////////////////////////////////////////////////////////////

}

// ���W�C���v�b�g�̍쐬���s���֐�
void CreateCoordFile(FILE *input)
{
	int nNumAtom;
	int alpha;

	// ���q���̋L�q
	fprintf(input, "%s\n", prot.Sequence[0]);

	// ���W�f�[�^�̋L�q
	for(nNumAtom=0; nNumAtom<prot.nNumAtom; ++nNumAtom)
	{
		for (alpha=0;alpha<3;++alpha)
		{
			fprintf(input, "%8.3lf  ", prot.atm[nNumAtom].coord[alpha]+0.1);
		}
		fprintf(input, "\n");
	}
}

// ���W�C���v�b�g�̍쐬���s���֐�
void CreatePDBFile(FILE *input)
{
	int nNumAtom=0,nNumAtom2, nNumRes;
	int nNumAtomdummy=0;
	int alpha;
	char *name;

	// ���q���̋L�q
//	fprintf(input, "%s\n", *prot.name);

	// ���W�f�[�^�̋L�q
//	for(nNumAtom=0; nNumAtom<prot.nNumAtom; ++nNumAtom)
//	{
	for (nNumRes=0; nNumRes<prot.nNumResidue; ++nNumRes)
	{
		for (nNumAtom2=0;nNumAtom2<prot.numsequence[nNumRes]-1;++nNumAtom2)
		{
		switch (prot.atm[nNumAtom].atomID.atomnum ) {
			case 1: name="N";
			case 2: name="H";
			case 3: name="C";
			case 4: name="H";
			case 5: name="C";
			case 6: name="O";
			default: name="C";
		}
			fprintf(input, "ATOM%7d  %-4s%-5s %4d%11.3lf%8.3lf%8.3lf\n",
					nNumAtom+1,
					name,
					prot.Sequence[nNumRes],
					nNumRes+1,
					prot.atm[nNumAtom].coord[0],prot.atm[nNumAtom].coord[1],prot.atm[nNumAtom].coord[2]);
			++nNumAtom;
		}
	}
}

// �c����C���v�b�g�̍쐬���s���֐�
void CreateSeqFile(FILE *input)
{
	int nNumRes;

	// �c��f�[�^�̋L�q
	for (nNumRes=0; nNumRes<prot.nNumResidue; ++nNumRes)
	{
		fprintf(input, "%s ", prot.Sequence[nNumRes]);
	}
	fprintf(input, "\n");
}

// ���̏��C���v�b�g�̍쐬���s���֐�
void CreateClustFile(FILE *input)
{
	int i;
	int nNumClust;
	int nNumBranch;

	// �^���p�N�����̍��̐��̋L�q
	fprintf(input, "%2d\n", prot.nNumClut);

	// ���̂̌��_���̋L�q
	for (nNumClust=0; nNumClust<prot.nNumClut; ++nNumClust)
	{
		fprintf(input, "%5d ", prot.clt[nNumClust].origin);
		if ((nNumClust+1)%10==0)
		fprintf(input, "\n ");
	}
	fprintf(input, "\n");

	// ���̂��I�_���ǂ���
	for (nNumClust=0; nNumClust<prot.nNumClut; ++nNumClust)
	{
		fprintf(input, "%5d ", prot.clt[nNumClust].termflag);
		if ((nNumClust+1)%10==0)
		fprintf(input, "\n ");
	}
	fprintf(input, "\n");

	// ���̒��̌��q���̋L�q
	for (nNumClust=0; nNumClust<prot.nNumClut; ++nNumClust)
	{
		fprintf(input, "%5d ", prot.clt[nNumClust].nNumAtom);
		if ((nNumClust+1)%10==0)
		fprintf(input, "\n ");
	}
	fprintf(input, "\n");

	// ���̂̎}�̐��̋L�q
	for (nNumClust=0; nNumClust<prot.nNumClut-1; ++nNumClust)
	{
		fprintf(input, "%5d ", prot.clt[nNumClust].nNumBranch);
		if ((nNumClust+1)%10==0)
		fprintf(input, "\n ");
	}
	fprintf(input, "1 ");
	fprintf(input, "\n");

	// ���̂̎��R�x�̋L�q
	for (nNumClust=0; nNumClust<prot.nNumClut; ++nNumClust)
	{
		fprintf(input, "%5d ", prot.clt[nNumClust].nDegOfFoo);
		if ((nNumClust+1)%10==0)
		fprintf(input, "\n ");
	}
	fprintf(input, "\n");

	// ���̂̏I�_���̋L�q
	for (nNumClust=0; nNumClust<prot.nNumClut-1; ++nNumClust)
	{
		for (  nNumBranch=0;
			   nNumBranch<prot.clt[nNumClust].nNumBranch;
		     ++nNumBranch)
		{
			fprintf(input, "%5d ", prot.clt[nNumClust].term[nNumBranch]);
		}
		if ((nNumClust+1)%10==0)
		fprintf(input, "\n ");
	}
	fprintf(input, "%5d ", prot.clt[nNumClust].term[0]);
	fprintf(input, "\n");

	// ���̂�"�e"�̏��̋L�q
	for (nNumClust=0; nNumClust<prot.nNumClut; ++nNumClust)
	{
		fprintf(input, "%5d ", prot.clt[nNumClust].nNumParentClust);
		if ((nNumClust+1)%10==0)
		fprintf(input, "\n ");
	}
	fprintf(input, "\n");

	// ���̂�"�q"�̏��̋L�q
	for (nNumClust=0; nNumClust<prot.nNumClut-1; ++nNumClust)
	{
		if (prot.clt[nNumClust].termflag==0)
		{
		  fprintf(input, "%5d ",-1);
		}
		else
		{
			for (  nNumBranch=0;
				   nNumBranch<prot.clt[nNumClust].nNumBranch;
			     ++nNumBranch)
			{
				fprintf(input, "%5d ", prot.clt[nNumClust].nNumChildClust[nNumBranch]);
			}
		}
		if ((nNumClust+1)%10==0)
		fprintf(input, "\n ");
	}
	fprintf(input, "%5d ",-1);
	fprintf(input, "\n");

	i = 1;
	// "�v�Z��"�̏��̋L�q
	for (nNumClust=0; nNumClust<prot.nNumClut; ++nNumClust)
	{
		fprintf(input, "%5d ", i);
		++i;
		if ((nNumClust+1)%10==0)
		fprintf(input, "\n ");
	}
	fprintf(input, "\n");

	// Ramachandren MAP �쐬�� flag
	fprintf(input, "%5d ", 0);
	fprintf(input, "\n");
}

// �͏���C���v�b�g�̍쐬���s���֐�
void CreateTopFile(FILE *input, int TotalAtomType)
{
	int i;
	int nNumAtom;
	int nNumAtom2;
	int nNumAtom3;
	int nNumClut;
	int nNumDihed;
	int TypeTotal;
	int mat[20/*10*/][20/*10*/];
	int CheckDihed;
	int numdihedsofar;

	double depth;
	double length;

///////////////////////////////////////////////////////////////////////////////
	// ��ʊp�f�[�^�̋L�q
		// ����ʊp���̋L�q
//		fprintf(input, "%d \n", prot.nNumDihedAng/2);
		fprintf(input, "0 \n");
		// ��ʊp�퐔�̋L�q
//		fprintf(input, "%d \n", prot.nNumDihedType);
		fprintf(input, "0 \n");

		// �e���̂̍\�����q�ԍ��̋L�q
		for (nNumDihed=0; nNumDihed<prot.nNumDihedAng; ++nNumDihed)
		{
			CheckDihed = 0;
			for(numdihedsofar=0;numdihedsofar<prot.nNumDihedAng;++numdihedsofar)
			{
				if(CheckDihedPair(prot.dihedang[numdihedsofar].atomnum,prot.dihedang[nNumDihed].atomnum) == 1)
				{
					CheckDihed = 1;
					break;
				}
			}
			if(CheckDihed == 0)
			{
			for (i=0;i<4;++i)
			{
//				fprintf(input, "%8d ", prot.dihedang[nNumDihed].atomnum[i]*3);
			}
//			fprintf(input, "%8d ", prot.dihedang[nNumDihed].atomnum[4]);
//			fprintf(input, "%8d ", prot.dihedang[nNumDihed].atomnum[5]+1);
//			fprintf(input, "\n");
			}
		}
//		fprintf(input, "\n");
		// �e��ʊp�� V_2 �̋L�q
		for (nNumDihed=0; nNumDihed<prot.nNumDihedType; ++nNumDihed)
		{
			if (nNumDihed%5 == 0 && nNumDihed != 0)
			{
//				fprintf(input, "\n ");
			}
//			fprintf(input, "%e ", dataofdihed[nNumDihed].dihedangp.V_2);
		}
//		fprintf(input, "\n");
		// �e��ʊp�̈ʑ��̋L�q
		for (nNumDihed=0; nNumDihed<prot.nNumDihedType; ++nNumDihed)
		{
			if (nNumDihed%5 == 0 && nNumDihed != 0)
			{
//				fprintf(input, "\n ");
			}
//			fprintf(input, "%6d ", dataofdihed[nNumDihed].dihedangp.num);
		}
//		fprintf(input, "\n");
		// �e��ʊp�̑��d�x�̋L�q
		for (nNumDihed=0; nNumDihed<prot.nNumDihedType; ++nNumDihed)
		{
			if (nNumDihed%5 == 0 && nNumDihed != 0)
			{
//				fprintf(input, "\n ");
			}
//			fprintf(input, "%e ", dataofdihed[nNumDihed].dihedangp.theta);
		}
//		fprintf(input, "\n");
//		fprintf(input, "\n");
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
		// �񌋍������ݍ�p�̗͏���C���v�b�g�̋L�q
//		// �e���q�̓d�ׂ̋L�q
		for (nNumAtom=0; nNumAtom<prot.nNumAtom;++nNumAtom)
		{
			 fprintf(input, "%e ", prot.atm[nNumAtom].atomID.e_elect*ConstScaleOfCharge);

			if (nNumAtom%5 == 0 && nNumAtom != 0)
			{
				fprintf(input, "\n ");
			}
		}
		fprintf(input, "\n");

		fprintf(input, "%d \n", TotalAtomType);

		for (nNumAtom=0; nNumAtom<prot.nNumAtom;++nNumAtom)
		{
			fprintf(input, "%d ", prot.atm[nNumAtom].atomID.atomnum);

			if (nNumAtom%5 == 0 && nNumAtom != 0)
			{
				fprintf(input, "\n ");
			}
		}
		fprintf(input, "\n");

		for (nNumAtom=0; nNumAtom<5;++nNumAtom)
		{
			for (nNumAtom2=0; nNumAtom2<5;++nNumAtom2)
			{
				mat[nNumAtom][nNumAtom2]=0;
			}
		}

		i = 0;
		for (nNumAtom=0; nNumAtom<TotalAtomType;++nNumAtom)
		{
			for (nNumAtom2=0; nNumAtom2<=nNumAtom;++nNumAtom2)
			{
				++i;
				mat[nNumAtom][nNumAtom2]=i;
				if (nNumAtom!=nNumAtom2)
				{
					mat[nNumAtom2][nNumAtom]=mat[nNumAtom][nNumAtom2];
				}
			}
		}
		TypeTotal = i;

		for (nNumAtom=0; nNumAtom<TotalAtomType;++nNumAtom)
		{
			for (nNumAtom2=0; nNumAtom2<TotalAtomType;++nNumAtom2)
			{
				fprintf(input, "%3d ",mat[nNumAtom][nNumAtom2]);
//				printf("%d ",mat[nNumAtom][nNumAtom2]);
			}
			fprintf(input, "\n");
		}
		fprintf(input, "%d \n",TypeTotal);

		nNumAtom3=0;
		// �e���q�� A �̋L�q
		for (nNumAtom=0; nNumAtom<TotalAtomType;++nNumAtom)
		{
			for (nNumAtom2=0/*nNumAtom*/; nNumAtom2/*<*/<=nNumAtom/*TotalAtomType*/;++nNumAtom2)
			{
				length = prot.parmA[nNumAtom2]+prot.parmA[nNumAtom];
				depth = sqrt(prot.parmB[nNumAtom2]*prot.parmB[nNumAtom]);
				fprintf(input, "%e ", depth*pow(length,12));
				++nNumAtom3;
			if (nNumAtom3%5 == 0 && nNumAtom3 != 0)
			{
				fprintf(input, "\n ");
			}
			}
		}
		fprintf(input, "\n");

		nNumAtom3=0;
		// �e���q�� B �̋L�q
		for (nNumAtom=0; nNumAtom<TotalAtomType;++nNumAtom)
		{
			for (nNumAtom2=0/*nNumAtom*/; nNumAtom2/*<*/<=nNumAtom/*TotalAtomType*/;++nNumAtom2)
			{
				length = prot.parmA[nNumAtom2]+prot.parmA[nNumAtom];
				depth = sqrt(prot.parmB[nNumAtom2]*prot.parmB[nNumAtom]);
				fprintf(input, "%e ", 2.0*depth*pow(length,6));
			++nNumAtom3;
			if (nNumAtom3%5 == 0 && nNumAtom3 != 0)
			{
				fprintf(input, "\n ");
			}
			}
		}
		fprintf(input, "\n");
///////////////////////////////////////////////////////////////////////////////

	fprintf(input, "\n");

///////////////////////////////////////////////////////////////////////////////
	// �񌋍������ݍ�p���Ȃ��g�̗͏���C���v�b�g�̋L�q
	for (nNumAtom=0; nNumAtom<prot.nNumAtom; ++nNumAtom)
	{
//		ascatomnum(prot.atm[nNumAtom].not_interacting,prot.atm[nNumAtom].num_not_interacting);
		for (nNumAtom2=0;nNumAtom2<prot.atm[nNumAtom].num_not_interacting+1;++nNumAtom2)
		{
			if (prot.atm[nNumAtom].not_interacting[nNumAtom2] != 0)
			{
				fprintf(input, "%d ", prot.atm[nNumAtom].not_interacting[nNumAtom2]);
			}
		}
		fprintf(input, "0\n");
	}
///////////////////////////////////////////////////////////////////////////////

	fprintf(input, "\n");

///////////////////////////////////////////////////////////////////////////////
	// 1-4 �񌋍������ݍ�p�̗͏���C���v�b�g�̋L�q
	for (nNumAtom=0; nNumAtom<prot.nNumAtom; ++nNumAtom)
	{
//		ascatomnum(prot.atm[nNumAtom].o_f_interacting,prot.atm[nNumAtom].num_o_f_interacting);
		for (nNumAtom2=0;nNumAtom2<prot.atm[nNumAtom].num_o_f_interacting;++nNumAtom2)
		{
			if (prot.atm[nNumAtom].o_f_interacting[nNumAtom2] != 0)
			{
				fprintf(input, "%d ", prot.atm[nNumAtom].o_f_interacting[nNumAtom2]);
			}
		}
		fprintf(input, "0\n");
	}
///////////////////////////////////////////////////////////////////////////////

}

int CheckDihedPair(int atoms1[4], int atoms2[4])
{
	int i;

	for (i=0;atoms2[3-i]==atoms1[i];++i)
	{
		if(i==3)
		{
			break;
		}
	}

	if (i<3)
		return 0;
	else if (i==3 && atoms2[1] < atoms2[2])
		return 1;
	else
		return 0;





}

void ascatomnum(int atomnum[20],int num)
{
	int i,j,k;
	int atomnumdummy[20];
	int numdummy=1;

	for (i=0;i<20;++i)
	{
		atomnumdummy[i]=0;
	}

	for (i=0;i<num;++i)
	{
		for (j=0;j<num;++j)
		{
			if(atomnumdummy[j] == 0)
			{
				atomnumdummy[j] = atomnum[i];
				break;
			}
			else if (atomnum[i] < atomnumdummy[j])
			{
				for (k=num-j-1;k>=j;--k)
				{
					atomnumdummy[k+1]=atomnumdummy[k];
				}
				atomnumdummy[j]=atomnum[i];
				break;
			}
		}
	}

	for (i=0;i<num;++i)
	{
		atomnum[i]=atomnumdummy[i];
	}



}

void CreateIsouFile(FILE *input)
{
	int nNumAtom;
	int nNumAtom2;

	fprintf(input, "%d\n",prot.nNumAtom);

	///////////////////////////////////////////////////////////////////////////////
	// �񌋍������ݍ�p���Ȃ��g�̗͏���C���v�b�g�̋L�q
	for (nNumAtom=0; nNumAtom<prot.nNumAtom; ++nNumAtom)
	{
		for (nNumAtom2=0;;++nNumAtom2)
		{
			if (prot.atm[nNumAtom].bond_interacting[nNumAtom2] != 0)
			{
				fprintf(input, "%d %d\n",nNumAtom, prot.atm[nNumAtom].bond_interacting[nNumAtom2]-1);
			}
			else
			{
				break;
			}
		}
	}
///////////////////////////////////////////////////////////////////////////////
}

void CreatemolFile(FILE *input)
{
	int nNumAtom=0,nNumAtom2, nNumRes;

///,///////////////////////////////////////////////////////////////////////////////
	for (nNumRes=0; nNumRes<prot.nNumResidue; ++nNumRes)
	{
		for (nNumAtom2=0;nNumAtom2<prot.numsequence[nNumRes]-1;++nNumAtom2)
		{
			fprintf(input, "%d\n",prot.atm[nNumAtom].atomID.atomnum);
		}
	}
///,///////////////////////////////////////////////////////////////////////////////
}


