#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "Vis_MD.h"

// �֐�
void atomscanf( FILE *dataf, /* residuedata dataofthisresidue */int NumRes);
void clustscanf(FILE *dataf, /* residuedata dataofthisresidue */int NumRes);
void dihedscanf(FILE *dataf, elementoftable *data, residuedata dataofthisresidue);

// �f�[�^�̎擾���s���֐�
int PickData(void) {
  int nNumRes;
  char	nameofthisresidue[MAXNRESIDUE][4];
  int	flag;
  FILE *datafile;
  
  // �W���c��f�[�^�x�[�X���J��
  datafile = efopen("stdres.dat","r");
  
  // �W���c��f�[�^�x�[�X�̓ǂݍ���
  for (nNumRes=0 ; nNumRes < NRESIDUE ; ++nNumRes){
    // �c���ǂ�
    fscanf(datafile, "%s", &(nameofthisresidue[nNumRes]));
    nameofthisresidue[nNumRes][3] = '\0';
    // �c��̌��q����ǂ�
    fscanf(datafile, "%d", &(dataofresidueontable[nNumRes].nNumAtom));
    // �c��̌��q�̃f�[�^��ǂ�
    atomscanf(datafile,/* dataofresidueontable[nNumRes] */nNumRes);
    // �c��̍��̐���ǂ�
    fscanf(datafile, "%d", &(dataofresidueontable[nNumRes].nNumClut));
    // �c��̍��̂̃f�[�^��ǂ�
    clustscanf(datafile,/* dataofresidueontable[nNumRes]*/nNumRes);
    // �e�[�u�����Ƀf�[�^��o�^
    LURData(nameofthisresidue[nNumRes],0,dataofresidueontable[nNumRes]);
    fscanf(datafile, "%d", &flag);

    if (++NumResInStdResData > NRESIDUE) {
      printf("too many residues in stdres.dat !!!\n");
      break;
    }
    
    if (flag == 0)
      break;
  }

  // �W���c��f�[�^�x�[�X�����
  fclose(datafile);
}

// �W���c��̌��q�f�[�^�̓ǂݍ��݂��s���֐�
void atomscanf(FILE *dataf, /* residuedata dataofthisresidue */int nNumRes ) {
  int alpha;
  int nNumAtom;
  int nNumAtom2;

  int nNum_not_interacting;
  int nNum_not_interacting_kminusone;
  int nNum_not_interacting_kplusone;
  
  int ThisnNumAtom;
  
  int nNumAtomOfResidue;
  
  double dData;
  double e_elect;
  double eata;
  double gamma;
  int nData;
  char sData;

  nNumAtomOfResidue = dataofresidueontable[nNumRes].nNumAtom;

  // �ʏ팴�q�̃f�[�^�̓ǂݍ���
  for (nNumAtom = 0; nNumAtom < nNumAtomOfResidue; ++nNumAtom) {
    // ���q�ԍ��̓ǂݍ���
    fscanf(dataf, "%d", &ThisnNumAtom);
    // �_�~�[���q�̏ꍇ�A���q���̓i�V
    if (ThisnNumAtom != nNumAtomOfResidue) 
      // ���q���̓ǂݍ���
      fscanf(dataf, "%d", &(dataofresidueontable[nNumRes].atm[nNumAtom].atomID.atomnum));
    
    // head or tail ���q�̓ǂݍ���
    fscanf(dataf, "%d", &nData);
    if (nData == 1)
      dataofresidueontable[nNumRes].HeadAtom = nNumAtom+1;
    else if (nData == -1)
	dataofresidueontable[nNumRes].TailAtom = nNumAtom+1;

    // ���W�f�[�^�̓ǂݍ���
    for (alpha = 0; alpha < 3; ++alpha)
      fscanf(dataf, "%lf", &(dataofresidueontable[nNumRes].atm[nNumAtom].coord[alpha]));

    // �_�~�[���q�̏ꍇ�A���q���̓i�V
    if (ThisnNumAtom != nNumAtomOfResidue) {
      // �d�ׂ̓ǂݍ���
      fscanf(dataf, "%lf", &(dataofresidueontable[nNumRes].atm[nNumAtom].atomID.e_elect));
      // �K���}�̓ǂݍ���
      fscanf(dataf, "%lf", &(dataofresidueontable[nNumRes].atm[nNumAtom].atomID.gamma));
      // �C�[�^�̓ǂݍ���
      fscanf(dataf, "%lf", &(dataofresidueontable[nNumRes].atm[nNumAtom].atomID.eata));      
      fscanf(dataf, "%d", &(dataofresidueontable[nNumRes].atm[nNumAtom].num_bond_interacting));
      // �����̌��q�̓ǂݍ���
      for (nNumAtom2=0;nNumAtom2<dataofresidueontable[nNumRes].atm[nNumAtom].num_bond_interacting;++nNumAtom2) {
	fscanf(dataf, "%d", &(dataofresidueontable[nNumRes].atm[nNumAtom].bond_interacting[nNumAtom2]));
      }
    }
  }
}


// �W���c��̍��̃f�[�^�̓ǂݍ��݂��s���֐�
void clustscanf(FILE *dataf, /* residuedata dataofthisresidue */int nNumRes ) {
  int nNumClt;
  int nNumBh;
  int nNumCltOfResidue;
  
  nNumCltOfResidue = dataofresidueontable[nNumRes].nNumClut;
  // ���̃f�[�^�̓ǂݍ���
  for (nNumClt = 0; nNumClt < nNumCltOfResidue; ++nNumClt) {
    // ���̒��̌��q���̓ǂݍ���
    fscanf(dataf, "%d", &(dataofresidueontable[nNumRes].clt[nNumClt].nNumAtom));    
    // ���̒��̎}�̐��̓ǂݍ���
    fscanf(dataf, "%d", &(dataofresidueontable[nNumRes].clt[nNumClt].nNumBranch));        
    // ���̖̂��[���̓ǂݍ���
    fscanf(dataf, "%d", &(dataofresidueontable[nNumRes].clt[nNumClt].termflag));    
    // ���̒��̌��_���̓ǂݍ���
    fscanf(dataf, "%d", &(dataofresidueontable[nNumRes].clt[nNumClt].origin));
    // ���̒��̏I�_���̓ǂݍ���
    for (nNumBh = 0; nNumBh < dataofresidueontable[nNumRes].clt[nNumClt].nNumBranch; ++nNumBh)
      fscanf(dataf, "%d", &(dataofresidueontable[nNumRes].clt[nNumClt].term[nNumBh]));	
    // "�e"�̍��̂̔ԍ��̓ǂݍ���
    fscanf(dataf, "%d", &(dataofresidueontable[nNumRes].clt[nNumClt].nNumParentClust));    
    // "�q"�̍��̂̔ԍ��̓ǂݍ���
    for (nNumBh = 0; nNumBh < dataofresidueontable[nNumRes].clt[nNumClt].nNumBranch; ++nNumBh)
      fscanf(dataf, "%d", &(dataofresidueontable[nNumRes].clt[nNumClt].nNumChildClust[nNumBh]));
  }
}

// �W���c��̓�ʊp�f�[�^�̓ǂݍ��݂��s���֐�
void dihedscanf(FILE *dataf,elementoftable *data, residuedata dataofthisresidue) {
  int nNumDihed;
  int nNumAtm;
  double dData;
  int nData;
  int nNumDihedOfResidue;

  fscanf(dataf, "%d", &(prot.nNumDihedType));

  // ��ʊp�f�[�^�̓ǂݍ���
  for (nNumDihed = 0; nNumDihed < nNumDihedOfResidue; ++nNumDihed) {
    fscanf(dataf, "%d", &nData);    
    // ��ʊp���ݍ�p��V_2�̓ǂݍ���
    fscanf(dataf, "%lf", &(dataofdihed[nNumDihed].dihedangp.V_2));
    // ��ʊp���ݍ�p�̑��d�x�̓ǂݍ���
    fscanf(dataf, "%d", &(dataofdihed[nNumDihed].dihedangp.num));
    // ��ʊp���ݍ�p�̈ʑ��̓ǂݍ���
    fscanf(dataf, "%lf", &(dataofdihed[nNumDihed].dihedangp.theta));
    for (nNumAtm=0;nNumAtm<4;++nNumAtm)	{
      fscanf(dataf, "%d", &(dataofdihed[nNumDihed].atomnum[nNumAtm]));
    }
  }
}

