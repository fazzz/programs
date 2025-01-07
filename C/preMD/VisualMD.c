#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "Vis_MD.h"

int TotalAtomType;

// �f�[�^�̎擾�A�C���v�b�g���̍쐬�A�C���v�b�g�t�@�C���̍쐬���s�����C���֐�
int main(int argc, char *argv[]) {
  int i;
  int j;
  int c;
  int num;
  int opts=3;
  int nNumTotalResdiue=0;
  char *nameofthisresidue;
  char nameofthisresidue2[MAXNRESIDUE][4];
  int opt;
  char *option;
  FILE *inputfile;
  char *progname;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];

  InpfilCOORD = "crd.in";
  InpfilCLUST = "clust.in";
  InpfilSEQ = "seq.in";
  InpfilTOP = "top.in";

  // �f�[�^�̎擾
  PickData();

  // �c��z��̏����擾
  if (argc > 1) {
    if (strcmp(argv[1],"-h")==0)
      ADDFLAG = 1;
    else if (strcmp(argv[1],"-l")==0)
      ADDFLAG = 0;
    else {
      usage(progname);
      exit(1);
    }
    if (strcmp(argv[2],"-i")==0) {
      InpfilCLUST=argv[3];
      opts=4;
    }
    else {
      opts=2;
    }
    if (strcmp(argv[opts],"-s")==0) {
      for (i=opts+1; i<argc; ++i) {
	prot.Sequence[i-opts-1] = argv[i];
	++nNumTotalResdiue;
      }
      prot.Sequence[i] = '\0';
    }
    else if (strcmp(argv[opts],"-f")==0) {
      inputfile=efopen(argv[opts+1],"r");
      for(i=0;fscanf(inputfile,"%3s",&nameofthisresidue2[i]) != EOF;++i) {
	nameofthisresidue2[i][3] = '\0';
	prot.Sequence[i] = nameofthisresidue2[i];
      }
      nNumTotalResdiue = i;
      fclose(inputfile);
    }
    else {
      usage(progname);
      exit(1);
    }
  }
  // �c��z��̏�񂪂Ȃ��ꍇ
  else {
    usage(progname);
    exit(1);
  }
	  
  // �^���p�N�����̎c��̎擾
  prot.nNumResidue = nNumTotalResdiue;
  // �C���v�b�g�t�@�C���̃f�[�^�̍쐬
  TotalAtomType=CreateData();
  // �C���v�b�g�t�@�C���̍쐬
  CreateInPut(TotalAtomType);
  
  for (num=0;num<prot.nNumResidue;++num) {
    printf("****");
  }
  printf("***************************");
  printf("\n");
  for (num=0;num<prot.nNumResidue;++num) {
    printf("**");
  }
  printf(" the input files are made: ");
  for (num=0;num<prot.nNumResidue;++num) {
    printf("**");
  }
  printf("\n");
  printf("***** the sequence is ");
  for (num=0;num<prot.nNumResidue;++num) {
    printf("%s ", prot.Sequence[num]);
  }
  printf("*****\n");
  for (num=0;num<prot.nNumResidue;++num) {
    printf("****");
  }
  printf("***************************");
  printf("\n");

  return 1;
}
