#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXNAC /*20*//*300*/600 /* ¹äÂÎ¿ô¤Î¾å¸Â*/
#define MAXA /*1000*/2000 /* ¸¶»Ò¿ô¤Î¾å¸Â*/
#define MAXDOF /*300*/600 /* ¼«Í³ÅÙ¤Î¾å¸Â*/
#define MAXNBH 4 /* Ê¬´ô¿ô¤Î¾å¸Â*/
#define MAXNDAC /*10*/20 /* ¹äÂÎÃæ¤Î¸¶»Ò¿ô¤Î¾å¸Â*/
#define MAXPDB 10000 /* ½ÐÎÏ PDB ¥Õ¥¡¥¤¥ë¿ô¤Î¾å¸Â*/
#define MAXPEPTIDE 10
#define ON 1
#define OFF 0

struct{
	int origin_atom_a;
	int terminal_atom_a[MAXNBH];
	int terminal;
	int num_clust;
	int num_atom_clust;
	int num_branch;
	int join;/*0410*/

	int nNumClutOfParent;
	int nNumClutOfChild[MAXNBH];

	int origin_xoord_a;
	int num_xoord_a;

	double xoord_clust[MAXNBH][MAXA][3];
	double mass_clust[MAXNAC];
	double sum_mass;
 	double Inertia_clust[3][3];
	double InertiaMatrix[6][6];
	double Coriolis_acc_Mat[3][3];
	double momentum_clust;
	double Inertia_clust_total;
	double PsedoInertia[4][4];
	double PsedoTransMatrix[4][4];
	double qCOM[3];
	double dihedang[MAXNBH];
	double ddihedang[MAXNBH];
	double dddihedang[MAXNBH];
	double correct_dihedang[6];
	double predict_dihedang[6];
	double now_deltadihedang[MAXNBH];
	double Hing[MAXNBH][6];
	double TransMatrix[MAXNBH][6][6];
	double trans_A_to_CN[MAXNBH][3][3];

	double sp_velo[6];
	double Coriolis_acc[6];
	double Coriolis_b[6];
	double sp_acc[6];
	double predict_alpha[6];
	double predict_velo[6];
} clust[MAXDOF];

struct protein{
	char *name_prot;
	int num_atom;
	int nNumPeptide;
	int nNumClutPeptide[MAXPEPTIDE];
	int inumrs;
	int num_clust;
	int DOF;
	int nNumDihedALL;
	int nNumDihedType;

        int nNumDihed_rest;

	double coord[MAXA][3];
	double velo[MAXA][3];
	double qCOM[3];
	double qCOM_Old[3];
	double veloCOM[3];
	double sumMass;
	int name_atom[MAXA];
}prot;

double delta_matrix[4][4];


/**********************/
/* char* InpfilCOORD; */
/* char* InpfilCLUST; */
/* char* InpfilDVELO; */
/**********************/

void calc_velo2(double *dotTransQ, double *velo);
void calc_dot_Pseduo_TransMatrix(int nNumClutI,int nNumClutK,double mat[4][4]);
void trans_A_to_CN(int nNumClut);
void sub_trans_A_to_CN(int nNumCltTar, int nNumCltCoo,int nNumBod, int nNumAtom);
void sub_trans_A_to_CN_Initial(void);
int coordscan(FILE *input);
void clustscan(FILE *input);
void pick_initial_velo(char *InpfilDVELO);
int setJoin(int nNumClut, int joinflag);
void set_trans_Matrix(int nNumClt,int nNumClutOrigBranch);
void sub_set_trans_Matrix(int nNumClt,int nNumCltminousone);
void set_pseduo_trans_matrix(int nNumClut);
void set_delts_matrix(void);

double ident[6][6];

int main(int argc, char *argv[]) {
  int i,j;
  FILE* input;
  char *line;
  size_t len=0;
  int joinflag;
  int nNumClut,nNumAtomClut,nNumAtom,alpha;
  double Energy_kinetic6;

  char* InpfilCOORD;
  char* InpfilCLUST;
  char* InpfilDVELO;

  double *dotTransQ;
  double *velo;

  if (argc < 3){
    printf("error : error inout option\n");
  }
  InpfilCOORD=*++argv;
  InpfilCLUST=*++argv;
  InpfilDVELO=*++argv;
  
  for (i=0;i<6;++i) {
    for (j=0;j<6;++j) {
    ident[i][j]=0.0;
    }
  }
  for (i=0;i<6;++i) {
    ident[i][i]=1.0;
  }
  
  if ((input=fopen(InpfilCOORD,"r")) == NULL) {
      printf("error %s cannot open\n", InpfilCOORD);
      exit(1);
  }
  getline(&line,&len,input);
  prot.num_atom = coordscan(input);
  fclose(input);
  
  if ((input=fopen(InpfilCLUST,"r")) == NULL) {
    printf("error %s cannot open\n", InpfilCLUST);
    exit(1);
  }
  clustscan(input);
  fclose(input);
  
  pick_initial_velo(InpfilDVELO);
  
  joinflag=0;
  for(nNumClut=0; nNumClut<prot.DOF; ++nNumClut) {
    joinflag=setJoin(nNumClut,joinflag);
  }
  
  for(nNumClut=0; nNumClut<prot.DOF; ++nNumClut){
      trans_A_to_CN(nNumClut);
  }
  double matrix[1000][20000][4];

  printf("m:162 dotTransQ %p \n",dotTransQ);
  printf("m:162 velo %p \n",velo);
  if((dotTransQ=(double *)malloc(sizeof(double)*prot.DOF*prot.num_atom*4)
// **10*20OK*//**100*100ERROR*//**50*100ERROR*//**50*100*//**203*520*/*4)
     )==NULL) {
    printf("error: malloc\n");
    exit(1);
  }
  printf("m:167 dotTransQ %p \n",dotTransQ);
  if((velo=(double *)malloc(sizeof(double)/**prot.num_atom*/*100*4))==NULL) {
    printf("error: malloc\n");
    exit(1);
  }
  printf("m:170 velo %p \n",velo);

  /*******************************/
  /* free(velo);		 */
  /* printf("velo is OK !!!\n"); */
  /* free(dotTransQ);		 */
  /*******************************/
  /************/
  /* exit(1); */
  /************/

  calc_velo2(dotTransQ,velo);
  Energy_kinetic6 = 0.0;
  for (nNumClut=0;nNumClut<prot.DOF;++nNumClut) {
    for (nNumAtomClut=0;nNumAtomClut<clust[nNumClut].num_atom_clust;++nNumAtomClut) {
      for (alpha=0;alpha<3;++alpha) {
  	Energy_kinetic6 += 0.5*clust[nNumClut].mass_clust[nNumAtomClut]*prot.velo[nNumAtom][alpha]*prot.velo[nNumAtom][alpha];
      }
      ++nNumAtom;
    }
  }
  Energy_kinetic6 = Energy_kinetic6/(4.18407*100.0);
}

void calc_velo2(double *dotTransQ, double *velo) {
  int i,ii,i_c,alpha,alpha2,alpha3,j,k,n;

  int nNumClut,nNumClutI,nNumClutJ,nNumCltCoo,nNumClutdummy;
  int nNumAtom,nNumAtomtotal,nNumAtomClut;
  int nNumOrigc, nNumAbsoc;

  double mat[4][4];
  double dummy[4][4];
  double dummy_xoord[4];
  double dummy_xoord2[4];
  //double *dotTransQ;
  //double *velo;

  printf("%p \n",dotTransQ);

  /***********************************************************************************************************************************************/
  /* if((dotTransQ=(double *)malloc(sizeof(double)*prot.DOF*prot.num_atom*4))==NULL) {                                printf("error: malloc\n"); */
  /*   exit(1);																	 */
  /*   }																	 */
  /* printf("%p \n",dotTransQ);															 */
  /* if((velo=(double *)malloc(sizeof(double)*prot.num_atom*4))==NULL) {									 */
  /*     printf("error: malloc\n");														 */
  /*     exit(1);																 */
  /*   }																	 */
  /***********************************************************************************************************************************************/

  for (nNumClutI=/*0*/1;nNumClutI<prot.DOF;++nNumClutI) {
    set_pseduo_trans_matrix(nNumClutI);
  }
  printf("a:1114\n");

  for (nNumClut=0;nNumClut<prot.DOF;++nNumClut) {
    for (nNumAtom=0;nNumAtom<prot.num_atom;++nNumAtom)	{
      for (i=0;i<4;++i) {
	printf("%d %d %d %d\n",nNumClut,nNumAtom,i,nNumClut*prot.num_atom*4+nNumAtom*4+i);
	printf("%p \n",dotTransQ);
	dotTransQ[nNumClut*prot.num_atom*4+nNumAtom*4+i] = 0.0;
      }
    }
  }
  printf("a:211 \n");

  for (nNumClutI=0;nNumClutI<prot.DOF;++nNumClutI) {
    nNumAtom=0;
    for (nNumClutJ=0;nNumClutJ<prot.DOF;++nNumClutJ) {
      if (/*nNumClutI != 0 && */nNumClutI<=nNumClutJ && clust[nNumClutI].join <= clust[nNumClutJ].join) {
	for (i=0;i<3;++i) {
	  for (j=0;j<3;++j) {
	    mat[i][j] = 0.0;
	  }
	}
	printf("a:222 \n");
	calc_dot_Pseduo_TransMatrix(nNumClutI,nNumClutJ,mat);
	printf("a:1135\n");
      }
      else {
	for (i=0;i<3;++i) {
	  for (j=0;j<3;++j) {
	    mat[i][j] = 0.0;
	  }
	}
      }
      printf("a:1143\n");
      for (nNumAtomClut=0;nNumAtomClut<clust[nNumClutJ].num_atom_clust;++nNumAtomClut) {
	if (/*nNumClutI != 0 &&*/ nNumClutI<=nNumClutJ && clust[nNumClutI].join <= clust[nNumClutJ].join) {

	  for (i=0;i<3;++i)	{
	    dummy_xoord[i] = clust[nNumClutJ].xoord_clust/*[0]*/[nNumAtomClut][i];
	  }
	  dummy_xoord[3] = 1.0;

	  for (i=0;i<4;++i) {
	    for (j=0;j<4;++j) {
	      dotTransQ[nNumClutI*prot.num_atom*4+nNumAtom*4+i]+= mat[i][j]*dummy_xoord[j];
	    }
	  }
	}
	else {
	  for (i=0;i<4;++i) {
	    for (j=0;j<4;++j) {
	      dotTransQ[nNumClutI*prot.num_atom*4+nNumAtom*4+i]=0.0;
	    }
	  }
	}
	++nNumAtom;
      }
    }
  }
  printf("a:1174\n");

  for (nNumAtom=0;nNumAtom<prot.num_atom;++nNumAtom){
    for (alpha=0;alpha<4;++alpha){
      velo[nNumAtom*4+alpha] = 0.0;
    }
  }

  for(nNumClutJ=0;nNumClutJ<prot.DOF;++nNumClutJ) {
    for (nNumAtom=0;nNumAtom<prot.num_atom;++nNumAtom) {
      for (alpha=0;alpha<3;++alpha) {
	velo[nNumAtom*4+alpha] += dotTransQ[nNumClutJ*prot.num_atom*4+nNumAtom*4+alpha]*clust[nNumClutJ].ddihedang[0];
      }
    }
  }
  printf("a:1189\n");

  printf("a:1196\n");

  free(dotTransQ);
  free(velo);
}

void calc_dot_Pseduo_TransMatrix(int nNumClutI,int nNumClutK,double mat[4][4]) {
  int i,j,k;
  int flag;
  int nNumClut,nNumClutdummy,nNumClutLast;
  double dummy_matrix[4][4], dummy_matrix2[4][4];

  flag=ON;
  if (clust[nNumClutI].join>0) {
    for (nNumClutdummy=nNumClutI;clust[nNumClutdummy].join!=clust[nNumClutI].join-1;++nNumClutdummy) {
      nNumClutLast = nNumClutdummy;
    }
    if (nNumClutLast < nNumClutK) {
      flag=OFF;
    }
  }
  for (i=0;i<4;++i) {
    for (j=0;j<4;++j) {
      mat[i][j] = 0.0;
    }
  }
  printf("Mat:nNumClutI=%d,nNumClutK=%d 138\n",nNumClutI,nNumClutK);

  if (flag==ON) {
    for (i=0;i<4;++i){
      for (j=0;j<4;++j){

dummy_matrix[i][j] = ident[i][j];
      }
    }

    for (nNumClut=1/*0*/;nNumClut<=nNumClutI;/* ++nNumClut */){
      if (clust[nNumClut].join <= clust[nNumClutI].join) {
	for (i=0;i<4;++i){
	  for (j=0;j<4;++j)	{
	    dummy_matrix2[i][j] = dummy_matrix[i][j];
	  }
	}
	printf("Mat:nNumClut=%d 155\n",nNumClut);

	for (i=0;i<4;++i){
	  for (j=0;j<4;++j){
	    dummy_matrix[i][j] = 0.0;
	  }
	}
	printf("Mat:161\n");

        /*************************************************************/
        /* printf("nNumClut=%d, nNumClutI=%d\n",nNumClut,nNumClutI); */
        /*************************************************************/

	for (i=0;i<4;++i){
	  for (j=0;j<4;++j){
	    for (k=0;k<4;++k){
	      dummy_matrix[i][j] += dummy_matrix2[i][k]*clust[nNumClut].PsedoTransMatrix[k][j];
	    }
	  }
	}
	printf("Mat:174\n");
      }
      printf("Mat:177\n");

      if (clust[nNumClut].num_branch >1 ) {
	if (clust[nNumClut].nNumClutOfChild[1]-1 <= nNumClutI) {
	  nNumClut = clust[nNumClut].nNumClutOfChild[1]-1-1;
	}
      }
      ++nNumClut;
    }
    printf("Mat:185\n");

    for (i=0;i<4;++i) {
      for (j=0;j<4;++j) {
	dummy_matrix2[i][j] = dummy_matrix[i][j];
      }
    }

    for (i=0;i<4;++i) {
      for (j=0;j<4;++j) {
	dummy_matrix[i][j] = 0.0;
      }
    }

    for (i=0;i<4;++i) {
      for (j=0;j<4;++j) {
	for (k=0;k<4;++k) {
	  dummy_matrix[i][j] += dummy_matrix2[i][k]*delta_matrix[k][j];
	}
      }
    }
    printf("Mat:204\n");

    nNumClut=nNumClutI;
    if (clust[nNumClut].num_branch >1 ) {
      if (clust[nNumClut].nNumClutOfChild[1]-1 <= nNumClutK) {
	nNumClut = clust[nNumClut].nNumClutOfChild[1]-1-1;
      }
    }
    ++nNumClut;
    for (/*nNumClut=nNumClutI+1*/;nNumClut<=nNumClutK;/* ++nNumClut */) {
      if (clust[nNumClut].join <= clust[nNumClutK].join) {
	for (i=0;i<4;++i) {
	  for (j=0;j<4;++j)	{
	    dummy_matrix2[i][j] = dummy_matrix[i][j];
	  }
	}

	for (i=0;i<4;++i){
	  for (j=0;j<4;++j){
	    dummy_matrix[i][j] = 0.0;
	  }
	}
	printf("Mat:226\n");

        /*************************************************************/
        /* printf("nNumClut=%d, nNumClutK=%d\n",nNumClut,nNumClutK); */
        /*************************************************************/

	for (i=0;i<4;++i){
	  for (j=0;j<4;++j)	{
	    for (k=0;k<4;++k) {
	      dummy_matrix[i][j] += dummy_matrix2[i][k]*clust[nNumClut].PsedoTransMatrix[k][j];
	    }
	  }
	}
	printf("Mat:239\n");
      }
      if (clust[nNumClut].num_branch >1 ) {
	if (clust[nNumClut].nNumClutOfChild[1]-1 <= nNumClutK) {
	  nNumClut = clust[nNumClut].nNumClutOfChild[1]-1-1;
	}
      }
      printf("Mat:246\n");
      ++nNumClut;
    }

    for (i=0;i<4;++i) {
      for (j=0;j<4;++j){
	mat[i][j] = dummy_matrix[i][j];
      }
    }
    printf("Mat:256\n");
  }
}

// 実験室系→局所座標系の変換を行う関数
void trans_A_to_CN(int nNumClut)
{
	int nNumClutOrigBranch;
	int nNumClutOfParent;
	int nNumClut2;
	int num_atom;
	int num,nNumClutdummy;
	int flag;

	nNumClutOfParent = clust[nNumClut].nNumClutOfParent-1;

	if (clust[nNumClut].num_branch > 1)
	{
		nNumClutOrigBranch = nNumClut;
	}

	// 0番目のクラスタでは、原点の移動を行う
	if (nNumClut == 0)
	{
		sub_trans_A_to_CN(0, -1,0,
					       prot.num_atom);
		sub_trans_A_to_CN_Initial();
	}
	// In Side Chain
	else if(clust[nNumClut].join > 0)
	{
	  num=0;
	  for (nNumClutdummy=nNumClut;clust[nNumClutdummy].join!=clust[nNumClut].join-1;++nNumClutdummy) {
	    num+=clust[nNumClutdummy].num_atom_clust;
	  }
	  sub_trans_A_to_CN(nNumClut, nNumClutOfParent,0,num);
	  }
	// In Main Chain
	else {
	  sub_trans_A_to_CN(nNumClut, nNumClutOfParent, 0,
			    prot.num_atom-clust[nNumClut].origin_atom_a+1);
	}
}

// 実験室系→局所座標系の変換の補助を行う関数
void sub_trans_A_to_CN(int nNumCltTar, int nNumCltCoo,
                       int nNumBod, int nNumAtom)
{
	int i;

	int alpha;

	int nNmAtomOfCN_A;
	int nNmAtomOfCN_1_A;

	double CN_A[3];
	double HN_A[3];
	double CN_1_A[3];
	double mat[MAXA][3];

	double ii[3];
	double jj[3];
	double kk[3];

	double jj_x;
	double jj_y;
	double jj_z;

	double DisOfCN_1_CN=0.0;
	double DisOfHN_CN=0.0;
	double SnCN_1_CN_HN=0.0;
	double CsCN_1_CN_HN=0.0;

	clust[nNumCltTar].num_xoord_a = nNumAtom;

	// 原子番号の取得
	nNmAtomOfCN_A = clust[nNumCltTar].origin_atom_a-1;
	if (nNumCltCoo != -1)
	{
		nNmAtomOfCN_1_A = clust[nNumCltCoo].terminal_atom_a[0]-1;
	}
	else
	{
		nNmAtomOfCN_1_A = 0;
	}

	clust[nNumCltTar].origin_xoord_a = 0/*nNmAtomOfCN_A-nNmAtomOfCN_1_A*/;

	// 実験室系での座標の取得
	for(alpha=0;alpha<3;++alpha)
	{
		CN_A[alpha]/*A*/=prot.coord[nNmAtomOfCN_A][alpha]/*A*/;
		HN_A[alpha]/*A*/=prot.coord[nNmAtomOfCN_A+1][alpha]/*A*/;
		if (nNumCltCoo != -1)
		{
			CN_1_A[alpha]/*A*/=prot.coord[nNmAtomOfCN_1_A][alpha]/*A*/;
		}
		else
		{
			CN_1_A[alpha]/*A*/=0.0/*A*/;
		}
	}


	for(alpha=0;alpha<3;++alpha)
	{
		// ^0k1の計算_1
		kk[alpha]/*A*/ = CN_1_A[alpha]-CN_A[alpha]/*A*/;
		// |CN_1-CN|の計算_1
		DisOfCN_1_CN
		 += (CN_1_A[alpha]-CN_A[alpha])*(CN_1_A[alpha]-CN_A[alpha]);
		// |CN_1-CN|の計算_1
		DisOfHN_CN
		 += (HN_A[alpha]-CN_A[alpha])*(HN_A[alpha]-CN_A[alpha]);
		// 角(CN_1-CN-HN)の計算_1
		CsCN_1_CN_HN += (CN_1_A[alpha]-CN_A[alpha])*(HN_A[alpha]-CN_A[alpha]);
	}

	// |CN_1-CN|の計算_2
	DisOfCN_1_CN/*A*/ = sqrt(DisOfCN_1_CN)/*A*/;
	// |CN_1-CN|の計算_2
	DisOfHN_CN/*A*/ = sqrt(DisOfHN_CN)/*A*/;
	// 角(CN_1-CN-HN)の計算_2
	CsCN_1_CN_HN = /*-*/CsCN_1_CN_HN/(DisOfCN_1_CN*DisOfHN_CN);
	SnCN_1_CN_HN = 1.0-CsCN_1_CN_HN*CsCN_1_CN_HN;
	SnCN_1_CN_HN = sqrt(SnCN_1_CN_HN);

	// ^0k1の計算_2
	for(alpha=0;alpha<3;++alpha)
	{
		kk[alpha]=kk[alpha]/DisOfCN_1_CN;
	}

	// ^0j1の計算
	jj[0]=((CN_1_A[1]-CN_A[1])*(HN_A[2]-CN_A[2])
		  -(CN_1_A[2]-CN_A[2])*(HN_A[1]-CN_A[1]))
		  /(DisOfCN_1_CN*DisOfHN_CN*SnCN_1_CN_HN);
	jj[1]=((CN_1_A[2]-CN_A[2])*(HN_A[0]-CN_A[0])
	      -(CN_1_A[0]-CN_A[0])*(HN_A[2]-CN_A[2]))
	      /(DisOfCN_1_CN*DisOfHN_CN*SnCN_1_CN_HN);
	jj[2]=((CN_1_A[0]-CN_A[0])*(HN_A[1]-CN_A[1])
	      -(CN_1_A[1]-CN_A[1])*(HN_A[0]-CN_A[0]))
	      /(DisOfCN_1_CN*DisOfHN_CN*SnCN_1_CN_HN);

	ii[0]=jj[1]*kk[2]-jj[2]*kk[1];
	ii[1]=jj[2]*kk[0]-jj[0]*kk[2];
	ii[2]=jj[0]*kk[1]-jj[1]*kk[0];

	// 座標系の変換_原点の重ね合わせ
	for(i=0; i<nNumAtom; ++i)
	{
		for(alpha=0;alpha<3;++alpha)
		{
		  mat[i][alpha] = prot.coord[/*nNmAtomOfCN_1_A*/nNmAtomOfCN_A+i][alpha]/*A*/
			              -CN_A[alpha]/*A*/;
		}
	}

	// 座標系の変換行列の作成
	for(alpha=0;alpha<3;++alpha)
	{
		clust[nNumCltTar].trans_A_to_CN[nNumBod][0][alpha]=ii[alpha];
		clust[nNumCltTar].trans_A_to_CN[nNumBod][1][alpha]=jj[alpha];
		clust[nNumCltTar].trans_A_to_CN[nNumBod][2][alpha]=kk[alpha];
	}

	// 座標系の変換_変換行列の乗算
	for(i=0; i < nNumAtom; ++i)
	{
		clust[nNumCltTar].xoord_clust[nNumBod][i][0]
		 =   clust[nNumCltTar].trans_A_to_CN[nNumBod][0][0]*mat[i][0]
		   + clust[nNumCltTar].trans_A_to_CN[nNumBod][0][1]*mat[i][1]
		   + clust[nNumCltTar].trans_A_to_CN[nNumBod][0][2]*mat[i][2];

		clust[nNumCltTar].xoord_clust[nNumBod][i][1]
		 =   clust[nNumCltTar].trans_A_to_CN[nNumBod][1][0]*mat[i][0]
		   + clust[nNumCltTar].trans_A_to_CN[nNumBod][1][1]*mat[i][1]
		   + clust[nNumCltTar].trans_A_to_CN[nNumBod][1][2]*mat[i][2];

		clust[nNumCltTar].xoord_clust[nNumBod][i][2]
		 =   clust[nNumCltTar].trans_A_to_CN[nNumBod][2][0]*mat[i][0]
		   + clust[nNumCltTar].trans_A_to_CN[nNumBod][2][1]*mat[i][1]
		   + clust[nNumCltTar].trans_A_to_CN[nNumBod][2][2]*mat[i][2];
	}
}

// 0番目の剛体の実験室系→局所座標系の変換の補助を行う関数
void sub_trans_A_to_CN_Initial(void)
{
	int i;
	int alpha, alpha2;

	// 単位行列の設定を行う
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

	// 座標系の変換
	for (i=0;i<prot.num_atom;++i)
	{
		for (alpha=0;alpha<3;++alpha)
		{
//			prot.coord[i][alpha] = clust[0].xoord_clust/*[0]*/[i][alpha];
		}
	}

	for(alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
			clust[0].trans_A_to_CN[0][alpha][alpha2]=0.0;
		}
	}

	// 座標系の変換行列の作成
	for(alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
			clust[0].trans_A_to_CN[0][alpha][alpha2]=ident[alpha][alpha2];
		}
	}

	clust[0].trans_A_to_CN[0][0][0]=1.0;
	clust[0].trans_A_to_CN[0][0][1]=0.0;
	clust[0].trans_A_to_CN[0][0][2]=0.0;
	clust[0].trans_A_to_CN[0][1][0]=0.0;
	clust[0].trans_A_to_CN[0][1][1]=1.0;
	clust[0].trans_A_to_CN[0][1][2]=0.0;
	clust[0].trans_A_to_CN[0][2][0]=0.0;
	clust[0].trans_A_to_CN[0][2][1]=0.0;
	clust[0].trans_A_to_CN[0][2][2]=1.0;
}

// 座標データの取得を行う関数
int coordscan(FILE *input)
{
	double x;
	int i,j;
	int flag = 0;
	int num = 0;

	fscanf(input,"%d",&i);

	for(i=0;i<MAXA;++i)
	{
		for(j=0;j<3;++j)
		{
			if (fscanf(input,"%lf",&x) == EOF)
			{
				flag = 1;
				break;
			}
			prot.coord[i][j]=x;
		}
		if (flag == 1)
			break;
		++num;
	}

	return num;
}

// 剛体データの取得を行う関数
void clustscan(FILE *input)
{
	int i,j,nb,x,k;

	// 剛体数の取得
	fscanf(input,"%d",&x);
	prot.DOF = x;

	for(k=0;k<prot.DOF;++k)
	{
		fscanf(input,"%d",&x);
		clust[k].origin_atom_a = x;
	}

	for(k=0;k<prot.DOF;++k)
	{
		fscanf(input,"%d",&x);
		clust[k].terminal = x;
	}

	for(k=0;k<prot.DOF;++k)
	{
		fscanf(input,"%d",&x);
		clust[k].num_atom_clust = x;
	}

	for(k=0;k<prot.DOF;++k)
	{
		fscanf(input,"%d",&x);
		clust[k].num_branch = x;
	}

	for(k=0;k<prot.DOF;++k)
	{
		fscanf(input,"%d",&x);
		/******************************/
                /* ABI[k].hingmat = x;	      */
                /******************************/
	}

	for(k=0;k<prot.DOF;++k)
	{
		for(nb=0;nb<clust[k].num_branch;++nb)
		{
			fscanf(input,"%d",&x);
			clust[k].terminal_atom_a[nb] = x;
		}
	}

	for (i=0;i<prot.DOF;++i)
	{
		fscanf(input, "%d", &x);
		clust[i].nNumClutOfParent = x;
	}

	for (i=0;i<prot.DOF;++i)
	{
		for(nb=0;nb<clust[i].num_branch;++nb)
		{
			fscanf(input, "%d", &x);
			clust[i].nNumClutOfChild[nb] = x;
		}
	}

	for (i=0;i<prot.DOF;++i)
	{
		fscanf(input, "%d", &x);
		//		IndexOfABICycle[i] = x;
	}

	fscanf(input, "%d", &x);
	/*******************/
        /* MAP = x;	   */
        /*******************/

/*************************************************************************************/
/* 	if (MAP == ON)								     */
/* 	{									     */
/* 		fscanf(input, "%d", &x);					     */
/* 		for(i=0;;)							     */
/* 		{								     */
/* 			fscanf(input,"%d",&x);					     */
/* 			if (x == -1)						     */
/* 			{							     */
/* 				break;						     */
/* 			}							     */
/* 			++i;							     */
/* 			/\**************************************\/		     */
/*                         /\* nNumAtomPeptide_N[i] = x-1;	      *\/	     */
/*                         /\**************************************\/		     */
/* 		}								     */
/* 		/\********************************\/				     */
/*                 /\* prot.nNumPeptide = i;        *\/				     */
/*                 /\********************************\/				     */
/* 		for(i=0;;++i)							     */
/* 		{								     */
/* 			if (fscanf(input,"%d",&x) == EOF)			     */
/* 			{							     */
/* 				break;						     */
/* 			}							     */
/* 			else							     */
/* 			{							     */
/* 				/\**************************************\/	     */
/*                                 /\* nNumAtomPeptide_C[i] = x-1;	      *\/    */
/*                                 /\**************************************\/	     */
/* 			}							     */
/* 		}								     */
/* 	}									     */
/*************************************************************************************/
}

void pick_initial_velo(char *InpfilDVELO)
{
	int nNumClut;
	double ddihedang;

	FILE *input;

	// velo.in を開く
	if ((input=fopen(/*"velo.in"*/InpfilDVELO,"r")) == NULL)
	{
		printf("error:cannot open velo.in \n");
		exit(1);
	}

	// 初期速度の取得を行う
	for (nNumClut=0;nNumClut<prot.DOF;++nNumClut)
	{
		fscanf(input,"%lf",&ddihedang);
		clust[nNumClut].ddihedang[0] = ddihedang;
	}
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

void set_trans_Matrix(int nNumClt,
					  int nNumClutOrigBranch)
{
	int alpha;
	int alpha2;

	int nNumClutParent;

	nNumClutParent = clust[nNumClt].nNumClutOfParent-1;

//	set_delts_matrix();

	// 0番目の剛体のとき
	if (nNumClt == 0)
	{
		for (alpha=0;alpha</*3*/6;++alpha)
		{
			for (alpha2=0;alpha2</*3*/6;++alpha2)
			{
				clust[nNumClt].TransMatrix[0][alpha][alpha2]= 0.0;
			}
		}
	}
	else
	{
		sub_set_trans_Matrix(nNumClt, nNumClutParent);
	}

	set_pseduo_trans_matrix(nNumClt);
}

// 座標系変換行列の作成の補助を行う関数
void sub_set_trans_Matrix(int nNumClt,
	                      int nNumCltminousone)
{
	int alpha,alpha2,alpha3,i,j,k;
	int nNumAtomOfClut,nNumAtomOfClut2;
	int nNumAtomOfClutminousone;
	int num,num2;
	double Coord[3],Coord2[3];
	double RotatnNumtonNumMiOn[3][3];
	double mat2[3][3];

	FILE *outtest;

	for (alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
			RotatnNumtonNumMiOn[alpha][alpha2] = 0.0;
		}
	}

	for (alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
			for (alpha3=0;alpha3<3;++alpha3)
			{
//				// Rot_0_to_N * Rot_n-1_to_0
				RotatnNumtonNumMiOn[alpha][alpha2]
				+=  clust[nNumCltminousone].trans_A_to_CN[0][alpha][alpha3]
			       *clust[nNumClt].trans_A_to_CN[0][alpha2][alpha3];
//				RotatnNumtonNumMiOn[alpha][alpha2]
//				+=  clust[nNumCltminousone].trans_A_to_CN[0][alpha3][alpha]
//			       *clust[nNumClt].trans_A_to_CN[0][alpha3][alpha2];
			}
		}
	}

	nNumAtomOfClut = clust[nNumClt].origin_atom_a-1;
	nNumAtomOfClutminousone = clust[nNumCltminousone].origin_atom_a-1;

	// 座標の取得
	for(alpha=0;alpha<3;++alpha)
	{
		Coord[alpha]/*m*/
		=( prot.coord[nNumAtomOfClut][alpha]
		  -prot.coord[nNumAtomOfClutminousone][alpha])/*A*/
/////////////////////////////////////////////////////////////////////
//		  +prot.coord[nNumAtomOfClutminousone][alpha]/*A*/
/////////////////////////////////////////////////////////////////////
		/**1.0e-10*/;
	}

	// 座標の取得
	for(alpha=0;alpha<3;++alpha)
	{
		Coord2[alpha]=0.0;
	}

	for(alpha=0;alpha<3;++alpha)
	{
		for(alpha2=0;alpha2<3;++alpha2)
		{
//			Coord2[alpha]
//				+=clust[nNumClt].trans_A_to_CN[0][alpha][alpha2]
//			     *Coord[alpha2];
			Coord2[alpha]
				+= clust[nNumCltminousone].trans_A_to_CN[0][alpha][alpha2]
			      *Coord[alpha2];
		}
	}

//	// 原子番号の取得
//	nNumAtomOfClut = clust[nNumClt].origin_atom_a-1;
//	nNumAtomOfClutminousone = clust[nNumCltminousone].origin_atom_a-1;
//	nNumAtomOfClut2 = clust[nNumCltminoustwo].terminal_atom_a[0]-1;
//
//	num = nNumAtomOfClut - nNumAtomOfClut2;
//	num2 = nNumAtomOfClutminousone - nNumAtomOfClut2;
//
//	// 座標の取得
//	for(alpha=0;alpha<3;++alpha)
//	{
//		Coord2[alpha]/*m*/
//		=( clust[nNumCltminousone].xoord_clust/*[0]*/[num2][alpha]
//		  -clust[nNumCltminousone].xoord_clust/*[0]*/[num][alpha])/*A*/
//         /**1.0e-10*/;
//	}

	// 座標系変換行列の左上の作成
	for(alpha=0;alpha<3;++alpha)
	{
		for(alpha2=0;alpha2<3;++alpha2)
		{
			clust[nNumClt].TransMatrix[0][alpha][alpha2]
			         =RotatnNumtonNumMiOn[alpha][alpha2];
//			clust[nNumClt].TransMatrix[0][alpha][alpha2]
//			         =ident[alpha][alpha2];
		}
	}

	// 座標系変換行列の右下の作成
	for(alpha=3;alpha<6;++alpha)
	{
		for(alpha2=3;alpha2<6;++alpha2)
		{
			clust[nNumClt].TransMatrix[0][alpha][alpha2]
			     =RotatnNumtonNumMiOn[alpha-3][alpha2-3];
//			clust[nNumClt].TransMatrix[0][alpha][alpha2]
//			     =ident[alpha-3][alpha2-3];
		}
	}

	// 座標系変換行列の左下の作成
	for(alpha=3;alpha<6;++alpha)
	{
		for(alpha2=0;alpha2<3;++alpha2)
		{
			clust[nNumClt].TransMatrix[0][alpha][alpha2]=0.0;
		}
	}

//	// 座標系変換行列の右上の作成
//	clust[nNumClt].TransMatrix[0][0][3]= 0.0;
//	clust[nNumClt].TransMatrix[0][1][3]=-Coord2[2]/*m*/;
//	clust[nNumClt].TransMatrix[0][2][3]= Coord2[1]/*m*/;
//	clust[nNumClt].TransMatrix[0][0][4]= Coord2[2]/*m*/;
//	clust[nNumClt].TransMatrix[0][1][4]= 0.0;
//	clust[nNumClt].TransMatrix[0][2][4]=-Coord2[0]/*m*/;
//	clust[nNumClt].TransMatrix[0][0][5]=-Coord2[1]/*m*/;
//	clust[nNumClt].TransMatrix[0][1][5]= Coord2[0]/*m*/;
//	clust[nNumClt].TransMatrix[0][2][5]= 0.0;

//	// 座標系変換行列の右上の作成
//	mat2[0][0]= 0.0;
//	mat2[1][0]=-Coord2[2]/*m*/;
//	mat2[2][0]= Coord2[1]/*m*/;
//	mat2[0][1]= Coord2[2]/*m*/;
//	mat2[1][1]= 0.0;
//	mat2[2][1]=-Coord2[0]/*m*/;
//	mat2[0][2]=-Coord2[1]/*m*/;
//	mat2[1][2]= Coord2[0]/*m*/;
//	mat2[2][2]= 0.0;

	// 座標系変換行列の右上の作成
	mat2[0][0]= 0.0;
	mat2[0][1]=-Coord2[2]/*m*/;
	mat2[0][2]= Coord2[1]/*m*/;
	mat2[1][0]= Coord2[2]/*m*/;
	mat2[1][1]= 0.0;
	mat2[1][2]=-Coord2[0]/*m*/;
	mat2[2][0]=-Coord2[1]/*m*/;
	mat2[2][1]= Coord2[0]/*m*/;
	mat2[2][2]= 0.0;

	for (i=0;i<3;++i)
	{
		for (j=0;j<3;++j)
		{
			clust[nNumClt].TransMatrix[0][i][j+3] = 0.0;
		}
	}

	for (i=0;i<3;++i)
	{
		for (j=0;j<3;++j)
		{
			for (k=0;k<3;++k)
			{
				clust[nNumClt].TransMatrix[0][i][j+3] += mat2[i][k]*RotatnNumtonNumMiOn[k][j];
			}
		}
	}

// test
//	clust[2].TransMatrix[0][2][2] = -1.0;
//
}

void set_pseduo_trans_matrix(int nNumClut) {
  int alpha,alpha2;
  int nNumAtomOfClut,nNumAtomOfClutminousone;
  int nNumCltminousone;
  double Coord[3],Coord2[3];

  for (alpha=0;alpha<3;++alpha) {
    for (alpha2=0;alpha2<3;++alpha2) {
      clust[nNumClut].PsedoTransMatrix[alpha][alpha2] = clust[nNumClut].TransMatrix[0][alpha][alpha2];
    }
  }

  nNumCltminousone = clust[nNumClut].nNumClutOfParent-1;
  nNumAtomOfClut = clust[nNumClut].origin_atom_a-1;
  nNumAtomOfClutminousone = clust[nNumCltminousone].origin_atom_a-1;

  // 座標の取得
  for(alpha=0;alpha<3;++alpha) {
    Coord[alpha]/*m*/ =(  prot.coord[nNumAtomOfClut][alpha] - prot.coord[nNumAtomOfClutminousone][alpha])/*A*/;
  }

  // 座標の取得
  for(alpha=0;alpha<3;++alpha) {
    Coord2[alpha]=0.0;
  }

  for(alpha=0;alpha<3;++alpha) {
    for(alpha2=0;alpha2<3;++alpha2) {
      Coord2[alpha] += clust[nNumCltminousone].trans_A_to_CN[0][alpha][alpha2]*Coord[alpha2];
    }
  }

  clust[nNumClut].PsedoTransMatrix[0][3] = Coord2[0];
  clust[nNumClut].PsedoTransMatrix[1][3] = Coord2[1];
  clust[nNumClut].PsedoTransMatrix[2][3] = Coord2[2];

  for (alpha2=0;alpha2<3;++alpha2) {
    clust[nNumClut].PsedoTransMatrix[3][alpha2] = 0.0;
  }

  clust[nNumClut].PsedoTransMatrix[3][alpha2] = 1.0;
}

void set_delts_matrix(void)
{
	int i,j;

	for (i=0;i</*3*/4;++i)
	{
		for (j=0;j</*3*/4;++j)
		{
			delta_matrix[i][j] = 0.0;
		}
	}

	delta_matrix[0][1] = -1.0;
	delta_matrix[1][0] = 1.0;

}

