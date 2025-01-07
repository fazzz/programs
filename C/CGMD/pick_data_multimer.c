#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "gener.h"
#include "ABA.h" // #include "ABA_multimer.h"
#include "physics.h"
#include "force.h"
#include "MD.h"
#include "BD.h"

#include "ParmTop.h"
#include "EF.h"

int coordscan(FILE *input);
void seqscan(FILE *input);
void inertiascan(FILE *input);
void topscan(FILE *input);
void ParmTopscan(FILE *input);
void clustscan(FILE *input);
void mdinputscan(FILE *input);
void pick_initial_velo(void);
void ParmTopscan(FILE *input);

double **coord_rest/*[MAXA][3]*/;

int refcoordscan(FILE *input);

double pick_dihed_one_clust3(int nNumAom1,
			     int nNumAom2,
			     int nNumAom3,
			     int nNumAom4,
			     int nNumDihed);

void optimize_q_value(double omega,double *q_value);

// インプットデータの取得を行う関数
int pick_data(int pflag) {
  int i,j;
  char *name;
  char name_prot[10];
  FILE *input;
  char *line;
  size_t len=0;

  printf("%s \n", InpfilCOORD);
  printf("%s \n", InpfilCLUST);
  printf("%s \n", InpfilTOP);
  if (pflag != AMBERMODE)
    printf("%s \n", InpfilSEQ);
  printf("%s \n", InpfilMDI);
  if (restflag == ON) {
    printf("%s \n", InpfilCOORDREF);
  }

  if ((input=fopen(InpfilCOORD,"r")) == NULL){
    printf("error %s cannot open\n", InpfilCOORD);
    exit(1);
  }
  /******************************************/
  /* fscanf(input, "%s", name_prot);	  */
  /* prot.name_prot = name_prot;		  */
  /******************************************/
  getline(&line,&len,input);
  prot.num_atom = coordscan(input);
  fclose(input);

  if ((input=fopen(InpfilCLUST,"r")) == NULL){
    printf("error %s cannot open\n", InpfilCLUST);
    exit(1);
  }
  clustscan(input);
  fclose(input);
  
  if (pflag != AMBERMODE) {
    if ((input=fopen(InpfilSEQ,"r")) == NULL) {
      printf("error %s cannot open\n", InpfilSEQ);
      exit(1);
    }
    seqscan(input);
    fclose(input);
    
    if ((input=fopen("atom.dat","r")) == NULL){
      printf("error atom.data cannot open\n");
      exit(1);
    }
    inertiascan(input);
    fclose(input);
  }
	
  if ((input=fopen(InpfilTOP,"r")) == NULL) {
    printf("error %s cannot open\n", InpfilTOP);
    exit(1);
  }
  if (pflag != AMBERMODE) {
    topscan(input);
  }
  else {
    ParmTopscan(input);
  }
  fclose(input);
  
  if ((input=fopen(InpfilMDI,"r")) == NULL){
    printf("error %s cannot open\n", InpfilMDI);
    exit(1);
  }
  mdinputscan(input);
  fclose(input);
  if(restflag ==ON) {
    if ((input=fopen(InpfilCOORDREF,"r")) == NULL) {
      printf("error %s cannot open\n", InpfilCOORDREF);
      exit(1);
    }
    refcoordscan(input);
    fclose(input);
  }

}

// 座標データの取得を行う関数
int coordscan(FILE *input) {
  double x;
  int i,j;
  int flag = 0;
  int num = 0;
  int numatom;
  int numtype;// 0811

  fscanf(input,"%d",&i);
  numatom=i;

  prot.coord=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i)
    prot.coord[i]=(double *)gcemalloc(sizeof(double)*3);

  prot.velo=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i)
    prot.velo[i]=(double *)gcemalloc(sizeof(double)*3);

  /****************************************************************************************************************************/
  /* numtype = AP.NTYPES*(AP.NTYPES+1)/2;										      */
  /* prot.L_J_parm.A=(double *)gcemalloc(sizeof(double)*numtype/\*numatom*\/); //0811					      */
  /* prot.L_J_parm.B=(double *)gcemalloc(sizeof(double)*numtype/\*numatom*\/); //0811					      */
  /* prot.L_J_parm.atomtype=(/\*double*\/int *)gcemalloc(sizeof(/\*double*\/int)*numatom); //0811			      */
  /* prot.L_J_parm.atomtypeIndex=(/\*double*\/int **)gcemalloc(sizeof(/\*double*\/int *)*\/\*numatom*\/prot.num_atom); //0811 */
  /* for (i=0;i<numatom;++i) {												      */
  /*   prot.L_J_parm.atomtypeIndex[i]=(/\*double*\/int *)gcemalloc(sizeof(/\*double*\/int )*\/\*numatom*\/AP.NTYPES); //0811  */
  /* }															      */
  /****************************************************************************************************************************/
  
  for(i=0;i<numatom;++i) {
    for(j=0;j<3;++j) {
      if (fscanf(input,"%lf",&x) == EOF) {
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

int refcoordscan(FILE *input) {
  double x;
  int i,j,nNumDihed,nNumDihed_now;
  int flag = 0;
  int num = 0;
  int numatom;

  FILE *log_rest;

  fscanf(input,"%d",&i);
  numatom=i;

  coord_rest=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i)
    coord_rest[i]=(double *)gcemalloc(sizeof(double)*3);

  for(i=0;i<MAXA;++i) {
    for(j=0;j<3;++j) {
      if (fscanf(input,"%lf",&x) == EOF) {
	flag = 1;
	break;
      }
      coord_rest[i][j]=x;
    }
    if (flag == 1)
      break;
    ++num;
  }


  for (nNumDihed=0;nNumDihed<prot.nNumDihed_rest;++nNumDihed) {
    theta_ref[nNumDihed]/*rad*/ = pick_dihed_one_clust/*2*/(dihed_rest[nNumDihed][0]-1,							    dihed_rest[nNumDihed][1]-1,
							    dihed_rest[nNumDihed][2]-1,
							    dihed_rest[nNumDihed][3]-1,
							    nNumDihed);
  }

  if ((log_rest=fopen("log_rest.txt","a"))==NULL) {
    printf("error: cannot open log_rest.txt\n");
    exit(1);
  }

  fprintf(log_rest,"log_rest.txt","a");
  for (i=0;i<prot.nNumDihed_rest;++i)
    fprintf(log_rest,"%e \n",theta_ref[i]);
  fclose(log_rest);
  
  
}


// 残基データの取得を行う関数
void seqscan(FILE *input) {
  int x;
  int i,j,jj,n=0;
  int ID;
  int num_res;
  int atom_res[30][50];
  char *tt1;
  char *tt2;
  FILE *input2;

  for (i=0;i<30;++i){
    for (j=0;j<50;++j){
      atom_res[i][j] = 0;
    }
  }

  if ((input2=fopen("stdres.dat","r")) == NULL){
    printf("error stdres.dat cannot open\n");
    exit(1);
  }

  for(i=0;i<30;++i){
    for(j=0; j</*MAXA_RES*/50; ++j){
      x=fgetc(input2);
      if (isdigit(x)){
	atom_res[i][j]=x-'0';
      }
      else if (x == '\n'){
	break;
      }
    }
    atom_res[i][j] = '\n';
  }
  
  fclose(input2);
  
  for(i=0; i<MAXRES; ++i){
    if (fscanf(input,"%d",&x) == EOF)
      break;
    residue[i].res_ID=x;
  }

  prot.inumrs = i/*+1*/;

  for(j=0; j<prot.inumrs ; ++j)	{
    ID = residue[j].res_ID-1;
    for (jj=0;;++jj){
      if (atom_res[ID][jj] != 0) {
	prot.name_atom[n]=atom_res[ID][jj];
	++n;
      }
      else
	break;
    }
    residue[0].num_atom=jj-1;
  }
}

// 慣性データの取得を行う関数
void inertiascan(FILE *input) {
  int i,j,k;
  int atom_ID;
  double mass_atom[KATS];
  int n=0;
  int x;
  double y;

  for(i=1;i<10;++i){
    fscanf(input, "%d", &x);
    fscanf(input, "%lf", &y);
    mass_atom[x]=y;
  }

  for(i=0; i<prot.DOF; ++i){
    for(j=0; j<clust[i].num_atom_clust; ++j){
      atom_ID = prot.name_atom[n];
      clust[i].mass_clust[j] = mass_atom[atom_ID];
      ++n;
    }
  }
  
}

// 力場データの取得を行う関数
void topscan(FILE *input) {
  int i,kk,l,k, nNumDihedType, n, nn;
  double x;
  int y;

  int num_atomtype;
  int numtype;


//	for(k=1;k<prot.DOF;++k)
//	{
//		fscanf(input,"%lf",&x);
//		clust[k].f_p_clust.num_dihed_clust = x;
//	}

//	for(k=1;k<prot.DOF;++k)
//	{
//		fscanf(input,"%lf",&x);
//		clust[k].f_p_clust.num_dihed_clust_1 = x;
//	}

//	for(k=1;k<prot.DOF;++k)
//	{
//		fscanf(input,"%lf",&x);
//		clust[k].f_p_clust.num_dihed_clust_2 = x;
//	}

//	for(k=1;k<prot.DOF;++k)
//	{
//		for(kk=0;kk<clust[k].f_p_clust.num_dihed_clust_1;++kk)
//		{
//			fscanf(input,"%lf",&x);
//			clust[k].f_p_clust.num_ATOM_N[kk] = x;
//		}
//	}

//	for(k=1;k<prot.DOF;++k)
//	{
//		for(kk=0;kk<clust[k].f_p_clust.num_dihed_clust_2;++kk)
//		{
//			fscanf(input,"%lf",&x);
//			clust[k].f_p_clust.num_ATOM_N_1[kk] = x;
//		}
//	}

//	for(k=1;k<prot.DOF;++k)
//	{
//		for(kk=0;kk<clust[k].f_p_clust.num_dihed_clust;++kk)
//		{
//			fscanf(input,"%lf",&x);
//			clust[k].f_p_clust.theta_dihed[kk] = x/180.0*PI;
//		}
//	}

//	for(k=1;k<prot.DOF;++k)
//	{
//		for(kk=0;kk<clust[k].f_p_clust.num_dihed_clust;++kk)
//		{
//			fscanf(input,"%d",&y);
//			clust[k].f_p_clust.n_dihed[kk] = y;
//		}
//	}

//	for(k=1;k<prot.DOF;++k)
//	{
//		for(kk=0;kk<clust[k].f_p_clust.num_dihed_clust;++kk)
//		{
//			fscanf(input,"%lf",&x);
//			clust[k].f_p_clust.V_dihed[kk] = x;
//		}
//	}

// 2面角項
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*	for(k=1;k<prot.DOF;++k)
	{
		fscanf(input,"%lf",&x);
		clust[k].f_p_clust.theta_dihed[0] = x/180.0*PI;
	}

	for(k=1;k<prot.DOF;++k)
	{
		fscanf(input,"%d",&y);
		clust[k].f_p_clust.n_dihed[0] = y;
	}

	for(k=1;k<prot.DOF;++k)
	{
		fscanf(input,"%lf",&x);
		clust[k].f_p_clust.V_dihed[0] = x;
	}*/

  // 総二面角数の取得
  fscanf(input,"%d",&y);
  prot.nNumDihedALL = y;

  // 総二面角タイプ数の取得
  fscanf(input,"%d",&y);
  prot.nNumDihedType = y;

  // 二面角中の原子の番号の取得
  for (k=0;k<prot.nNumDihedALL;++k){
    // improper dihed


    for (i=0;i<4;++i){
      fscanf(input,"%lf",&x);
      atom_dihed_pair[k/*][*/*6+i] = abs(x)/3+1;
    }
    fscanf(input,"%lf",&x);
    atom_dihed_pair[k/*][*/*6+4] = x;
    fscanf(input,"%lf",&x);
    atom_dihed_pair[k/*][*/*6+5] = x;
  }

  for(nNumDihedType=0;nNumDihedType<prot.nNumDihedType;++nNumDihedType){
    fscanf(input,"%lf",&x);
    V_dihed[nNumDihedType] = x*2.0/*modified 22II10*/;
  }

  for(nNumDihedType=0;nNumDihedType<prot.nNumDihedType;++nNumDihedType){
    fscanf(input,"%d",&y);
    n_dihed[nNumDihedType] = y;
    //		fscanf(input,"%lf",&x);
    //		n_dihed[nNumDihedType] = x;
  }

  for(nNumDihedType=0;nNumDihedType<prot.nNumDihedType;++nNumDihedType){
    fscanf(input,"%lf",&x);
    theta_dihed[nNumDihedType] = x/*/180.0*PI*/;
  }

// 静電相互作用項
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

  for(k=0;k<prot.DOF;++k){
    for(kk=0;kk<clust[k].num_atom_clust;++kk){
      fscanf(input,"%lf",&x);
      clust[k].f_p_clust.e_f[kk] = x/*/18.2223*/;
    }
  }

// VDW互作用項
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

  fscanf(input,"%d",&y);
  num_atomtype = y;

  for(n=0;n<prot.num_atom;++n){
    fscanf(input,"%d",&y);
    prot.L_J_parm.atomtype[n] = y;
  }

  for(n=0;n<num_atomtype;++n){
    for(nn=0;nn<num_atomtype;++nn){
      fscanf(input,"%d",&y);
      prot.L_J_parm.atomtypeIndex[n][nn] = y;
    }
  }

  fscanf(input,"%d",&y);
  numtype = y;

  for(k=0;k<numtype;++k){
    fscanf(input,"%lf",&x);
    prot.L_J_parm.A[k] = x;
  }

  for(k=0;k<numtype;++k){
    fscanf(input,"%lf",&x);
    prot.L_J_parm.B[k] = x;
  }

// 非結合性相互作用の組
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

  for(k=0;k<prot.DOF;++k){
    for(kk=0;kk<clust[k].num_atom_clust;++kk){
      for(l=0;l<200;++l){
	fscanf(input,"%d",&y);
	if (y != 0){
	  clust[k].pairs.not_interacting[kk][l] = y;
	}
	else{
	  break;
	}
      }
    }
  }
// 1-4非結合性相互作用の組
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

  for(k=0;k<prot.DOF;++k){
    for(kk=0;kk<clust[k].num_atom_clust;++kk){
      for(l=0;l<20;++l){
	fscanf(input,"%d",&y);
	if (y != 0){
	  clust[k].o_f_pairs.o_f_not_interacting[kk][l] = y;
	}
	else{
	  break;
	}
      }
    }
  }
  //f	cndUmbSanscan();
}

// 剛体データの取得を行う関数
void clustscan(FILE *input){
  int i,j,nb,x,k;
  
  // 剛体数の取得
  fscanf(input,"%d",&x);
  prot.DOF = x;

  ABI=(struct Articulated_Body_Inertia *)gcemalloc(sizeof(struct Articulated_Body_Inertia)*prot.DOF);
  zzz=(struct Bias_force *)gcemalloc(sizeof(struct Bias_force)*prot.DOF);
  clust=(struct clustdata *)gcemalloc(sizeof(struct clustdata)*prot.DOF);
  for (i=0;i<prot.DOF;++i) {
    clust[i].xoord_clust=(double **)gcemalloc(sizeof(double *)*prot.num_atom);
    for (j=0;j<prot.num_atom;++j) {
      clust[i].xoord_clust[j]=(double *)gcemalloc(sizeof(double)*3);
    }
  }
  for (i=0;i<prot.DOF;++i) {
    clust[i].f_p_clust.e_f=(double *)gcemalloc(sizeof(double)*prot.num_atom);
    /************************************************************************************/
    /* clust[i].f_p_clust.V_dihed=(double *)gcemalloc(sizeof(double)*prot.numatom);     */
    /* clust[i].f_p_clust.theta_dihed=(double *)gcemalloc(sizeof(double)*prot.numatom); */
    /* clust[i].f_p_clust.n_dihed=(double *)gcemalloc(sizeof(double)*prot.numatom);     */
    /************************************************************************************/
  }

  old_dddihedang=(double *)gcemalloc(sizeof(double)*prot.DOF);

  potential_pro.p_elesta=(double *)gcemalloc(sizeof(double)*prot.num_atom);
  potential_pro.p_L_J=(double *)gcemalloc(sizeof(double)*prot.num_atom);
  potential_pro.p_1_4_elesta=(double *)gcemalloc(sizeof(double)*prot.num_atom);
  potential_pro.p_1_4_L_J=(double *)gcemalloc(sizeof(double)*prot.num_atom);
  potential_pro.p_at=(double *)gcemalloc(sizeof(double)*prot.num_atom);
  potential_pro.p_dihedc=(double *)gcemalloc(sizeof(double)*prot.DOF);
  potential_pro.p_rest=(double *)gcemalloc(sizeof(double)*prot.DOF);

  for(k=0;k<prot.DOF;++k){
    fscanf(input,"%d",&x);
    clust[k].origin_atom_a = x;
  }

  for(k=0;k<prot.DOF;++k){
    fscanf(input,"%d",&x);
    clust[k].terminal = x;
  }

  for(k=0;k<prot.DOF;++k){
    fscanf(input,"%d",&x);
    clust[k].num_atom_clust = x;
  }

  for(k=0;k<prot.DOF;++k){
    fscanf(input,"%d",&x);
    clust[k].num_branch = x;
  }

  for(k=0;k<prot.DOF;++k){
    fscanf(input,"%d",&x);
    ABI[k].hingmat = x;
  }

  for(k=0;k<prot.DOF;++k){
    for(nb=0;nb<clust[k].num_branch;++nb){
      fscanf(input,"%d",&x);
      clust[k].terminal_atom_a[nb] = x;
    }
  }

  for (i=0;i<prot.DOF;++i){
    fscanf(input, "%d", &x);
    clust[i].nNumClutOfParent = x;
  }

  for (i=0;i<prot.DOF;++i){
    for(nb=0;nb<clust[i].num_branch;++nb){
      fscanf(input, "%d", &x);
      clust[i].nNumClutOfChild[nb] = x;
    }
  }

  IndexOfABICycle=(int *)gcemalloc(sizeof(int)*prot.DOF);
  for (i=0;i<prot.DOF;++i){
    fscanf(input, "%d", &x);
    IndexOfABICycle[i] = x;
  }

  fscanf(input, "%d", &x);
  MAP = x;

  if (MAP == ON){
    fscanf(input, "%d", &x);
    for(i=0;;)	{
      fscanf(input,"%d",&x);
      if (x == -1){
	break;
      }
      ++i;
      nNumAtomPeptide_N[i] = x-1;
    }
    prot.nNumPeptide = i;
    for(i=0;;++i){
      if (fscanf(input,"%d",&x) == EOF){
	break;
      }
      else{
	nNumAtomPeptide_C[i] = x-1;
      }
    }
  }


  for (i=0;i<prot.DOF;++i) {
    clust[i].f_c.f_elesta=(double **)gcemalloc(sizeof(double *)*clust[i].num_atom_clust);
    clust[i].f_c.f_L_J=(double **)gcemalloc(sizeof(double *)*clust[i].num_atom_clust);
    clust[i].f_c.f_1_4_elesta=(double **)gcemalloc(sizeof(double *)*clust[i].num_atom_clust);
    clust[i].f_c.f_1_4_L_J=(double **)gcemalloc(sizeof(double *)*clust[i].num_atom_clust);
    clust[i].f_c.f_atoms=(double **)gcemalloc(sizeof(double *)*clust[i].num_atom_clust);
    for (j=0;j<clust[i].num_atom_clust;++j) {
      clust[i].f_c.f_elesta[j]=(double *)gcemalloc(sizeof(double)*3);
      clust[i].f_c.f_L_J[j]=(double *)gcemalloc(sizeof(double)*3);
      clust[i].f_c.f_1_4_elesta[j]=(double *)gcemalloc(sizeof(double)*3);
      clust[i].f_c.f_1_4_L_J[j]=(double *)gcemalloc(sizeof(double)*3);
      clust[i].f_c.f_atoms[j]=(double *)gcemalloc(sizeof(double)*3);
    }
    clust[i].pairs.not_interacting=(int/*double*/ **)gcemalloc(sizeof(int/*double*/ *)*prot.num_atom);
    clust[i].o_f_pairs.o_f_not_interacting=(int/*double*/ **)gcemalloc(sizeof(int/*double*/ *)*prot.num_atom);
    for (j=0;j<50;++j) {
      clust[i].pairs.not_interacting[j]=(int/*double*/ *)gcemalloc(sizeof(int/*double*/)*50);
      clust[i].o_f_pairs.o_f_not_interacting[j]=(int/*double*/ *)gcemalloc(sizeof(int/*double*/)*50);
    }
  }


/*	for (i=prot.DOF-1;i<0;--i)
	{
		if (i==prot.DOF-1)
		{
			Next = prot.DOF-1;
			IndexforABICycle[nNumClut]=Next;
		}
		if (clust[i].num_branch > 1)
		{
			Next = prot.DOF-1;
			IndexforABICycle[nNumClut]=Next;
		}
		if (i==prot.DOF-1)
		{
			Next = prot.DOF-1;
			IndexforABICycle[nNumClut]=Next;
		}
	}*/
}

// 計算条件の取得を行う関数
void mdinputscan(FILE *input) {
  int i,j,l,k;
  double x;
  int y;
  double omega;
  int m[MAXDOF];
  int temp[MAXDOF];
  FILE *log_rest;

  // 統計集団
  fscanf(input,"%d",&y);
  if (y == 0 || y == 1){
    MODE = y;
  }
  else{
    printf("error; ANNSANNBLE MODE SET KEY MUST BE 0 OR 1 !!!\n");
    exit(1);
  }

  // 速度の設定
  fscanf(input,"%d",&y);
  if (y == 0 || y == 1)	{
    MODEV = y;
  }
  else{
    printf("error; VELOCITY SET KEY MUST BE 0 OR 1 !!! \n");
    exit(1);
  }

//	if (MODE == NVT && MODEV == OFF)
//	{
//		printf("error; conflict option ! \n");
//		exit(1);
//	}

  // 温度
  fscanf(input,"%lf",&x);
  T_Kelvin = x;

  // 初期温度
  fscanf(input,"%lf",&x);
  T_Kelvin_Initial = x;

  // タイムステップの幅
  fscanf(input,"%lf",&x);
  deltat = x/**1e-12*/;

  // タイムステップ数
  fscanf(input,"%d",&y);
  TIME_LIMIT = y;

  // カットオフ半径 // iteration num
  fscanf(input,"%d",&y);
  NUM_IETRATION=y;
  //	Cut_Off = x;

  
  // 体積 // GearMode flag
  fscanf(input,"%d",&y);
  GearMODE = y;
  //  V_system = x;

  // ステップ_1(構造、速度)
  fscanf(input,"%d",&y);
  out_put_steps = y;

  // ステップ_2(エネルギー、 etc)
  fscanf(input,"%d",&y);
  out_put_steps_thomo = y;

  // 再スタート
  fscanf(input,"%d",&y);
  if (y == 0 || y == 1){
    Reatat = y;
  }
  else {
    printf("error; VELOCITY SET KEY MUST BE 0 OR 1 !!! \n");
    exit(1);
  }

  // LD
  fscanf(input,"%d",&y);
  DYNMODE = y;

  //	if (MODE == NVT)
  //	{
  for(i=0;i<6;++i){
    s_NVT_correct[i] = 0.0;			
  }
  for(i=0;i<5;++i){
    xi_NVT_correct[i] = 0.0;			
  }
  // vv
  fscanf(input,"%lf",&x);
//		xi_NVT = x;
//		s_NVT_correct[1] = xi_NVT*deltat;
  s_NVT_correct[1] = x*deltat;
  dot_s_NVT = s_NVT_correct[1]/deltat;
//		xi_NVT_correct[1] = 0.0;
  
  // vvop
  fscanf(input,"%lf",&x);
  //  tau_NVT = x;
  omega = x;

  // vvop
  fscanf(input,"%lf",&x);
  s_NVT=x;
  s_NVT_correct[0] = s_NVT;
  xi_NVT = s_NVT_correct[1]/deltat/s_NVT;
//		dot_xi_NVT = s_NVT_correct[1]/(s_NVT_correct[0]*s_NVT_correct[0]);
  xi_NVT_correct[0] = s_NVT_correct[1]/deltat/s_NVT;
  xi_NVT_correct[1] = s_NVT_correct[1]/(s_NVT_correct[0]*s_NVT_correct[0]);
  /*************************************************************************************************************************/
  /* q_NVT = tau_NVT*deltat*tau_NVT*deltat*(prot.DOF-1)*k_B*T_Kelvin*2.3889e-4/\*Kcal/J*\/\*6.022142e23/\*\/mol*\/;	 */
  /*************************************************************************************************************************/
  q_NVT = tau_NVT*deltat*tau_NVT*deltat*(prot.DOF-1)*k_B_kcm*T_Kelvin*4.18407*100.0;
  //	}
  
  s_NVT_correct[0]=1.0;
  for (i=1;i<6;++i) {
    s_NVT_correct[i]=0.0;
  }


  dhstopflag=0;
  nbstopflag=0;
  fscanf(input,"%d",&dhstopflag);
  fscanf(input,"%d",&nbstopflag);
  veloutflag=0;
  fscanf(input,"%d",&veloutflag);

  restflag=0;
  if(fscanf(input,"%d",&restflag)== EOF) {
    restflag=0;
  }
  if (restflag!=3) {
    fscanf(input,"%d",&prot.nNumDihed_rest);
  }
  if (restflag==1) {
    for (i=0;i<prot.nNumDihed_rest;++i) {
      fscanf(input,"%d",&y);
      for (j=0;j<prot.nNumDihedALL;++j) {
	if (y==atom_dihed_pair[j/*][*/*6+5]-1) {
	  for (k=0;k<4;++k) {
	    dihed_rest[i][k]=atom_dihed_pair[j/*][*/*6+k];
	  }
	  dihed_rest[i][4]=atom_dihed_pair[j/*][*/*6+5];
	  break;
	}
      }
    }
  }
  else if (restflag==2) {
    for (i=0;i<prot.nNumDihed_rest;++i) {
      fscanf(input,"%d",&m[i]);
    }
    j=0;
    for (i=0;i<prot.DOF;++i) {
      if (i!=m[i]-1) {
	dihed_rest[j][4]=i+2;
	++j;
      }
    }
    prot.nNumDihed_rest=j;
    for (i=0;i<prot.nNumDihed_rest;++i) {
      for (j=0;j<prot.nNumDihedALL;++j) {
	if (dihed_rest[i][4]==atom_dihed_pair[j/*][*/*6+5]) {
	  for (k=0;k<4;++k) {
	    dihed_rest[i][k]=atom_dihed_pair[j/*][*/*6+k];
	  }
	  break;
	}
      }
    }
    restflag=ON;
  }
  else if (restflag==3) {
    j=0;
    for (i=0;i<prot.DOF;++i) {
      if (      strncmp(AP.IGRAPH[clust[i].origin_atom_a-1],"CA\0",3) == 0 
		|| strncmp(AP.IGRAPH[clust[i].origin_atom_a-1],"C\0\0",3) == 0
		|| strncmp(AP.IGRAPH[clust[i].origin_atom_a-1],"N\0\0",3) == 0
		) {
	temp[j] = i;
	++j;
      }
    }
    prot.nNumDihed_rest = j;
    for (i=0;i<prot.nNumDihed_rest;++i) {
      y=temp[i];
      for (j=0;j<prot.nNumDihedALL;++j) {
	if (y==atom_dihed_pair[j/*][*/*6+5]-1) {
	  for (k=0;k<4;++k) {
	    dihed_rest[i][k]=atom_dihed_pair[j/*][*/*6+k];
	  }
	  dihed_rest[i][4]=atom_dihed_pair[j/*][*/*6+5];
	  break;
	}
      }
    }
    restflag=ON;
    if ((log_rest=fopen("log_rest.txt","w"))==NULL) {
      printf("error: cannot open log_rest.txt\n");
      exit(1);
    }
    for (i=0;i<prot.nNumDihed_rest;++i) {
      for (j=0;j<4;++j) {
	fprintf(log_rest,"%d  ",dihed_rest[i][j]);
      }
      fprintf(log_rest,"\n ");
    }
    fclose(log_rest);
  }
  fscanf(input,"%lf",&x);
  for (i=0;i<prot.nNumDihed_rest;++i) {
    V_rest[i]=x;
  }

  /********************/
  /* TermMoveMode=ON; */
  /********************/
  fscanf(input,"%d",&TermMoveMode);
  TermMoveMode2=0;
  fscanf(input,"%d",&TermMoveMode2);
  fscanf(input,"%d",&MASSIVEOUT);
  fscanf(input,"%d",&PEPCAACCMODE);
    

  /*****************************/
  /* for (i=0;i<6;++i)	       */
  /*   correct_Term[i][1]=1.0; */
  /*****************************/

  //	if (restflag==ON) {
  //	  fscanf(input,"%d",&prot.nNumDihed_rest);
  //	  for (i=0;i<prot.nNumDihed_rest;++i) {
  //	    fscanf(input,"%d %d %d %d ",&dihed_rest[i][0],&dihed_rest[i][1],&dihed_rest[i][2],&dihed_rest[i][3]);
  //	  }
  //	  for (i=0;i<prot.nNumDihed_rest;++i) {
  //	    fscanf(input,"%lf",&x);
  //    V_rest[i]=x;
  //  }
  //	}
  //
  //	for (i=0;i<prot.nNumDihed_rest;++i) {
  //  l = 1;
  //  if (dihed_rest[i][1] < dihed_rest[i][2]) {
  //    l = 2;
  //  }
  //  for (j=0;j<prot.DOF;++j){
  //    if (dihed_rest[i][l]==clust[j].origin_atom_a){
  //      dihed_rest[i][4] = j+1;
  //    }	
  //  }
  //	}

  optimize_q_value(omega,&q_NVT);
}

// 初期速度の取得を行う関数
void pick_initial_velo(void)
{
  int i,j;
	int nNumClut;
	double ddihedang;

	FILE *input;

	// velo.in を開く
	if ((input=fopen(/*"velo.in"*/InpfilDVELO,"r")) == NULL)
	{
		printf("error:cannot open velo.in \n");
		exit(1);
	}

	if (TermMoveMode2==4 || TermMoveMode2==12){
	  for (i=0;i<6;++i){
	    fscanf(input,"%lf",&vel_Term[i]);
	  }
	}

	// 初期速度の取得を行う
	if (TermMoveMode2==4 || TermMoveMode2==12){
	  for (nNumClut=1;nNumClut<prot.DOF;++nNumClut)   {
	    if (ABI[nNumClut].hingmat!=6) {
	      fscanf(input,"%lf",&ddihedang);
	      clust[nNumClut].ddihedang[0] = ddihedang;
	    }
	    else {
	      for (j=0;j<6;++j) {
		fscanf(input,"%lf",&ddihedang);
		clust[nNumClut].ddihedang_six[j] = ddihedang;
		clust[nNumClut].correct_dihedang_six[j][1] = ddihedang*deltat;
	      }
	    }
	  }
	}
	else {
	  for (nNumClut=0;nNumClut<prot.DOF;++nNumClut)   {
	    if (ABI[nNumClut].hingmat!=6) {
	      fscanf(input,"%lf",&ddihedang);
	      clust[nNumClut].ddihedang[0] = ddihedang;
	    }
	    else {
	      for (j=0;j<6;++j) {
		fscanf(input,"%lf",&ddihedang);
		clust[nNumClut].ddihedang_six[j] = ddihedang;
		clust[nNumClut].correct_dihedang_six[j][1] = ddihedang*deltat;
	      }
	    }
	  }
	}

	if (MODE==NVT) {
	  fscanf(input,"%lf",&s_NVT_correct[0]);
	  fscanf(input,"%lf",&s_NVT_correct[1]);
	  s_NVT=s_NVT_correct[0];
	  dot_s_NVT = s_NVT_correct[1]/deltat;
	  xi_NVT = s_NVT_correct[1]/deltat/s_NVT;
	  xi_NVT_correct[0] = s_NVT_correct[1]/deltat/s_NVT;
	  xi_NVT_correct[1] = s_NVT_correct[1]/(s_NVT_correct[0]*s_NVT_correct[0]);

	}
}

// 力場データの取得を行う関数
void ParmTopscan(FILE *input)
{
  int i,j,l,m,ll,k, nNumDihedType, n, nn,na,p;
  double x;
  int y;
  int flag,flag2;
  //  double d[100];
  int not_1_4[40];

  int num_atomtype;
  int numtype;

  int sum;
  FILE *debug;
  FILE *out;

  int inpindexH[100],inpindexA[100];
  int inpnumH=0,inpnumA=0;


  readParmtop(input);

  // 0811
  numtype = AP.NTYPES*(AP.NTYPES+1)/2;
  prot.L_J_parm.A=(double *)gcemalloc(sizeof(double)*numtype/*numatom*/); //0811
  prot.L_J_parm.B=(double *)gcemalloc(sizeof(double)*numtype/*numatom*/); //0811
  prot.L_J_parm.atomtype=(/*double*/int *)gcemalloc(sizeof(/*double*/int)*prot.num_atom); //0811
  prot.L_J_parm.atomtypeIndex=(/*double*/int **)gcemalloc(sizeof(/*double*/int *)*/*numatom*/prot.num_atom); //0811
  for (i=0;i<prot.num_atom;++i) {
    prot.L_J_parm.atomtypeIndex[i]=(/*double*/int *)gcemalloc(sizeof(/*double*/int )*/*numatom*/AP.NTYPES); //0811
  }
  // 0811

  prot.num_atom = AP.NATOM;

  // 総二面角数の取得
  prot.nNumDihedALL = AP.NPHIH+AP.MPHIA;

  atom_dihed_pair=(int *)malloc(sizeof(int)*prot.nNumDihedALL*6);

  for (k=0;k<prot.nNumDihedALL;++k) {
    atom_dihed_pair[k/*][*/*6+5] = 0;
  }

  // 総二面角タイプ数の取得
  prot.nNumDihedType = AP.NPTRA;

  // 二面角中の原子の番号の取得
  for (k=0;k<AP.NPHIH;++k) {
    // check for improper dihed
    flag=OFF;
    for (i=0;i<AP.NBONH;++i) {
      if (((AP.BH[i][0] == abs(AP.PH[k][1]) && AP.BH[i][1] == abs(AP.PH[k][0])) || (AP.BH[i][0] == abs(AP.PH[k][0]) && AP.BH[i][1] == abs(AP.PH[k][1])))) {
	flag = ON;
	break;
      }
    }
    if (flag==OFF) {
      for (i=0;i<AP.NBONA;++i) {
	if (((AP.BA[i][0] == abs(AP.PH[k][1]) && AP.BA[i][1] == abs(AP.PH[k][0])) || (AP.BA[i][0] == abs(AP.PH[k][0]) && AP.BA[i][1] == abs(AP.PH[k][1])))) {
	  flag = ON;
	  break;
	}
      }
    }
    if (flag==ON) {
      flag=OFF;
      for (i=0;i<AP.NBONH;++i) {
	if (((AP.BH[i][0] == abs(AP.PH[k][2]) && AP.BH[i][1] == abs(AP.PH[k][3])) || (AP.BH[i][0] == abs(AP.PH[k][3]) && AP.BH[i][1] == abs(AP.PH[k][2]))) ) {
	  flag = ON;
	  break;
	}
      }
      if (flag==OFF ) {
	for (i=0;i<AP.NBONA;++i) {
	  if (((AP.BA[i][0] == abs(AP.PH[k][2]) && AP.BA[i][1] == abs(AP.PH[k][3])) || (AP.BA[i][0] == abs(AP.PH[k][3]) && AP.BA[i][1] == abs(AP.PH[k][2]))) ) {
	    flag = ON;
	    break;
	  }
	}
      }
    }
    if (flag==OFF) {
      inpindexH[inpnumH]=k;
      ++inpnumH;
    }
    //    else {

    for (i=0;i<4;++i) {
      atom_dihed_pair[k/*][*/*6+i] = abs(AP.PH[k][i])/3+1;
    }
    atom_dihed_pair[k/*][*/*6+4] = AP.PH[k][i];
    l = 1;
    ll = 2;
    if (atom_dihed_pair[k/*][*/*6+1] > atom_dihed_pair[k/*][*/*6+2]) {
      l = 2;
      ll = 1;
    }
    //    else {
    for (i=0;i<prot.DOF;++i){
      p = clust[i].nNumClutOfParent-1;
      if (atom_dihed_pair[k/*][*/*6+l]==clust[p].terminal_atom_a[0] && atom_dihed_pair[k/*][*/*6+ll]==clust[i].origin_atom_a) {
	  /******************************************************************/
          /* flag=OFF;							    */
	  /* flag2=OFF;							    */
	  /* 								    */
	  /* if (atom_dihed_pair[k][0] < atom_dihed_pair[k][1] ) {	    */
	  /*   for (j=0;j<AP.NBONH;++j) {				    */
	  /*     if (AP.BH[j][1] == (atom_dihed_pair[k][1]-1)*3) {	    */
	  /* 	if ((atom_dihed_pair[k][0]-1)*3 == AP.BH[j][0]) {	    */
	  /* 	  flag = ON;						    */
	  /* 	}							    */
	  /*     }							    */
	  /*   }							    */
	  /*   for (j=0;j<AP.NBONA;++j) {				    */
	  /*     if (AP.BA[j][1] == (atom_dihed_pair[k][1]-1)*3) {	    */
	  /* 	if ((atom_dihed_pair[k][0]-1)*3 == AP.BA[j][0]) {	    */
	  /* 	  flag = ON;						    */
	  /* 	}							    */
	  /*     }							    */
	  /*   }							    */
	  /* }								    */
	  /* else {							    */
	  /*   for (j=0;j<AP.NBONH;++j) {				    */
	  /*     if (AP.BH[j][0] == (atom_dihed_pair[k][1]-1)*3) {	    */
	  /* 	if ((atom_dihed_pair[k][0]-1)*3 == AP.BH[j][1]) {	    */
	  /* 	  flag = ON;						    */
	  /* 	}							    */
	  /*     }							    */
	  /*   }							    */
	  /*   for (j=0;j<AP.NBONA;++j) {				    */
	  /*     if (AP.BA[j][0] == (atom_dihed_pair[k][1]-1)*3) {	    */
	  /* 	if ((atom_dihed_pair[k][0]-1)*3 == AP.BA[j][1]) {	    */
	  /* 	  flag = ON;						    */
	  /* 	}							    */
	  /*     }							    */
	  /*   }							    */
	  /* }								    */
	  /* 								    */
	  /* if (atom_dihed_pair[k][2] < atom_dihed_pair[k][3] ) {	    */
	  /*   for (j=0;j<AP.NBONH;++j) {				    */
	  /*     if (AP.BH[j][0] == (atom_dihed_pair[k][2]-1)*3) {	    */
	  /* 	if ((atom_dihed_pair[k][3]-1)*3 == AP.BH[j][1]) {	    */
	  /* 	  flag2 = ON;						    */
	  /* 	}							    */
	  /*     }							    */
	  /*   }							    */
	  /*   for (j=0;j<AP.NBONA;++j) {				    */
	  /*     if (AP.BA[j][0] == (atom_dihed_pair[k][2]-1)*3) {	    */
	  /* 	if ((atom_dihed_pair[k][3]-1)*3 == AP.BA[j][1]) {	    */
	  /* 	  flag2 = ON;						    */
	  /* 	}							    */
	  /*     }							    */
	  /*   }							    */
	  /* }								    */
	  /* else {							    */
	  /*   for (j=0;j<AP.NBONH;++j) {				    */
	  /*     if (AP.BH[j][1] == (atom_dihed_pair[k][2]-1)*3) {	    */
	  /* 	if ((atom_dihed_pair[k][3]-1)*3 == AP.BH[j][0]) {	    */
	  /* 	  flag2 = ON;						    */
	  /* 	}							    */
	  /*     }							    */
	  /*   }							    */
	  /*   for (j=0;j<AP.NBONA;++j) {				    */
	  /*     if (AP.BA[j][1] == (atom_dihed_pair[k][1]-1)*3) {	    */
	  /* 	if ((atom_dihed_pair[k][3]-1)*3 == AP.BA[j][0]) {	    */
	  /* 	  flag2 = ON;						    */
	  /* 	}							    */
	  /*     }							    */
	  /*   }							    */
	  /* }								    */
	  /* 								    */
	  /* if (flag == ON && flag2 == ON)				    */
          /******************************************************************/
	atom_dihed_pair[k/*][*/*6+5] = i+1;
      }	
    }
    //  }
  }
  //  }
  for (k=0;k<AP.NPHIA;++k) {
    // check for improper dihed
    flag=OFF;
    for (i=0;i<AP.NBONH;++i) {
      if (   ((AP.BH[i][0] == abs(AP.PA[k][1]) && AP.BH[i][1] == abs(AP.PA[k][0])) || (AP.BH[i][0] == abs(AP.PA[k][0]) && AP.BH[i][1] == abs(AP.PA[k][1])))) {
	flag = ON;
	break;
      }
    }
    if (flag==OFF) {
      for (i=0;i<AP.NBONA;++i) {
	if (   ((AP.BA[i][0] == abs(AP.PA[k][1]) && AP.BA[i][1] == abs(AP.PA[k][0])) || (AP.BA[i][0] == abs(AP.PA[k][0]) && AP.BA[i][1] == abs(AP.PA[k][1])))) {
	  flag = ON;
	  break;
	}
      }
    }
    if (flag==ON) {
      flag=OFF;
      for (i=0;i<AP.NBONH;++i) {
	if (((AP.BH[i][0] == abs(AP.PA[k][2]) && AP.BH[i][1] == abs(AP.PA[k][3])) || (AP.BH[i][0] == abs(AP.PA[k][3]) && AP.BH[i][1] == abs(AP.PA[k][2]))) ) {
	  flag = ON;
	  break;
	}
      }
      if (flag==OFF ) {
	for (i=0;i<AP.NBONA;++i) {
	  if (((AP.BA[i][0] == abs(AP.PA[k][2]) && AP.BA[i][1] == abs(AP.PA[k][3])) || (AP.BA[i][0] == abs(AP.PA[k][3]) && AP.BA[i][1] == abs(AP.PA[k][2]))) ) {
	    flag = ON;
	    break;
	  }
	}
      }
    }
    if (flag==OFF) {
      inpindexA[inpnumA]=k;
      ++inpnumA;
    }
    //    else {

    for (i=0;i<4;++i) {
      atom_dihed_pair[(k+AP.NPHIH)/*][*/*6+i] = abs(AP.PA[k][i])/3+1;
    }
    atom_dihed_pair[(k+AP.NPHIH)/*][*/*6+4] = AP.PA[k][4];
    l = 1;
    ll = 2;
    if (atom_dihed_pair[(k+AP.NPHIH)/*][*/*6+1] > atom_dihed_pair[(k+AP.NPHIH)/*][*/*6+2]) {
      l = 2;
      ll = 1;
    }
    //    else {
    for (i=0;i<prot.DOF;++i){
      p = clust[i].nNumClutOfParent-1;
      if (atom_dihed_pair[(k+AP.NPHIH)/*][*/*6+l]==clust[p].terminal_atom_a[0] && atom_dihed_pair[(k+AP.NPHIH)/*][*/*6+ll]==clust[i].origin_atom_a){
	  /**********************************************************************************/
          /* flag=OFF;									    */
	  /* flag2=OFF;									    */
	  /* 										    */
	  /* if (atom_dihed_pair[k+AP.NPHIH][0] < atom_dihed_pair[k+AP.NPHIH][1] ) {	    */
	  /*   for (j=0;j<AP.NBONH;++j) {						    */
	  /*     if (AP.BH[j][1] == (atom_dihed_pair[k+AP.NPHIH][1]-1)*3) {		    */
	  /* 	if ((atom_dihed_pair[k+AP.NPHIH][0]-1)*3 == AP.BH[j][0]) {		    */
	  /* 	  flag = ON;								    */
	  /* 	}									    */
	  /*     }									    */
	  /*   }									    */
	  /*   for (j=0;j<AP.NBONA;++j) {						    */
	  /*     if (AP.BA[j][1] == (atom_dihed_pair[k+AP.NPHIH][1]-1)*3) {		    */
	  /* 	if ((atom_dihed_pair[k+AP.NPHIH][0]-1)*3 == AP.BA[j][0]) {		    */
	  /* 	  flag = ON;								    */
	  /* 	}									    */
	  /*     }									    */
	  /*   }									    */
	  /* }										    */
	  /* else {									    */
	  /*   for (j=0;j<AP.NBONH;++j) {						    */
	  /*     if (AP.BH[j][0] == (atom_dihed_pair[k+AP.NPHIH][1]-1)*3) {		    */
	  /* 	if ((atom_dihed_pair[k+AP.NPHIH][0]-1)*3 == AP.BH[j][1]) {		    */
	  /* 	  flag = ON;								    */
	  /* 	}									    */
	  /*     }									    */
	  /*   }									    */
	  /*   for (j=0;j<AP.NBONA;++j) {						    */
	  /*     if (AP.BA[j][0] == (atom_dihed_pair[k+AP.NPHIH][1]-1)*3) {		    */
	  /* 	if ((atom_dihed_pair[k+AP.NPHIH][0]-1)*3 == AP.BA[j][1]) {		    */
	  /* 	  flag = ON;								    */
	  /* 	}									    */
	  /*     }									    */
	  /*   }									    */
	  /* }										    */
	  /* 										    */
	  /* if (atom_dihed_pair[k+AP.NPHIH][2] < atom_dihed_pair[k+AP.NPHIH][3] ) {	    */
	  /*   for (j=0;j<AP.NBONH;++j) {						    */
	  /*     if (AP.BH[j][0] == (atom_dihed_pair[k+AP.NPHIH][2]-1)*3) {		    */
	  /* 	if ((atom_dihed_pair[k+AP.NPHIH][3]-1)*3 == AP.BH[j][1]) {		    */
	  /* 	  flag2 = ON;								    */
	  /* 	}									    */
	  /*     }									    */
	  /*   }									    */
	  /*   for (j=0;j<AP.NBONA;++j) {						    */
	  /*     if (AP.BA[j][0] == (atom_dihed_pair[k+AP.NPHIH][2]-1)*3) {		    */
	  /* 	if ((atom_dihed_pair[k+AP.NPHIH][3]-1)*3 == AP.BA[j][1]) {		    */
	  /* 	  flag2 = ON;								    */
	  /* 	}									    */
	  /*     }									    */
	  /*   }									    */
	  /* }										    */
	  /* else {									    */
	  /*   for (j=0;j<AP.NBONH;++j) {						    */
	  /*     if (AP.BH[j][1] == (atom_dihed_pair[k+AP.NPHIH][2]-1)*3) {		    */
	  /* 	if ((atom_dihed_pair[k+AP.NPHIH][3]-1)*3 == AP.BH[j][0]) {		    */
	  /* 	  flag2 = ON;								    */
	  /* 	}									    */
	  /*     }									    */
	  /*   }									    */
	  /*   for (j=0;j<AP.NBONA;++j) {						    */
	  /*     if (AP.BA[j][1] == (atom_dihed_pair[k+AP.NPHIH][1]-1)*3) {		    */
	  /* 	if ((atom_dihed_pair[k+AP.NPHIH][3]-1)*3 == AP.BA[j][0]) {		    */
	  /* 	  flag2 = ON;								    */
	  /* 	}									    */
	  /*     }									    */
	  /*   }									    */
	  /* }										    */
	  /* 										    */
	  /* 										    */
	  /* if (flag == ON && flag2 == ON)						    */
          /**********************************************************************************/
	atom_dihed_pair[(k+AP.NPHIH)/*][*/*6+5] = i+1;
      }	
    }
    //    }
  }
  //  }
  /*************************************************/
  /* out=fopen("pair.txt","w");			   */
  /* for (i=0;i<prot.nNumDihedALL;++i) {	   */
  /*   for (j=0;j<6;++j) {			   */
  /*     fprintf(out,"%d ",atom_dihed_pair[i][j]); */
  /*   }					   */
  /*   fprintf(out,"\n");			   */
  /* }						   */
  /* fclose(out);				   */
  /*************************************************/

  for(nNumDihedType=0;nNumDihedType<prot.nNumDihedType;++nNumDihedType){
    n_dihed[nNumDihedType] = AP.PN[nNumDihedType];
  }

  for(nNumDihedType=0;nNumDihedType<prot.nNumDihedType;++nNumDihedType){
    theta_dihed[nNumDihedType] = AP.PHASE[nNumDihedType]/*/180.0*PI*/;
  }

  //  for(i=0;i<prot.nNumDihedType;++i){
  //    V_dihed[i] = AP.RK[i]*2.0/*modified 22II10*/;
  //  }


// 静電相互作用項
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

  na=0;
  for(k=0;k<prot.DOF;++k){
    for(kk=0;kk<clust[k].num_atom_clust;++kk){
      clust[k].f_p_clust.e_f[kk] = AP.CHRG[na]/*18.2223*/;
      ++na;
    }
  }

// VDW互作用項
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

  num_atomtype = AP.NTYPES;

  for(n=0;n<prot.num_atom;++n){
    prot.L_J_parm.atomtype[n] = AP.IAC[n];
  }

  for(n=0;n<prot.num_atom;++n){
    for(nn=0;nn<num_atomtype;++nn){
      prot.L_J_parm.atomtypeIndex[n][nn] = AP.ICO[n*num_atomtype+nn];
    }
  }

  numtype = AP.NTYPES*(AP.NTYPES+1)/2;

  for(k=0;k<numtype;++k){
    prot.L_J_parm.A[k] = AP.CN1[k];
  }

  for(k=0;k<numtype;++k){
    prot.L_J_parm.B[k] = AP.CN2[k];
  }

// 非結合性相互作用の組
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

  na=0;
  sum=0;
  for(k=0;k<prot.DOF;++k){
    for(kk=0;kk<clust[k].num_atom_clust;++kk){
      for(l=0;l<AP.NUMEX[na];++l){
	clust[k].pairs.not_interacting[kk][l] = AP.NATEX[sum+l];
      }
      clust[k].pairs.not_interacting[kk][l] = 0;
      ++na;
      sum+=l;
    }
  }
// 1-4非結合性相互作用の組
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*  na=0;
  for(k=0;k<prot.DOF;++k){
    for(kk=0;kk<clust[k].num_atom_clust;++kk){
      ll=0;
      for(l=0;l<prot.nNumDihedALL;++l){
	if (na==(atom_dihed_pair[l][0]-1)){
	  flag=0;
	  for (m=0;m<ll;++m){
	    if (atom_dihed_pair[l][3]==clust[k].o_f_pairs.o_f_not_interacting[kk][m]){
	      flag=1;
	    }
	  }
	  if (flag==0 && na+1 < atom_dihed_pair[l][3]){
	    clust[k].o_f_pairs.o_f_not_interacting[kk][ll] = atom_dihed_pair[l][3];
	    ++ll;
	  }
	}
 	else if (na==(atom_dihed_pair[l][3]-1)){
	  flag=0;
	  for (m=0;m<ll;++m){
	    if (atom_dihed_pair[l][0]==clust[k].o_f_pairs.o_f_not_interacting[kk][m]){
	      flag=1;
	    }
	  }
	  if (flag==0 && na+1 < atom_dihed_pair[l][3]){
	    clust[k].o_f_pairs.o_f_not_interacting[kk][ll] = atom_dihed_pair[l][0];
	    ++ll;
	  }
	}
      }
      ++na;
    }
    }*/

  na=0;
  for(k=0;k<prot.DOF;++k){
    for(kk=0;kk<clust[k].num_atom_clust;++kk){
      ll=0;
      for (i=0;i<AP.NUMEX[na];++i) {
	not_1_4[i]=1;
      }
      for(l=0;l<AP.NBONH;++l){
	if (na==(AP.BH[l][0])/3){
	  for (m=0;m<AP.NUMEX[na];++m) {
	    if ((AP.BH[l][1])/3+1 == clust[k].pairs.not_interacting[kk][m]) {
	      not_1_4[m]=0;
	    }
	    else {
	      //	      not_1_4[m]=1;
	    }
	  }
	}
      }
      for(l=0;l<AP.NBONA;++l){
	if (na==((AP.BA[l][0])/3)){
	  for (m=0;m<AP.NUMEX[na];++m) {
	    if ((AP.BA[l][1])/3+1 == clust[k].pairs.not_interacting[kk][m]) {
	      not_1_4[m]=0;
	    }
	    else {
	      //	      not_1_4[m]=1;
	    }
	  }
	}
      }
      for(l=0;l<AP.NTHETH;++l){
	if (na==((AP.TH[l][0])/3)){
	  for (m=0;m<AP.NUMEX[na];++m) {
	    if ((AP.TH[l][2])/3+1 == clust[k].pairs.not_interacting[kk][m]) {
	      not_1_4[m]=0;
	    }
	    else {
	      //	      not_1_4[m]=1;
	    }
	  }
	}
      }
      for(l=0;l<AP.NTHETA;++l){
	if (na==(AP.TA[l][0])/3){
	  for (m=0;m<AP.NUMEX[na];++m) {
	    if ((AP.TA[l][2])/3+1 == clust[k].pairs.not_interacting[kk][m]) {
	      not_1_4[m]=0;
	    }
	    else {
	      //	      not_1_4[m]=1;
	    }
	  }
	}
      }
      ll=0;
      for (m=0;m<AP.NUMEX[na];++m) {
	if (not_1_4[m]!=0) {
	  clust[k].o_f_pairs.o_f_not_interacting[kk][ll] = clust[k].pairs.not_interacting[kk][m];
	  ++ll;
	}
      }
      ++na;
    }
  }


  na=0;
  for (i=0;i<prot.DOF;++i) {
    for (j=0;j<clust[i].num_atom_clust;++j) {
      clust[i].mass_clust[j] = AP.AMASS[na];
      ++na;
    }
  }

  /************************************************************************************/
   if ((debug=fopen("debug.txt","w"))==NULL){
     printf("error!!\n");
     exit(1);
  }
  /* 										      */
  /* 										      */
  /* 										      */
  /* na=0;									      */
  /* for (i=0;i<prot.DOF;++i){							      */
  /*   for (j=0;j<clust[i].num_atom_clust;++j){					      */
  /*     ++na;									      */
  /*     fprintf(debug,"%d",na);						      */
  /*     for(k=0;clust[i].o_f_pairs.o_f_not_interacting[j][k]!=0;++k) {		      */
  /* 	fprintf(debug,"%5d ",clust[i].o_f_pairs.o_f_not_interacting[j][k]);	      */
  /*     }									      */
  /*     fprintf(debug,"\n");							      */
  /*   }									      */
  /* }										      */
  /* 										      */
  /* fclose(debug);								      */
  /************************************************************************************/
   
  	// 総二面角数の取得
	fprintf(debug,"%d",prot.nNumDihedALL);
	fprintf(debug,"\n");

	// 総二面角タイプ数の取得
	fprintf(debug,"%d",prot.nNumDihedType);
	fprintf(debug,"\n");

	// 二面角中の原子の番号の取得
	for (k=0;k<prot.nNumDihedALL;++k)
	{
	  for (i=0;i<6;++i)
	    {
	      fprintf(debug,"%5d ",atom_dihed_pair[k*6+i]);
	    }
	  fprintf(debug,"\n");
	}

	for (k=0;k<inpnumH;++k) {
	  fprintf(debug,"%5d ",inpindexH[k]);
	}
	fprintf(debug,"\n");
	for (k=0;k<inpnumA;++k) {
	  fprintf(debug,"%5d ",inpindexA[k]);
	}
	fprintf(debug,"\n");

	fclose(debug);

	for (i=0;i<inpnumH;++i)
	  inpindex[i]=inpindexH[i];
	for (i=0;i<inpnumA;++i)
	  inpindex[i+inpnumH]=inpindexA[i]+AP.NPHIH;
	inpnum=inpnumH+inpnumA;

  //	for(nNumDihedType=0;nNumDihedType<prot.nNumDihedType;++nNumDihedType)
  //{
  //	fprintf(debug,"%e ",V_dihed[nNumDihedType]/2.0);
  //	if ((nNumDihedType+1)%5 ==0){
  //	  fprintf(debug,"\n");		  
  //	}
  //}
  //fprintf(debug,"\n");

  //	for(nNumDihedType=0;nNumDihedType<prot.nNumDihedType;++nNumDihedType)
	  //	{
  //	fprintf(debug,"%d ",n_dihed[nNumDihedType]);
  //	if ((nNumDihedType+1)%5 ==0){
  //	  fprintf(debug,"\n");		  
  //	}
  //	//		fprintf(debug,"%lf",&x);
  //	//		n_dihed[nNumDihedType] = x;
  //}
  //fprintf(debug,"\n");

  //	for(nNumDihedType=0;nNumDihedType<prot.nNumDihedType;++nNumDihedType)
  //{
  //	fprintf(debug,"%lf ",theta_dihed[nNumDihedType]);
  //	if ((nNumDihedType+1)%5 ==0){
  //	  fprintf(debug,"\n");		  
  //	}
  //}
  //fprintf(debug,"\n");

// 静電相互作用項
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//	na=0;
  //for(k=0;k<prot.DOF;++k)
  //{
  //	for(kk=0;kk<clust[k].num_atom_clust;++kk)
  //	{
  //	  ++na;
  //		fprintf(debug,"%e ",clust[k].f_p_clust.e_f[kk]);
  //		if ((na+1)%5 ==0){
  //		  fprintf(debug,"\n");		  
  //		}
  //	}
  //}
  //fprintf(debug,"\n");

// VDW互作用項
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

  //fprintf(debug,"%d ",num_atomtype);
  //fprintf(debug,"\n");

  //for(n=0;n<prot.num_atom;++n)
  //{
  //  fprintf(debug,"%5d ",prot.L_J_parm.atomtype[n]);
  //  if ((n+1)%5 ==0){
  //    fprintf(debug,"\n");		  
  //  }
  //}
  //fprintf(debug,"\n");

  //na=0;
  //for(n=0;n<num_atomtype;++n){
  //  for(nn=0;nn<num_atomtype;++nn){
  //      ++na;
  //      fprintf(debug,"%5d ",prot.L_J_parm.atomtypeIndex[n][nn]);
  //      if ((na+1)%10 ==0){
  //	fprintf(debug,"\n");		  
  //      }
  //  }
  //}

  //	fprintf(debug,"%d \n",numtype);

  //for(k=0;k<numtype;++k)
  //
  //fprintf(debug,"%10.8e ",prot.L_J_parm.A[k]);
  //		if ((k+1)%5 ==0){
  //		  fprintf(debug,"\n");		  
  //		}
  //}
  //fprintf(debug,"\n");

  //for(k=0;k<numtype;++k)
  //{
  //  fprintf(debug,"%10.8e ",prot.L_J_parm.B[k]);
  //  if ((k+1)%5 ==0){
  //    fprintf(debug,"\n");		  
  //  }
  //}
  //fprintf(debug,"\n");

// 非結合性相互作用の組
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//	na=0;
//	for(k=0;k<prot.DOF;++k)
  //{
  //	for(kk=0;kk<clust[k].num_atom_clust;++kk)
  //	{
  //	  //		  fprintf(debug," %d ",na+1);
  //	  for(l=0;l<AP.NUMEX[na]/*+1*/;++l)
  //		{
  //			fprintf(debug," %d ",clust[k].pairs.not_interacting[kk][l]);
  //
  //		}
  //		fprintf(debug,"\n ");
  //	  ++na;
  //	}
  //}
// 1-4非結合性相互作用の組
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

  //for(k=0;k<prot.DOF;++k)
  //{
  //	for(kk=0;kk<clust[k].num_atom_clust;++kk)
  //	{
  //		for(l=0;clust[k].o_f_pairs.o_f_not_interacting[kk][l]!=0;++l)
  //		{
  //			fprintf(debug,"%5d ",clust[k].o_f_pairs.o_f_not_interacting[kk][l]);
  //		}
  //		fprintf(debug,"0");
  //		fprintf(debug,"\n");
  //	}
  //}

  DOFOFPROT=3*AP.NATOM-AP.NBONH-AP.NBONA-AP.NTHETH-AP.NTHETA-6;
  if (DOFOFPROT <= 0)
    DOFOFPROT+=6;


  //fclose(debug);
}

// 剛体中の全ての 2 面角の計算を行う関数
double pick_dihed_one_clust3(int nnumaom1,int nnumaom2,int nnumaom3,int nnumaom4,int nNumDihed) {
  int alpha;
  
  double atom_i[3],atom_j[3],atom_k[3],atom_l[3];
  
  double vec_ij[3],vec_jk[3],vec_kl[3];
  double out_ij_jk[3],out_jk_kl[3],out_ijkl_jkkl[3];

  double d_ij_jk=0.0,d_jk_kl=0.0,d_ijkl_jkkl=0.0,d_jk=0.0;

  double det=0.0;

  double cs=0.0;
  double theta;
  double pi;

  double out_ijjk_jk[3];
  double in_ijjkjk_jkkl=0.0,in_ijjk_jkkl=0.0;

  // 座標の取得
  for(alpha=0;alpha<3;++alpha){
    atom_i[alpha] = coord_rest[nnumaom1][alpha];
    atom_j[alpha] = coord_rest[nnumaom2][alpha];
    atom_k[alpha] = coord_rest[nnumaom3][alpha];
    atom_l[alpha] = coord_rest[nnumaom4][alpha];
  }

  for (alpha=0;alpha<3;++alpha) {
    vec_ij[alpha] = atom_j[alpha]-atom_i[alpha];
    vec_jk[alpha] = atom_k[alpha]-atom_j[alpha];
    vec_kl[alpha] = atom_l[alpha]-atom_k[alpha];
  }

  out_ij_jk[0]=vec_ij[1]*vec_jk[2]-vec_ij[2]*vec_jk[1];
  out_ij_jk[1]=vec_ij[2]*vec_jk[0]-vec_ij[0]*vec_jk[2];
  out_ij_jk[2]=vec_ij[0]*vec_jk[1]-vec_ij[1]*vec_jk[0];
			
  out_jk_kl[0]=vec_jk[1]*vec_kl[2]-vec_jk[2]*vec_kl[1];
  out_jk_kl[1]=vec_jk[2]*vec_kl[0]-vec_jk[0]*vec_kl[2];
  out_jk_kl[2]=vec_jk[0]*vec_kl[1]-vec_jk[1]*vec_kl[0];

  d_ij_jk += out_ij_jk[0]*out_ij_jk[0]+out_ij_jk[1]*out_ij_jk[1]+out_ij_jk[2]*out_ij_jk[2];
  d_jk_kl += out_jk_kl[0]*out_jk_kl[0]+out_jk_kl[1]*out_jk_kl[1]+out_jk_kl[2]*out_jk_kl[2];

  d_ij_jk = sqrt(d_ij_jk);
  d_jk_kl = sqrt(d_jk_kl);

  for (alpha=0;alpha<3;++alpha) {
    out_ij_jk[alpha] = out_ij_jk[alpha]/d_ij_jk;
    out_jk_kl[alpha] = out_jk_kl[alpha]/d_jk_kl;
  }

  for(alpha=0;alpha<3;++alpha) {
    cs += out_ij_jk[alpha]*out_jk_kl[alpha];
  }

  if (cs < -1.0 ){
    cs = -1.0;
  }
  else if (cs > 1.0 ) {
    cs = 1.0;
  }

  out_ijkl_jkkl[0] = out_ij_jk[1]*out_jk_kl[2]-out_ij_jk[2]*out_jk_kl[1];
  out_ijkl_jkkl[1] = out_ij_jk[2]*out_jk_kl[0]-out_ij_jk[0]*out_jk_kl[2];
  out_ijkl_jkkl[2] = out_ij_jk[0]*out_jk_kl[1]-out_ij_jk[1]*out_jk_kl[0];

  det = out_ijkl_jkkl[0]*vec_jk[0]+out_ijkl_jkkl[1]*vec_jk[1]+out_ijkl_jkkl[2]*vec_jk[2];
  
  d_ijkl_jkkl += out_ijkl_jkkl[0]*out_ijkl_jkkl[0]+out_ijkl_jkkl[1]*out_ijkl_jkkl[1]+out_ijkl_jkkl[2]*out_ijkl_jkkl[2];
  d_ijkl_jkkl = sqrt(d_ijkl_jkkl);
  d_jk += vec_jk[0]*vec_jk[0]+vec_jk[1]*vec_jk[1]+vec_jk[2]*vec_jk[2];
  d_jk = sqrt(d_jk);

  det = det/(d_ijkl_jkkl*d_jk);
  theta = acos(cs)*det;
	
  pi = acos(-1.0);
  if (det<0) {
    theta = 2.0*pi+theta;
  }

  ///////////////////////////////////////////////////////////////
  /*******************************************************************/
  /* out_ijjk_jk[0] = out_ij_jk[1]*vec_jk[2]-out_ij_jk[2]*vec_jk[1]; */
  /* out_ijjk_jk[1] = out_ij_jk[2]*vec_jk[0]-out_ij_jk[0]*vec_jk[2]; */
  /* out_ijjk_jk[2] = out_ij_jk[0]*vec_jk[1]-out_ij_jk[1]*vec_jk[0]; */
  /* 								     */
  /* for (alpha=0;alpha<3;++alpha) {				     */
  /*   in_ijjkjk_jkkl += out_ijjk_jk[alpha]*out_jk_kl[alpha];	     */
  /* }								     */
  /* 								     */
  /* for (alpha=0;alpha<3;++alpha) {				     */
  /*   in_ijjk_jkkl += out_ij_jk[alpha]*out_jk_kl[alpha];	     */
  /* }								     */
  /* 								     */
  /* theta = atan2(in_ijjkjk_jkkl,in_ijjk_jkkl);		     */
  /*******************************************************************/
	
  return theta/*rad*/;

}

void optimize_q_value(double omega,double *q_value) {
  int dof;
  FILE *debug;

  if (TermMoveMode2==12)
    dof=prot.DOF-1+6;
  else
    dof=prot.DOF-1;
  
  *q_value=dof*k_B_kcm*T_Kelvin/4.18407/100.0/(omega)/(omega);

  debug=fopen("debug_op.txt","w");
  fprintf(debug,"q=%e\nomega=%e\n",*q_value,omega);
  fclose(debug);
}
