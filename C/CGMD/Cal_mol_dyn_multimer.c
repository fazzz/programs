#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gener.h"
#include "ABA.h" // #include "ABA_multimer.h"
#include "physics.h"
#include "MD.h"
#include "force.h"
#include "BD.h"


void pre_dyn(void);
void scaling_momentum(int nNumStep);

// ¥¤¥ó¥×¥Ã¥È¥Õ¥¡¥¤¥ë¤ÎÆÉ¤ß¹þ¤ß¡¢¥·¥ß¥å¥ì¡¼¥·¥ç¥ó¤Î¼Â¹Ô¡¢
// ·×»»·ë²Ì¤Î½ÐÎÏ¤ò¹Ô¤¦¥á¥¤¥ó´Ø¿ô
int main(int argc, char *argv[]) {
  int i,j,k; 
  int alpha;
  
  int time,i_c;
  
  int opt;
  
  int nNumClut, nNumClut_dammy;
  
  int nNumClutOrigBranch=0;
  int nNumClutOrigBranch2=0;
  
  int pflag;
  int cflag;
  int numDOF;
  
  FILE *output_c;
  FILE *output_v;
  FILE *output2_c[MAXPDB];
  FILE *outputMM;
  FILE *outputthomo;
  FILE *outputMap;
  FILE *coord_in_restat;
  FILE *velo_in_restat;
  FILE *outtest_dihed;
  FILE *rst_c;
  FILE *rst_v;
  FILE *rst_a;
  
  FILE *db;
  FILE *d;

  char filename[MAXPDB];
  
  char *option;
  
  double T_Kelvin_Now;
  double inertiaofvertialvariable,inertiaofvertialvariable2;
	//	double vel_db[100][3];
  double vel_Term_bh[6],vel_Term_bh_dummy[6],vel_Term_boh[6],acc_Term_bh[6];
  double predict_Term[6][6],correct_Term[6][6],delta_acc_Term[6];
  double predict_Term3[6][6],correct_Term3[6][6],delta_acc_Term3[6];
  double q_Term[4],predict_q_Term[4][5],correct_q_Term[4][5],delta_q_Term[4],len;
  
  int joinflag;
  
  double *ele, *ALJ, *BLJ;
  double *p_e, *p_1_4_e; 
  double *p_LJ,*p_1_4_LJ;
  double *f_e, *f_1_4_e; 
  double *f_LJ, *f_1_4_LJ;
  int numnb, num14;
  int *indexnb, *index14 ;
  double *vel,*vel_bh,*vel_boh,*vel_bh_dummy,dot_s_NVT_bh,dot_s_NVT_bh_dummy,dot_s_NVT_boh,*acc_bh;
  //  double inertiaofvertialvariable2;
  
  double *cord;
  double *coord_rst;
  int tnum_atom_clust;
  int tDOF;

  double*S0,*omega,*iniV;

  double *eig,*eig_14,*eig_dihed;
  FILE *eiginputfile;

  InpfilCOORD = "crd.in";
  InpfilCOORDREF = "crd.in";
  InpfilCLUST = "clust.in";
  InpfilSEQ = "seq.in";
  InpfilTOP = "top.in";
  InpfilMDI = "md.in";
  InpfilDVELO = "velo.in";
  InpfilEig= "eig1.txt";
  OutfilTHMO = "thermo_dyn_properties.out";
  OutfilCOORD = "coo_pro.out";
  OutfilVELO = "velo_pro.out";
  OutfilRESTAC = "crd_rst.in";
  OutfilRESTAV = "velo_new.in";
  OutfilRESTAMFORM = "prot.rst";
//	MM = 6;

  // ¥ª¥×¥·¥ç¥ó¤òÆÉ¤à
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  // ¥ª¥×¥·¥ç¥ó¤¬Â¸ºß
  if (argc > 1){
    i = 1;
    // ¥ª¥×¥·¥ç¥ó¤ò¤¹¤Ù¤ÆÆÉ¤ó¤À¤«
    while (++i < argc){
      // ¥ª¥×¥·¥ç¥ó¤«
      if ( (*++argv)[0] == '-'){
	// "-x"¤Î"x"¤òÆÉ¤à
	option = *argv;
	opt = *++option;
	if (opt == 'c'){
	  ++i;
	  // ºÂÉ¸¥Õ¥¡¥¤¥ëÌ¾¤òÆÉ¤à
	  if (*++argv != NULL){
	    InpfilCOORD = *argv;
	  }
	}
	else if (opt == 'r'){
	  ++i;
	  // ¹äÂÎ¥Õ¥¡¥¤¥ëÌ¾¤òÆÉ¤à
	  if (*++argv != NULL){
	    InpfilCLUST = *argv;
	  }
	}
	else if (opt == 'p'){
	  ++i;
	  // ¥È¥Ý¥í¥¸¡¼¥Õ¥¡¥¤¥ëÌ¾¤òÆÉ¤à
	  if (*++argv != NULL){
	    InpfilTOP = *argv;
	    pflag = AMBERMODE;
	  }
	}
	else if (opt == 't'){
	  ++i;
	  // ¥È¥Ý¥í¥¸¡¼¥Õ¥¡¥¤¥ëÌ¾¤òÆÉ¤à
	  if (*++argv != NULL){
	    InpfilTOP = *argv;
	  }
	}
	else if (opt == 's'){
	  ++i;
	  // »Ä´ð¥Õ¥¡¥¤¥ëÌ¾¤òÆÉ¤à
	  if (*++argv != NULL){
	    InpfilSEQ  = *argv;
	  }
	}
	else if (opt == 'm'){
	  ++i;
	  // »Ä´ð¥Õ¥¡¥¤¥ëÌ¾¤òÆÉ¤à
	  if (*++argv != NULL){
	    InpfilMDI = *argv;
	  }
	}
	else if (opt == 'f'){
	  ++i;
	  // »²¾È¹½Â¤¥Õ¥¡¥¤¥ëÌ¾¤òÆÉ¤à
	  if (*++argv != NULL){
	    InpfilCOORDREF = *argv;
	    restflag = ON;
	  }
	}
	else if (opt == 'd'){
	  ++i;
	  // Â®ÅÙ¥Õ¥¡¥¤¥ëÌ¾¤òÆÉ¤à
	  if (*++argv != NULL){
	    InpfilDVELO = *argv;
	  }
	}				  
	else if (opt == 'e'){
	  ++i;
	  // Â®ÅÙ¥Õ¥¡¥¤¥ëÌ¾¤òÆÉ¤à
	  if (*++argv != NULL){
	    InpfilEig = *argv;
	  }
	}				  
	else if (opt == 'o'){
	  ++i;
	  // ¥¢¥¦¥È¥×¥Ã¥È¥Õ¥¡¥¤¥ëÌ¾¤òÆÉ¤à
	  if (*++argv != NULL){
	    OutfilTHMO = *argv;
	  }
	}				  
	else if (opt == 'x'){
	  ++i;
	  // ¥¢¥¦¥È¥×¥Ã¥È¥Õ¥¡¥¤¥ëÌ¾¤òÆÉ¤à
	  if (*++argv != NULL) {
	    OutfilCOORD = *argv;
	  }
	}				  
	else if (opt == 'v'){
	  ++i;
	  // ¥¢¥¦¥È¥×¥Ã¥È¥Õ¥¡¥¤¥ëÌ¾¤òÆÉ¤à
	  if (*++argv != NULL){
	    OutfilVELO = *argv;
	  }
	}
	else if (opt == 'j'){
	  ++i;
	  // ¥¢¥¦¥È¥×¥Ã¥È¥Õ¥¡¥¤¥ëÌ¾¤òÆÉ¤à
	  if (*++argv != NULL){
	    OutfilRESTAC = *argv;
	  }
	}				  
	else if (opt == 'l'){
	  ++i;
	  // ¥¢¥¦¥È¥×¥Ã¥È¥Õ¥¡¥¤¥ëÌ¾¤òÆÉ¤à
	  if (*++argv != NULL){
	    OutfilRESTAV = *argv;
	  }
	}				  
	else if (opt == 'k'){
	  ++i;
	  // ¥¢¥¦¥È¥×¥Ã¥È¥Õ¥¡¥¤¥ëÌ¾¤òÆÉ¤à
	  if (*++argv != NULL){
	    OutfilRESTAMFORM = *argv;
	  }
	}				  
	else{
	  // ´Ö°ã¤¤¤ò»ØÅ¦
	  printf("illegal option !!! \n");
	  exit(1);
	}
      }
    }
  }
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

  // ºÂÉ¸¥Õ¥¡¥¤¥ë¡¢¥È¥Ý¥í¥¸¡¼¥Õ¥¡¥¤¥ë¡¢·×»»¾ò·ï¥Õ¥¡¥¤¥ë¤ÎÆÉ¤ß¹þ¤ß
  //  printf("yes\n");
  pick_data(pflag);
  //  printf("yes\n");
  joinflag=0;
  for(nNumClut=0; nNumClut<prot.DOF; ++nNumClut)	{
    joinflag=setJoin(nNumClut,joinflag);
  }
  
  // ºÂÉ¸¥Õ¥¡¥¤¥ë¤Î¥Ç¡¼¥¿¤ò¼Â¸³¼¼·Ï¤«¤é¥¯¥é¥¹¥¿·Ï¤ØÊÑ´¹
  for(nNumClut=0; nNumClut<prot.DOF; ++nNumClut){
    trans_A_to_CN(nNumClut);
  }
  q_Term[0]=1.0;   	
  for (i = 1; i < 3; ++i) {
    q_Term[i]=0.0;   	
  }

  // ´·À­¥â¡¼¥á¥ó¥È¤Ê¤É¤Î·×»»
  set();
  
  // ¥¹¥¿¡¼¥È¤«¡¢ºÆ¥¹¥¿¡¼¥È¤«
  if (MODEV == ON){
    // ½é´üÂ®ÅÙ¤Î·×»»
    set_initial_velo();
  }
  else{
    // ½é´üÂ®ÅÙ¤ÎÆÉ¤ß¹þ¤ß
    pick_initial_velo();
  }
  
  // ½ÐÎÏ¥Õ¥¡¥¤¥ë¤Î¥ª¡¼¥×¥ó_1
  if ((output_c = fopen(/*"coo_pro.out"*/OutfilCOORD,"a")) == NULL){
    printf("error coo_pro.out cannot open\n");
    exit(1);
  }
  // ½é´ü¹½Â¤¤Î out ¥Õ¥¡¥¤¥ë½ÐÎÏ
  output_file_coo_pro(output_c);
  fclose(output_c);
  
  // ±¿Æ°¥¨¥Í¥ë¥®¡¼¡¢¥Ý¥Æ¥ó¥·¥ã¥ë¥¨¥Í¥ë¥®¡¼¡¢Á´¥¨¥Í¥ë¥®¡¼¤Î½é´ü²½
  Energy_potential_stac = 0.0;
  Energy_kinetic_stac   = 0.0;
  Energy_Internal_stac  = 0.0;
  //    printf("yes2\n");
  set_delts_matrix();
  
  for (nNumClut=1;nNumClut<prot.DOF;++nNumClut){
    for (i=0;i<6;++i){
      clust[nNumClut].correct_dihedang[i] = 0.0;
    }
  }
  
  // ÊÑ´¹¹ÔÎó¤ÎÀßÄê
  for(nNumClut=0; nNumClut<prot.DOF; ++nNumClut){
    if (clust[nNumClut].num_branch > 1){
      nNumClutOrigBranch2 = nNumClut;
    }
    set_trans_Matrix(nNumClut, nNumClutOrigBranch2);
  }

  /**************************/
  for (nNumClut=1;nNumClut<prot.DOF;++nNumClut)	{
    clust[nNumClut].correct_dihedang[0] = clust[nNumClut].dihedang[0]/**deltat*/;
    clust[nNumClut].correct_dihedang[1] = clust[nNumClut].ddihedang[0]*deltat;
  }
  //	GearMODE = OFF;
  if (GearMODE == OFF)	{
    //////////////////////////////////////////////////////////////////////////
    for (nNumClut=1;nNumClut<prot.DOF;++nNumClut){
      /*clust[nNumClut].*/old_dddihedang[nNumClut] = 0.0;
    }
    if (MODEV == ON){
      // ½é´üÂ®ÅÙ¤Î·×»»
      set_initial_velo();
    }
    for (nNumClut=1;nNumClut<prot.DOF;++nNumClut) {
      //			old_deltadihedang[nNumClut] = deltat*clust[nNumClut].ddihedang[0];
    }
    //		for(nNumClut=1; nNumClut<prot.DOF; ++nNumClut)
    //		{
    //			if (clust[nNumClut].num_branch > 1)
    //			{
    //				nNumClutOrigBranch = nNumClut;
    //			}
    //			trans_CN_to_A(nNumClut, nNumClutOrigBranch);
    //		}
  }
  //////////////////////////////////////////////////////////////////////////
  
  GearsConstant[0] = 3.0/16.0;
  GearsConstant[1] = 251.0/360.0;
  GearsConstant[2] = 1.0;
  GearsConstant[3] = 11.0/18.0;
  GearsConstant[4] = 1.0/6.0;
  GearsConstant[5] = 1.0/60.0;

  GearsConstant2[0] = 3.0/20.0;
  GearsConstant2[1] = 251.0/360.0;
  GearsConstant2[2] = 1.0;
  GearsConstant2[3] = 11.0/18.0;
  GearsConstant2[4] = 1.0/6.0;
  GearsConstant2[5] = 1.0/60.0;

  GearsConstant3[0] = 251.0/720.0;
  GearsConstant3[1] = 1.0;
  GearsConstant3[2] = 11.0/12.0;
  GearsConstant3[3] = 1.0/3.0;
  GearsConstant3[4] = 1.0/24.0;
  
  for (i=0;i<6;++i){
    for (j=0;j<6;++j) {
      if (i != j){
	Telar_Matrix[i][j] = 0.0;
      }
      else{
	Telar_Matrix[i][i] = 1.0;
      }
    }
  }
  
  for (i=0;i<5;++i){
    for (j=0;j<5;++j){
      if (i != j){
	Telar_Matrix2[i][j] = 0.0;
      }
      else{
	Telar_Matrix2[i][i] = 1.0;
      }
    }
  }
	
  Telar_Matrix[0][1] = 1.0;
  Telar_Matrix[0][2] = 1.0;
  Telar_Matrix[0][3] = 1.0;
  Telar_Matrix[0][4] = 1.0;
  Telar_Matrix[0][5] = 1.0;
  Telar_Matrix[1][2] = 2.0;
  Telar_Matrix[1][3] = 3.0;
  Telar_Matrix[1][4] = 4.0;
  Telar_Matrix[1][5] = 5.0;
  Telar_Matrix[2][3] = 3.0;
  Telar_Matrix[2][4] = 6.0;
  Telar_Matrix[2][5] = 10.0;
  Telar_Matrix[3][4] = 4.0;
  Telar_Matrix[3][5] = 10.0;
  Telar_Matrix[4][5] = 5.0;
  
  Telar_Matrix2[0][1] = 1.0;
  Telar_Matrix2[0][2] = 1.0;
  Telar_Matrix2[0][3] = 1.0;
  Telar_Matrix2[0][4] = 1.0;
  Telar_Matrix2[1][2] = 2.0;
  Telar_Matrix2[1][3] = 3.0;
  Telar_Matrix2[1][4] = 4.0;
  Telar_Matrix2[2][3] = 3.0;
  Telar_Matrix2[2][4] = 6.0;
  Telar_Matrix2[3][4] = 4.0;
	
  GearsConstant[0] = 3.0/16.0;
  GearsConstant[1] = 251.0/360.0;
  GearsConstant[2] = 1.0;
  GearsConstant[3] = 11.0/18.0;
  GearsConstant[4] = 1.0/6.0;
  GearsConstant[5] = 1.0/60.0;
  
  GearsConstant2[0] = 3.0/20.0;
  GearsConstant2[1] = 251.0/360.0;
  GearsConstant2[2] = 1.0;
  GearsConstant2[3] = 11.0/18.0;
  GearsConstant2[4] = 1.0/6.0;
  GearsConstant2[5] = 1.0/60.0;

  GearsConstant3[0] = 251.0/720.0;
  GearsConstant3[1] = 1.0;
  GearsConstant3[2] = 11.0/12.0;
  GearsConstant3[3] = 1.0/3.0;
  GearsConstant3[4] = 1.0/24.0;

//	dot_dot_s_NVT_old = 0.0;
//	dot_s_NVT = 0.0;


  nNumStep=0;
  if ((ele=malloc(sizeof(double)*prot.num_atom)) == NULL) {
    printf("error : alocate of ele\n");
  }
  if ((ALJ=malloc(sizeof(double)*prot.num_atom*prot.num_atom)) == NULL) {
    printf("error : alocate of ALJ\n");
  }
  if ((BLJ=malloc(sizeof(double)*prot.num_atom*prot.num_atom)) == NULL) {
    printf("error : alocate of BLJ\n");
  }
  set_non_bonding_parameters(ele,ALJ,BLJ);
  if ((p_e=(double *)calloc(prot.num_atom,sizeof(double))) == NULL) {
    printf("error : alocate of p_e\n");
  }
  if ((p_1_4_e=(double *)calloc(prot.num_atom,sizeof(double))) == NULL) {
    printf("error : alocate of p_1_4_e\n");
  }
  if ((p_LJ=(double *)calloc(prot.num_atom,sizeof(double))) == NULL) {
    printf("error : alocate of p_LJ\n");
  }
  if ((p_1_4_LJ=(double *)calloc(prot.num_atom,sizeof(double))) == NULL) {
    printf("error : alocate of p_1_4_LJ\n");
  }
  if ((f_e=(double *)calloc(prot.num_atom*3,sizeof(double))) == NULL) {
    printf("error : alocate of f_e\n");
  }
  if ((f_1_4_e=(double *)calloc(prot.num_atom*3,sizeof(double))) == NULL) {
	  printf("error : alocate of f_1_4_e\n");
  }
  if ((f_LJ=(double *)calloc(prot.num_atom*3,sizeof(double))) == NULL) {
    printf("error : alocate of f_1_4_e\n");
  }
  if ((f_1_4_LJ=(double *)calloc(prot.num_atom*3,sizeof(double))) == NULL) {
    printf("error : alocate of f_1_4_LJ\n");
  }

  if (GearMODE == OFF) {
    if ((vel_bh=(double *)calloc(prot.DOF,sizeof(double))) == NULL) {
      printf("error : alocate of vel_bh\n");
    }
    if ((acc_bh=(double *)calloc(prot.DOF,sizeof(double))) == NULL) {
      printf("error : alocate of acc_bh\n");
    }
    if ((vel_boh=(double *)calloc(prot.DOF,sizeof(double))) == NULL) {
      printf("error : alocate of vel_boh\n");
    }
    if ((vel_bh_dummy=(double *)calloc(prot.DOF,sizeof(double))) == NULL) {
      printf("error : alocate of vel_boh\n");
    }
    for (i=0;i<prot.DOF;++i){
      vel_bh[i] = clust[i].ddihedang[0];
      vel_boh[i]= clust[i].ddihedang[0];
      if (TermMoveMode==ON) {
	for (j=0;j<6;++j) {
	  acc_Term[j]=0.0;
	  acc_Term_bh[j]=0.0;
	  vel_Term[j]=0.0;
	  vel_Term_bh[j]=0.0;
	  vel_Term_boh[j]=0.0;
	  delta_Term[j]=0.0;
	}
      }
    }
    if (MODE == NVT) {
      dot_s_NVT_bh=dot_s_NVT;
      dot_s_NVT_boh=dot_s_NVT;
    }
  }
  
  set_non_bonding_index();


  numnb=gnumnb;
  num14=gnum14;
  indexnb=(int *)calloc(numnb*2,sizeof(int));
  index14=(int *)calloc(num14*2,sizeof(int));
  for (i=0;i<numnb;++i){
    indexnb[i*2]=gindexnb[i*2];
    indexnb[i*2+1]=gindexnb[i*2+1];
  }
  for (i=0;i<num14;++i){
    index14[i*2]=gindex14[i*2];
    index14[i*2+1]=gindex14[i*2+1];
  }

  printf("%d %d %d\n",numnb,num14,prot.nNumDihedALL);
  printf("%d \n",PEPCAACCMODE);
  
  if (PEPCAACCMODE >= 1/*== ON*/) {
    if ((eig=(double *)malloc(sizeof(double)*numnb*2))==NULL) {
      exit(1);
      printf("error of allocation\n");
    }
    if ((eig_14=(double *)malloc(sizeof(double)*num14*2))==NULL) {
      exit(1);
      printf("error of allocation\n");
    }
    if ((eig_dihed=(double *)malloc(sizeof(double)*prot.nNumDihedALL))==NULL) {
      exit(1);
      printf("error of allocation\n");
    }
    if ((eiginputfile=fopen(InpfilEig,"r"))==NULL) {
      exit(1);
      printf("error of file open\n");
    }
    for (i=0;i<numnb;++i) {
      fscanf(eiginputfile,"%lf",&eig[i*2]);
    }
    for (i=0;i<numnb;++i) {
      fscanf(eiginputfile,"%lf",&eig[i*2+1]);
    }
    for (i=0;i<num14;++i) {
      fscanf(eiginputfile,"%lf",&eig_14[i*2]);
    }
    for (i=0;i<num14;++i) {
      fscanf(eiginputfile,"%lf",&eig_14[i*2+1]);
    }
    for (i=0;i<prot.nNumDihedALL;++i) {
      fscanf(eiginputfile,"%lf",&eig_dihed[i]);
    }
    fscanf(eiginputfile,"%lf",&fact);
    fclose(eiginputfile);
    printf("%lf\n",fact);
    for (i=0;i<numnb;++i) {
      printf("eig_es_%d= %lf\n",i,eig[i*2]);
    }
    for (i=0;i<numnb;++i) {
      printf("eig_LJ_%d=%lf\n",i,eig[i*2+1]);
    }
    for (i=0;i<num14;++i) {
      printf("eig_es_1-4_%d=%lf\n",i,eig_14[i*2]);
    }
    for (i=0;i<num14;++i) {
      printf("eig_LJ_1-4_%d=%lf\n",i,eig_14[i*2+1]);
    }
    for (i=0;i<prot.nNumDihedALL;++i) {
      printf("eig_di_%d=%lf\n",i,eig_dihed[i]);
    }
    printf("\n");
  }
  

  if (dhstopflag==0/*ON && dhstopflag!=2*/) {
    // 2 ÌÌ³ÑÁê¸ßºîÍÑ¤Î·×»»
    Calc_dihed_Potential_Force(eig_dihed);
    //  printf("yes-0\n");
  }
  else if (dhstopflag==2) {
    // 2 ÌÌ³ÑÁê¸ßºîÍÑ¤Î·×»»
    Calc_dihed_Potential_Force_for_db();
  }
  //  printf("yes1\n");  
  if(nbstopflag!=ON) {
    // VDW Áê¸ßºîÍÑ¤Î·×»»
    if (pflag!=AMBERMODE) {
      Calc_L_J_PotentialandForce();
    }
    else {
      cord=(double *)malloc(sizeof(double)*prot.num_atom*3);
      for (i=0;i<prot.num_atom;++i) {
	for (alpha=0;alpha<3;++alpha) {
	  cord[i*3+alpha]=prot.coord[i][alpha];
	}
      }
      
      Calc_L_J_PotentialandForce2(ele, ALJ, BLJ,
				  p_e,  p_1_4_e, 
				  p_LJ, p_1_4_LJ,
				  f_e, f_1_4_e, 
				  f_LJ,f_1_4_LJ,
				  numnb, indexnb,
				  num14, index14, cord,eig,eig_14);
      //  printf("yes1-2\n");
	    /******************************************************************************************/
            /* num_a_prot=0;									      */
	    /* for(nNumClut = 0;nNumClut < tDOF; ++nNumClut) {					      */
	    /*   tnum_atom_clust = clust[nNumClut].num_atom_clust;				      */
	    /*   for(i=0;i < tnum_atom_clust;++i) {						      */
	    /* 	potential_pro.p_elesta[num_a_prot] = p_e[num_a_prot];				      */
	    /* 	potential_pro.p_1_4_elesta[num_a_prot] = p_1_4_e[num_a_prot];			      */
	    /* 	potential_pro.p_L_J[num_a_prot] = p_LJ[num_a_prot];				      */
	    /* 	potential_pro.p_1_4_L_J[num_a_prot] = p_1_4_LJ[num_a_prot];			      */
	    /* 	for(alpha=0; alpha<3; ++alpha) {						      */
	    /* 	  clust[nNumClut].f_c.f_L_J[i][alpha] = f_LJ[num_a_prot*3+alpha];		      */
	    /* 	  clust[nNumClut].f_c.f_1_4_L_J[i][alpha] = f_1_4_LJ[num_a_prot*3+alpha];	      */
	    /* 	  clust[nNumClut].f_c.f_elesta[i][alpha] = f_e[num_a_prot*3+alpha];		      */
	    /* 	  clust[nNumClut].f_c.f_1_4_elesta[i][alpha] = f_1_4_e[num_a_prot*3+alpha];	      */
	    /* 	}										      */
	    /* 	++num_a_prot;									      */
	    /*   }										      */
	    /* }										      */
            /******************************************************************************************/


    }
    if (restflag==ON) {
      // À©¸Â¤ò¤«¤±¤ë
      //Calc_restraintForce();
    }
  }
  
  T_Kelvin_Now = analthermo_dyn_properties(outputthomo,vel_Term);
  //  printf("yes1-3\n");
	if (TermMoveMode2 > 6) {
          /**********************************************************************************/
          /* for (i=0;i<3;++i)								    */
	  /*   vel_Term[i]=1.0;								    */
	  /* vel_Term[3]=-clust[0].qCOM[2]*vel_Term[1]+clust[0].qCOM[1]*vel_Term[2];	    */
	  /* vel_Term[4]= clust[0].qCOM[2]*vel_Term[0]-clust[0].qCOM[0]*vel_Term[2];	    */
	  /* vel_Term[5]=-clust[0].qCOM[1]*vel_Term[0]+clust[0].qCOM[0]*vel_Term[1];	    */
          /**********************************************************************************/
	}


	if (TermMoveMode==ON) {
	  q_Term[0]=1.0;
	  q_Term[1]=0.0;
	  q_Term[2]=0.0;
	  q_Term[3]=0.0;
	  vel_q_Term[i]=0.0;
	  for (i=0;i<3;++i)
	    for (j=0;j<3;++j)
	      Rot_Term[i][j]=0.0;
	  for (i=0;i<3;++i)
	      Rot_Term[i][i]=1.0;
	  for (i=0;i<4;++i)
	    for (j=0;j<5;++j)
	      correct_q_Term[i][0]=0.0;
	  for (i=0;i<4;++i)
	    correct_q_Term[i][0]=q_Term[i];
	  for (i=0;i<6;++i)
	    correct_Term[i][1]=vel_Term[i]*deltat;
	}
	if (TermMoveMode2 == 2 /*|| TermMoveMode2 == 4*/) {
	  S0=malloc(sizeof(double)*3);
	  omega=malloc(sizeof(double)*3*3);
	  iniV=malloc(sizeof(double)*3);
	  for (i=0;i<3;++i)
	    correct_Term[i][1]=0.01;
	  omega[0]= 0.0;  
	  omega[1]=-correct_Term[2][1]/deltat;  
	  omega[2]= correct_Term[1][1]/deltat;
	  omega[3]= correct_Term[2][1]/deltat;  
	  omega[4]= 0.0;  
	  omega[5]=-correct_Term[0][1]/deltat;
	  omega[6]=-correct_Term[1][1]/deltat;  
	  omega[7]= correct_Term[0][1]/deltat;  
	  omega[8]= 0.0;
	  for (i=0;i<3;++i)
	    S0[i]=clust[0].qCOM[i];
	  mvmult(omega,S0,iniV,3);
	  correct_Term[3][1]=-iniV[0]*deltat;
	  correct_Term[4][1]=iniV[1]*deltat;
	  correct_Term[5][1]=-iniV[2]*deltat;
	  for(i=0;i<3;++i)
	    Euler[i]=0.0;

	}
	/* for (i=3;i<6;++i)		     */
	//  printf("yes1-4\n");
	/*   correct_Term[i][1]=1.00;	     */
        /*************************************/
	// ¥·¥ß¥å¥ì¡¼¥·¥ç¥ó¤Î³«»Ï ////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	for(nNumStep = 1; nNumStep <= TIME_LIMIT; ++nNumStep) {
	  T_Kelvin_Now = analthermo_dyn_properties(outputthomo,vel_Term);
	  /***********************************************************/
          /* if ((d = fopen("db.txt","a")) == NULL){		     */
	  /*   printf("error db.txt cannot open\n");		     */
	  /*   exit(1);						     */
	  /* }							     */
	  /* fprintf(d,"%d ",nNumStep);				     */
	  /* for(nNumClut=0;nNumClut<prot.DOF;++nNumClut){	     */
	  /*   fprintf(d,"%e ",clust[nNumClut].ddihedang[0]);	     */
	  /* }							     */
	  /* fprintf(d,"\n");					     */
          /***********************************************************/
	  //  printf("yes2\n");
	  
	  if (DYNMODE == MD || DYNMODE == LD){
	    if (GearMODE == OFF){
	      numDOF=prot.DOF;
	      if (TermMoveMode == ON) {
		for (i=0;i<6;++i) {
		  vel_Term[i]=1.5*vel_Term_bh[i]-0.5*vel_Term_boh[i]; // 0811
		  //		  vel_Term[i]+=deltat*acc_Term[i];                    // 0811
		}
	      }
	      for (nNumClut=1;nNumClut<numDOF;++nNumClut) {
                clust[nNumClut].ddihedang[0] = 1.5*vel_bh[nNumClut]-0.5*vel_boh[nNumClut]; // 0811
		//                clust[nNumClut].ddihedang[0] = vel_bh[nNumClut]+deltat*clust[nNumClut].dddihedang[0]; // 0811
		//                clust[nNumClut].ddihedang[0] += deltat*clust[nNumClut].dddihedang[0];
		//		vel_bh_dummy[nNumClut] = clust[nNumClut].ddihedang[0];
	      }
	      if (MODE==NVT) {
		dot_s_NVT = 1.5*dot_s_NVT_bh-0.5*dot_s_NVT_boh;
	      }
              /***********************************************************/
              /* if ((d = fopen("db.txt","a")) == NULL){		 */
	      /* 	printf("error db.txt cannot open\n");		 */
	      /* 	exit(1);					 */
	      /* }							 */
	      /* fprintf(d,"predict value it=0");			 */
	      /* for(nNumClut=0;nNumClut<numDOF;++nNumClut){		 */
	      /* 	fprintf(d,"%e ",clust[nNumClut].ddihedang[0]);	 */
	      /* }							 */
	      /* fprintf(d,"\n");					 */
              /***********************************************************/


	      for (i=0;i<NUM_IETRATION;++i){
		pre_dyn();
		if (MODE==NVT) {
		  //CalcT(vel_Term);
		  //inertiaofvertialvariable2=2.0*Energy_kinetic8/q_NVT*s_NVT*(4.18407*100.0);
		  //acc_s_NVT=inertiaofvertialvariable2-(prot.DOF-1)*k_B_kcm*T_Kelvin/q_NVT*s_NVT*4.18407*100.0+dot_s_NVT*dot_s_NVT/s_NVT;
		  CalcT(vel_Term);
		  inertiaofvertialvariable2=2.0*Energy_kinetic8/q_NVT*s_NVT*(4.18407*100.0);

		  acc_s_NVT = inertiaofvertialvariable2-/*(prot.DOF-1)*//*degF*/DOFOFPROT*k_B_kcm*/*1.38065e-23*/T_Kelvin/q_NVT*s_NVT
		    *4.18407*100.0
		    /*/1.660539e-27/1.0e-20*1.0e-24*/
		    +dot_s_NVT*dot_s_NVT/s_NVT;
		}
		if (i==0)
		  cflag=ON;
		else
		  cflag=OFF;
		//  printf("yes3-1\n");
		calc_dd_theta_cycle(pflag,cflag,
				    ele, ALJ, BLJ,
				    p_e, p_1_4_e,
				    p_LJ,p_1_4_LJ,
				    f_e,f_1_4_e,
				    f_LJ,f_1_4_LJ,
				    numnb,indexnb,
				    num14,index14,eig,eig_14,eig_dihed);
		//  printf("yes3\n");
		for (nNumClut=0;nNumClut<numDOF;++nNumClut) {
                  vel_bh_dummy[nNumClut]=vel_bh[nNumClut]+deltat*clust[nNumClut].dddihedang[0]; // 0811
		  clust[nNumClut].ddihedang[0]=0.5*vel_bh[nNumClut]+0.5*vel_bh_dummy[nNumClut]; //0811
		  //		  clust[nNumClut].ddihedang[0]=vel_bh[nNumClut]+0.5*deltat*clust[nNumClut].dddihedang[0]+0.5*deltat*acc_bh[nNumClut];
		  if (MODE==NVT) {
		    dot_s_NVT_bh_dummy=dot_s_NVT_bh+deltat*acc_s_NVT;
		    //		    dot_s_NVT=dot_s_NVT_bh+0.5*deltat*acc_s_NVT;
		    dot_s_NVT=0.5*dot_s_NVT_bh+0.5*dot_s_NVT_bh_dummy;
		  }
		  if (TermMoveMode==ON) {
		    for (j=0;j<6;++j) {
		      vel_Term_bh_dummy[j]=vel_Term_bh[j]+deltat*acc_Term/*3*//*0911*/[j]; // 0811
		      vel_Term[j]=0.5*vel_Term_bh[j]+0.5*vel_Term_bh_dummy[j]; // 0811
		      //vel_Term[j]=vel_Term_bh[j]+0.5*deltat*(acc_Term[j]+acc_Term_bh[j]); // 0811
		    }
		  }
		}

		//                fprintf(d,"collect value it=%d",i);
		//		for(nNumClut=0;nNumClut<numDOF;++nNumClut){
		//		  fprintf(d,"%20.12e ",clust[nNumClut].ddihedang[0]);
		//		}
		//		fprintf(d,"\n");

	      }

	      //              fclose(d);

	      for (nNumClut=0;nNumClut<numDOF;++nNumClut) {
		vel_boh[nNumClut]=vel_bh[nNumClut]; // 0811
		vel_bh[nNumClut] = vel_bh_dummy[nNumClut]; //0811
		clust[nNumClut].now_deltadihedang[0]=deltat*vel_bh[nNumClut]; //0811
                /************************************************************************/
		//		clust[nNumClut].now_deltadihedang[0]=deltat*clust[nNumClut].ddihedang[0]+0.5*deltat*deltat*clust[nNumClut].dddihedang[0]; 0811
		//		vel_bh[nNumClut] = clust[nNumClut].ddihedang[0]; 0811
		//		acc_bh[nNumClut]=clust[nNumClut].dddihedang[0]; 0811
		if (TermMoveMode==ON) {
		  for (i=0;i<6;++i) {
		    vel_Term_boh[i]=vel_Term_bh[i]; // 0811
		    vel_Term_bh[i]=vel_Term_bh_dummy[i]; // 0811
		    //delta_Term[i]=deltat*vel_Term_bh[i]; // 0811
		    //                    delta_Term[i]=deltat*vel_Term[i]+0.5*deltat*deltat*acc_Term[i]; // 0811
		    //		    vel_Term_boh[i]=vel_bh[i]; // 0811
		    //		    vel_Term_bh[i]=vel_Term[i]; // 0811
		    //		    acc_Term_bh[i]=acc_Term[i]; // 0811
		  }
		}
	      }
	      if (MODE==NVT) {
		dot_s_NVT_boh=dot_s_NVT_bh;
		s_NVT+=deltat*dot_s_NVT_bh_dummy;
		dot_s_NVT_bh=dot_s_NVT_bh_dummy;
	      }
	    }
	    else if (GearMODE == ON) {
	      time = velocity_verlet_step1(nNumStep,pflag,
					   ele, ALJ, BLJ,
					   p_e, p_1_4_e, 
					   p_LJ,p_1_4_LJ,
					   f_e,f_1_4_e, 
					   f_LJ,f_1_4_LJ,
					   numnb, indexnb,
					   num14, index14);
	      //  printf("yes3\n");

	      if (TermMoveMode == ON) {
		for (i=0;i<6;++i) {
		  for (j=0;j<6;++j) {
		    predict_Term[i][j] = 0.0;
		    predict_Term3[i][j] = 0.0;
		  }
		}
		for (i=0;i<4;++i) {
		  for (j=0;j<5;++j) {
		    predict_q_Term[i][j]=0.0;
		  }
		}
		for (i=0;i<6;++i) {
		  for (j=0;j<6;++j) {
		    for (k=0;k<6;++k) {
		      predict_Term[i][j] += Telar_Matrix[j][k]*correct_Term[i][k];
		      predict_Term3[i][j] += Telar_Matrix[j][k]*correct_Term3[i][k];
		    }
		  }
		}
		for (i=0;i<4;++i) {
		  for (j=0;j<5;++j) {
		    for (k=0;k<5;++k) {
		      predict_q_Term[i][j] += Telar_Matrix[j][k]*correct_q_Term[i][k];
		    }
		  }
		}
		for (i=0;i<4;++i) {
		  q_Term[i]=predict_q_Term[i][0];
		  vel_q_Term[i]=predict_q_Term[i][1]/deltat;
		}
		len=0.0;
		for (i=0;i<4;++i) {
		  len+=q_Term[i]*q_Term[i];
		}
		for (i=0;i<4;++i) {
		  q_Term[i]=q_Term[i]/len;
		}


		for (i=1;i<6;++i) {
		  delta_Term[i]= 0.0;
		}
		for (i=0;i<6;++i) {
		  for (j=1;j<6;++j){
		    delta_Term[i]+= Telar_Matrix[0][j]*correct_Term3[i][j];
		  }
		  vel_Term[i]=predict_Term[i][1]/deltat;
		}

	      }
	    }
	  }
//BD///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
	  else if (DYNMODE == BD) {
	    time = verlet(nNumStep);
	  }

	  if (GearMODE == ON ) {
	    for(nNumClut=/*1*//*0*/1; nNumClut<prot.DOF; ++nNumClut) {
	      if (clust[nNumClut].num_branch > 1) {
		nNumClutOrigBranch = nNumClut;
	      }
	      trans_CN_to_A(nNumClut, nNumClutOrigBranch, q_Term);
	    }
	  }
	  //  printf("yes4\n");
		
	  if(GearMODE == ON) {
	    for (i=0;i<NUM_IETRATION;++i){
	      pre_dyn();
	      //  printf("yes4-2\n");
	      velocity_verlet_step2(pflag,
				    ele, ALJ, BLJ,
				    p_e, p_1_4_e, 
				    p_LJ,p_1_4_LJ,
				    f_e,f_1_4_e, 
				    f_LJ,f_1_4_LJ,
				    numnb, indexnb,
				    num14, index14,
				    vel_Term,eig,eig_14,eig_dihed);
	      //  printf("yes5\n");

	      if (TermMoveMode == ON) {
		for (j=0;j<6;++j) {
		  delta_acc_Term[j] = 0.5*deltat*deltat*acc_Term[j]-predict_Term[j][2];
		  delta_acc_Term3[j] = 0.5*deltat*deltat*acc_Term2/*3*0911*/[j]-predict_Term3[j][2];
		}
		for (j=0;j<6;++j) {
		  for (k=0;k<6;++k) {
		    correct_Term[j][k] = 0.0;
		    correct_Term3[j][k] = 0.0;
		  }
		}
		for (j=0;j<4;++j) {
		  for (k=0;k<5;++k) {
		    correct_q_Term[j][k] = 0.0;
		  }
		}
		for (j=0;j<6;++j) {
		  for (k=0;k<6;++k) {
		    correct_Term[j][k] = predict_Term[j][k]+GearsConstant[k]*delta_acc_Term[j];
		    correct_Term3[j][k] = predict_Term3[j][k]+GearsConstant[k]*delta_acc_Term3[j];
		  }
		}    

		for (j=0;j<6;++j) {
		  delta_Term[j]=GearsConstant[0]*delta_acc_Term3[j];
		  vel_Term[j]=correct_Term[j][1]/deltat;
		}

		vel_q_Term[3]=0.5*(-q_Term[2]*vel_Term[0]-q_Term[0]*vel_Term[1]+q_Term[1]*vel_Term[2]);
		vel_q_Term[1]=0.5*( q_Term[0]*vel_Term[0]-q_Term[2]*vel_Term[1]-q_Term[3]*vel_Term[2]);
		vel_q_Term[2]=0.5*( q_Term[3]*vel_Term[0]+q_Term[1]*vel_Term[1]+q_Term[0]*vel_Term[2]);
		vel_q_Term[0]=0.5*(-q_Term[1]*vel_Term[0]+q_Term[3]*vel_Term[1]-q_Term[2]*vel_Term[2]);

		for (j=0;j<5;++j) {
		  delta_q_Term[j] = deltat*vel_q_Term[j]-predict_q_Term[j][1];
		}
		for (j=0;j<4;++j) {
		  for (k=0;k<5;++k) {
		    correct_q_Term[j][k] = predict_q_Term[j][k]+GearsConstant3[k]*delta_q_Term[j];
		  }
		  q_Term[j]=correct_q_Term[j][0];
		  vel_q_Term[j]=correct_q_Term[j][1]/deltat;
		}    
		len=0.0;
		for (i=0;i<4;++i) {
		  len+=q_Term[i]*q_Term[i];
		}
		for (i=0;i<4;++i) {
		  q_Term[i]=q_Term[i]/len;
		}


	      }

	      // ¼¡¤Î¥¹¥Æ¥Ã¥×¤Î 2 ÌÌ³Ñ¤ò»ý¤Ã¤¿¹½Â¤¤ò·×»»
	      for(nNumClut=/*1*/0; nNumClut<prot.DOF; ++nNumClut) {
		if (clust[nNumClut].num_branch > 1) {
		  nNumClutOrigBranch = nNumClut;
		}
		trans_CN_to_A(nNumClut, nNumClutOrigBranch,q_Term);
	      }
	    }
	  }
	  else if (GearMODE == OFF) {
	    // ¼¡¤Î¥¹¥Æ¥Ã¥×¤Î 2 ÌÌ³Ñ¤ò»ý¤Ã¤¿¹½Â¤¤ò·×»»
	    for(nNumClut=1; nNumClut<prot.DOF; ++nNumClut) {
	      if (clust[nNumClut].num_branch > 1) {
		nNumClutOrigBranch = nNumClut;
	      }
	      trans_CN_to_A(nNumClut, nNumClutOrigBranch,q_Term);
	    }
	  }
	    

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
	  // °Ê²¼¡¢½ÐÎÏ ////////////////////////////////////////////////////////
	  //////////////////////////////////////////////////////////////////////

	  // ¥¹¥Æ¥Ã¥×¿ô¤Î½ÐÎÏ
	  //			printf("nstep = %d\n",time);

	  if ((nNumStep % out_put_steps) == 0)  {
	    
	    // ½ÐÎÏ¥Õ¥¡¥¤¥ë¤Î¥ª¡¼¥×¥ó_1
	    if ((output_c = fopen(/*"coo_pro.out"*/OutfilCOORD,"a")) == NULL) {
	      printf("error coo_pro.out cannot open\n");
	      exit(1);
	    }

	    // Ãà¼¡¹½Â¤¤Î out ¥Õ¥¡¥¤¥ë½ÐÎÏ
	    output_file_coo_pro(output_c);
	    /*****************************/
	    /* printf("m:776\n");	       */
	    /*****************************/
	    
	    fclose(output_c);

	    if (veloutflag == ON || veloutflag == 3) {
	      if ((output_v = fopen(/*"velo_pro.out"*/OutfilVELO,"a")) == NULL) {
		printf("error velo_pro.out cannot open\n");
		exit(1);
	      }
		      
	      if (veloutflag==ON) {
		// Ãà¼¡Â®ÅÙ¤Î out ¥Õ¥¡¥¤¥ë½ÐÎÏ
		output_file_vel_pro(output_v);
	      }
	      else if (veloutflag==3) {
		output_file_dihed_vel_pro(output_v);
	      }
	      fclose(output_v);
	    }
	  }
		  
	  
	  
	}
	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////

	// ½ÐÎÏ¥Õ¥¡¥¤¥ë¤Î¥¯¥í¡¼¥º
	///////////////////////////////////////////////////////
	// ½ÐÎÏ¥Õ¥¡¥¤¥ë¤Î¥ª¡¼¥×¥ó_1
	if ((output_c = fopen(/*"coo_pro.out"*/OutfilCOORD,"a")) == NULL) {
	  printf("error coo_pro.out cannot open\n");
	  exit(1);
	}
	fprintf(output_c,"111111");
	fclose(output_c);
	///////////////////////////////////////////////////////

	if ((rst_c = fopen(/*"crd_rst.in"*/OutfilRESTAC,"w")) == NULL) {
	  printf("cnt open crd_rst.in\n");
	  exit(1);
	}
	fprintf(rst_c,"\n");
	fprintf(rst_c,"%d\n",prot.num_atom);
	for (i=0;i<prot.num_atom;++i) {
	  fprintf(rst_c,"%12.8lf %12.8lf %12.8lf\n",prot.coord[i][0],prot.coord[i][1],prot.coord[i][2]);
	}
	fclose(rst_c);
	
	if ((rst_v = fopen("velo.in"/*OutfilRESTAV*/,"w")) == NULL) {
	  printf("cnt open velo.in\n");
	  exit(1);
	}
	if (TermMoveMode2==12) {
	  for (i=0;i<6;++i) {
	    fprintf(rst_v,"%lf\n",vel_Term[i]);
	  }
	  for (i=1;i<prot.DOF;++i) {
	    fprintf(rst_v,"%lf\n",clust[i].ddihedang[0]);
	  }
	}
	else {
	  for (i=0;i<prot.DOF;++i) {
	    fprintf(rst_v,"%lf\n",clust[i].ddihedang[0]);
	  }
	}
	if (GearMODE==ON) {
	  fprintf(rst_v,"%12.8lf\n",s_NVT_correct[0]);
	  fprintf(rst_v,"%12.8lf\n",s_NVT_correct[1]);
	}
	else {
	  fprintf(rst_v,"%12.8lf\n",s_NVT);
	  fprintf(rst_v,"%12.8lf\n",dot_s_NVT*deltat);
	}
	fclose(rst_v);

	coord_rst=(double *)calloc(sizeof(double),prot.num_atom*3);
	for (i=0;i<prot.num_atom;++i)
	  for (j=0;j<3;++j)
	    for (k=0;k<3;++k)
	      coord_rst[i*3+j]+=clust[0].trans_A_to_CN[0][j][k]*prot.coord[i][k];
	if ((rst_a = fopen(OutfilRESTAMFORM,"w")) == NULL) {
	    printf("cnt open velo.in\n");
	    exit(1);
	}
	fprintf(rst_a,"\n");
	fprintf(rst_a,"%5d  0.1000000E+04\n",prot.num_atom);
	k=0;
	for (i=0;i<prot.num_atom;++i) {
	  for (j=0;j<3;++j) {
	    ++k;
	    fprintf(rst_a,"%12.7lf",coord_rst[i*3+j]);
	    if (k ==6 ) {
	      fprintf(rst_a,"\n");
	      k=0;
	    }
	  }
	}
	if (k !=6 )
	  fprintf(rst_a,"\n");
	k=0;
	for (i=0;i<prot.num_atom;++i) {
	  for (j=0;j<3;++j) {
	    ++k;
	    fprintf(rst_a,"%12.7lf",prot.velo[i][j]/20.455);
	      if (k ==6 ) {
		fprintf(rst_a,"\n");
		k=0;
	      }
	  }
	}
	fclose(rst_a);


	/******************************/ // 0811
        /* free(ele);		      */
	/* free(ALJ);		      */
	/* free(BLJ);		      */
	/* free(p_e);		      */
	/* free(p_1_4_e);	      */
	/* free(p_LJ);		      */
	/* free(p_1_4_LJ);	      */
	/* free(f_e); 		      */
	/* free(f_1_4_e); 	      */
	/* free(f_LJ);		      */
	/* free(f_1_4_LJ);	      */
	/* free(vel);		      */
	/* free(vel_bh);	      */
	/* free(vel_bh_dummy);	      */
	/* free(vel_boh);	      */
        /******************************/ // 0811
	
	return 0;
}

void pre_dyn(void) {
  int i;
  int nNumClut;
  int alpha;
  int nNumClutOrigBranch2;

  /***************************/
  Initialize_Variables();
  /***************************/

  // ÊÑ´¹¹ÔÎó¤ÎÀßÄê
  for(nNumClut=0; nNumClut<prot.DOF; ++nNumClut){
    if (clust[nNumClut].num_branch > 1){
      nNumClutOrigBranch2 = nNumClut;
    }
    set_trans_Matrix(nNumClut, nNumClutOrigBranch2);
  }

  // sptial velocity¡¢corioli ²ÃÂ®ÅÙ¡¢corioli ÎÏ¤Î·×»»
  for(nNumClut=0; nNumClut<prot.DOF; ++nNumClut){
    if (clust[nNumClut].num_branch > 1) {
      nNumClutOrigBranch2 = nNumClut;
    }
    set_sp_velo(nNumClut, nNumClutOrigBranch2);
  }

  // sptial velocity¡¢corioli ²ÃÂ®ÅÙ¡¢corioli ÎÏ¤Î·×»»
  for(nNumClut=0; nNumClut<prot.DOF; ++nNumClut) {
    // corioli ²ÃÂ®ÅÙ
    set_coriolis_acc(nNumClut);
    // corioli ÎÏ
    set_coriolis_force(nNumClut);
  }
  scaling_momentum(0);


}

void scaling_momentum(int nNumStep) {
  int i,j,nNumClut;
  double momentum[6];
  double sum;

  for(nNumClut=0; nNumClut<prot.DOF; ++nNumClut){
    for (i=0; i<6; ++i) {
      for (j=0; j<6; ++j) {
	momentum[i] += clust[nNumClut].InertiaMatrix[i][j]*clust[nNumClut].sp_velo[j];
      }
    }
  }

}

