#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "netcdf.h"
#include <getopt.h>

#include "ABAb.h"  // 2014-06-18

#include "PTL.h"
#include "FFL.h"
#include "EF.h"

#include "ECEPE.h" // 2015-02-16

#include "netcdf_mineL.h"

#include "f2c.h"     // 2014-08-13
#include "clapack.h" // 2014-08-13

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,d;
  int numatom,numclut,numstep=100000;
  int interval=1000,intervalout=1000,intervalnc=1000,intervalflag;
  double dt=0.001; double pot_KE_KEv_PEv=0.0; // 2014-08-13
  CLT *clt;
  double *Q,*frc,pot;
  double *q;
  double *qacc,*qvel,*qrot;
  double *predict,*correct;
  double KE,KEv,PEv;
  double dummy;

  int MODE=NVT,MODEV=OFF,TERMMODE=OFF;
  double Tobj=300,KEobj;
  double k_B=1.98723e-3;

  double s=1.0,s_vel=0.0,s_acc,gzi=0.0,gzi_vel,predict_gzi[5],correct_gzi[5],predict_s[5],correct_s[5];
  double Q_NH,tau=0.1,tau2;
  double T;
  int DOF;

  int dfalg=ON,efalg=ON,LJfalg=ON,e14falg=ON,LJ14falg=ON;
  double coff=1.0;
  int *numclutparent,*terminal,*origin;
  double *trans_A_to_CN_terminal/*[3][3]*/,*l_Term/*[3]*/; // 2014-06-29

  double /**delta_Term, 2014-08-13*/*vel_Term,*acc_Term,*acc_Term2,**predict_Term,**predict_Term2,/***correct_Term,2014-08-13*/**correct_Term2;
  double **correct_Term_2014_08_13; // 2014-08-13
  //  double delta_Term[6]; // 2014-08-13
  double *delta_Term; // 2014-08-13
  //  double correct_Term[6][6]; // 2014-08-13

  double /***crd_nc,*/*energy;
  double crd_nc[MAXATOM][3]; // 2014-08-13 CHECK
  struct my_netcdf_out_id_AMBER nc_id;

  double *crd,*mass,*crd_Term/*2014-06-30*/;
  struct potential e;
  struct force f;

  char *line="line"; // 2014-08-13
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double maxabs_force,maxabs_Q;
  int index_atom,index_clust;
  double maxabs_force_trj,maxabs_Q_trj;
  int index_atom_trj,index_clust_trj;

  char *inputfilename,*reffilename,*outputfilename,*outputfilename2,*parmfilename,*clustfilename,*inputvelofilename;
  char *trjfilename;
  char *rstfilename="rstfile",*rstvelfilename="rstvelfile";
  FILE *inputfile,*reffile,*inputfilegoal,*outputfile,*outputfile2,*parmfile,*clustfile,*inputvelofile;
  FILE *velfile,*rstfile,*rstvelfile;

  char *preofilename;                  // 2015-02-16
  char *bd8filename;                   // 2015-02-16
  char *coofilename;                   // 2015-02-16
  char *coofilename_for_sflag;         // 2015-02-16
  char *massfilename;                  // 2015-02-16
  char *angfilename;                   // 2015-02-16
  char *pairsfilename;                 // 2015-02-16

  int x;
  int *IndexOfABICycle; // 2014-07-02

  int joinflag; // 2014-08-13

  double *qvel_b1,*qvel_b2;           // 2015-01-19
  double *qvel_b1_temp,*qvel_b2_temp; // 2015-01-19

  double *vel_Term_b1,*vel_Term_b2;           // 2015-01-19
  double *vel_Term_b1_temp,*vel_Term_b2_temp; // 2015-01-19

  double zeta=0.0;                            // 2015-01-19

  double KE_1_2,KEv_1_2,PEv_1_2,T_1_2;        // 2015-01-19

  int leapfrog_flag=OFF;              // 2015-01-19

  int AmberON=ON;                           // 2015-02-16

  // for ECEPE                               // 2015-02-16
  double *dihed;                             // 2015-02-16
  double *delta_dihed;                       // 2015-02-16
  //  double *p_d,*p_e, *p_LJ;                   // 2015-02-16
  double p_d_t,p_e_t,p_LJ_t;                 // 2015-02-16
  double *f_e, *f_LJ;                        // 2015-02-16
  int **pairs1_5,**pairs1_4,*num1_5,*num1_4; // 2015-02-16
  double *charge;                            // 2015-02-16
  int *nbtype;                               // 2015-02-16
  struct ECEPE_parms ECEPE_p;                // 2015-02-16
  struct pnb  nb_p;                          // 2015-02-16
  // for ECEPE                               // 2015-02-16
  //  int oliflag=OFF;                           // 2015-02-16

  flagnb14=ON/*OFF*/;    // 2015-02-16
  flagnb15=ON/*OFF*/;    // 2015-02-16
  flages14=ON/*OFF*/;    // 2015-02-16
  flages15=ON/*OFF*/;    // 2015-02-16
  flagang=OFF;     // 2015-02-16
  flagpairs=OFF;   // 2015-02-16

  preofilename = "preo.in";              // 2015-02-16
  bd8filename  =  "bd8";                 // 2015-02-16
  coofilename  =  "coo.in";              // 2015-02-16
  angfilename  =  "ang.in";              // 2015-02-16
  coofilename_for_sflag = "coo_s.in";    // 2015-02-16
  massfilename =  "mass.in";             // 2015-02-16
  pairsfilename = "pairs.in";            // 2015-02-16

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"f",0,NULL,'f'},
    {"nve",0,NULL,'*'},
    {"termon",0,NULL,'+'},
    {"dih",0,NULL,'d'},
    {"els",0,NULL,'e'},
    {"lj",0,NULL,'l'},
    {"e14",0,NULL,'1'},
    {"l14",0,NULL,'4'},
    {"h",0,NULL,'h'},
    {"kc",1,NULL,'k'},
    {"cof",1,NULL,'c'},
    {"temp",1,NULL,'t'},
    {"dt",1,NULL,'@'},
    {"omg",1,NULL,'o'},
    {"nums",1,NULL,'s'},
    {"int",1,NULL,'i'},
    {"intnc",1,NULL,'j'},
    {"intout",1,NULL,'k'},
    {"vel",1,NULL,'v'},
    {"tau",1,NULL,'?'},
    {"rst",1,NULL,'{'},
    {"rstv",1,NULL,'}'},
    {"leapfrog",0,NULL,'L'},       // 2015-01-19
    {"predictcorrect",0,NULL,'P'}, // 2015-01-19
    {"ECEPE",0,NULL,'E'}, // 2015-02-16
    {"preo",1,NULL,'p'}, // 2015-02-16
    {"bd8",1,NULL,'b'}, // 2015-02-16
    {"coo",1,NULL,'2'}, // 2015-02-16
    {"angfil",1,NULL,'3'}, // 2015-02-16
    {"coosf",1,NULL,'u'}, // 2015-02-16
    {"pairs",1,NULL,'q'}, // 2015-02-16
    {"mass",1,NULL,'w'}, // 2015-02-16
    //    {"oliflagON",0,NULL,'0'}, // 2015-02-16
    {0,0,0,0}
  };

  /***************************************************************************/
  /* static long int m=3,n=3,lda=3,info,piv[3],lwork=3; // 2014-08-13	     */
  /* static double work[3];			     // 2014-08-13	     */
  /* double A[9];                                       // 2014-08-13	     */
  /* A[0]=1.;A[1]=2.;A[2]=3.;			     // 2014-08-13	     */
  /* A[3]=4.;A[4]=9.;A[5]=14.;			     // 2014-08-13	     */
  /* A[6]=6.;A[7]=14.;A[8]=23.;			     // 2014-08-13	     */
  /* dgetrf_( &m, &n, A, &lda, piv, &info);	     // 2014-08-13	     */
  /* dgetri_(&n, A, &lda, piv, work, &lwork, &info);    // 2014-08-13	     */
  /* //  invm3(A,B);                                    // 2014-08-13	     */
  /***************************************************************************/

  while((c=getopt_long(argc,argv,"*+del14h@:t:k:c:t:o:s:i:j:k:v:?:{:}:L:P:E:p:b:2:3:u:q:w:",long_opt,&opt_idx))!=-1) { // 2015-02-16
    switch(c) {
    case '*':
      MODE=NVE;
      break;
    case '+':
      TERMMODE=ON;
      break;
    case 'd':
      dfalg=OFF;
      break;
    case 'e':
      efalg=OFF;
      break;
    case 'l':
      LJfalg=OFF;
      break;
    case '1':
      e14falg=OFF;
      break;
    case '4':
      LJ14falg=OFF;
      break;
    case 't':
      Tobj=atof(optarg);
      break;
    case 'o':
      tau=atof(optarg);
      break;
    case '@':
      dt=atof(optarg);
      break;
    case 'i':
      interval=atoi(optarg);
      break;
    case 'j':
      intervalnc=atoi(optarg);
      break;
    case 'k':
      intervalout=atoi(optarg);
      break;
    case 's':
      numstep=atoi(optarg);
      break;
    case 'c':
      coff=atof(optarg);
      break;
    case 'v':
      inputvelofilename=optarg;
      MODEV=ON;
      break;
    case '?':
      tau=atof(optarg);
      break;
    case '{':
      rstfilename=optarg;
      break;
    case '}':
      rstvelfilename=optarg;
      break;
    case 'L':
      leapfrog_flag=ON;
      break;
    case 'P':
      leapfrog_flag=OFF;
      break;
    case 'E':            // 2015-02-16
      AmberON=OFF;       // 2015-02-16
      break;             // 2015-02-16
    case 'p':              // 2015-02-16
      preofilename=optarg; // 2015-02-16
      break;               // 2015-02-16
    case 'b':              // 2015-02-16
      bd8filename=optarg;  // 2015-02-16
      break;               // 2015-02-16
    case '2':              // 2015-02-16
      coofilename=optarg;  // 2015-02-16
      break;               // 2015-02-16
    case '3':              // 2015-02-16
      flagang=ON;          // 2015-02-16
      angfilename=optarg;  // 2015-02-16
      break;               // 2015-02-16
    case 'u':              // 2015-02-16
      coofilename_for_sflag=optarg; // 2015-02-16
      break;               // 2015-02-16
    case 'q':               // 2015-02-16
      flagpairs=ON;         // 2015-02-16
      pairsfilename=optarg; // 2015-02-16
      break;                // 2015-02-16
    case 'w':               // 2015-02-16
      massfilename=optarg;  // 2015-02-16
      break;                // 2015-02-16
      //    case '0':               // 2015-02-16
      //      oliflag=ON;           // 2015-02-16
      //      break;                // 2015-02-16
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (AmberON==ON) {                // 2015-02-16
    if (argc < 6) {
      USAGE(progname);
      exit(1);
    }
    inputfilename     = *argv;
    clustfilename     = *++argv;
    parmfilename      = *++argv;
    outputfilename    = *++argv;
    outputfilename2   = *++argv;
    trjfilename       = *++argv;
  }                                // 2015-02-16
  else {                           // 2015-02-16
    if (argc < 6) {               // 2015-02-16
      USAGE(progname);             // 2015-02-16
      exit(1);                     // 2015-02-16
    }                              // 2015-02-16
    inputfilename     = *argv;     // 2015-02-16
    clustfilename     = *++argv;   // 2015-02-16
    parmfilename      = *++argv;   // 2015-02-16
    outputfilename    = *++argv;   // 2015-02-16
    outputfilename2   = *++argv;   // 2015-02-16
    trjfilename       = *++argv;   // 2015-02-16

  }

  //  correct_Term_2014_08_13=(double **)emalloc(sizeof(double *)*6);   // 2014-08-13 // 2014-09-05
  correct_Term_2014_08_13=(double **)calloc(6,sizeof(double *));   // 2014-09-05
  //  correct_Term_2014_08_13[1]=(double *)emalloc(sizeof(double)*6);   // 2014-08-13
  //  correct_Term_2014_08_13[0]=(double *)emalloc(sizeof(double)*6);   // 2014-08-13
  for (i=0/*2*/;i<6;++i) {                                    // 2014-08-13
    //    correct_Term_2014_08_13[i]=(double *)emalloc(sizeof(double)*6); // 2014-08-13 // 2014-09-05
    correct_Term_2014_08_13[i]=(double *)calloc(6,sizeof(double)); // 2014-09-05
  }                                             // 2014-08-13

  if (AmberON==ON) {                                                                  // 2015-02-16
    parmfile=efopen(parmfilename,"r");                                                
    readParmtopL(parmfile);                                                           
    fclose(parmfile);                                                                 
    //  exit(1); // 2014-07-15                                                        
    numatom=AP.NATOM;                                                                 
    //  mass=(double *)gcemalloc(sizeof(double)*numatom); // 2014-07-22               
    //  mass=(double *)emalloc(sizeof(double)*numatom); // 2014-07-22 // 2014-09-05   
    mass=(double *)calloc(numatom,sizeof(double)); // 2014-09-05                      
    for (i=0;i<numatom;++i) mass[i]=AP.AMASS[i];                                      
    //  crd=(double *)gcemalloc(sizeof(double)*numatom*3); // 2014-07-22
    //  crd=(double *)emalloc(sizeof(double)*numatom*3); // 2014-07-22 // 2014-09-05
    crd=(double *)calloc(numatom*3,sizeof(double)); // 2014-09-05
    //  crd_Term=(double *)gcemalloc(sizeof(double)*numatom*3); // 2014-06-30 // 2014-07-22
    //  crd_Term=(double *)emalloc(sizeof(double)*numatom*3); // 2014-07-22 // 2014-09-05
    crd_Term=(double *)calloc(numatom*3,sizeof(double)); // 2014-07-22 // 2014-09-05
    inputfile=efopen(inputfilename,"r");
    getline(&line,&len,inputfile);
    fscanf(inputfile,"%d",&d);
    for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
    fclose(inputfile);

  }                                                                                   // 2015-02-16
  else {                                                                              // 2015-02-16
    read_ECEPE_parm(preofilename,bd8filename,&ECEPE_p,&nb_p);                         // 2015-02-16
    numatom=ECEPE_p.NUMATM;                                                           // 2015-02-16
    pairs1_5=(int **)calloc((numatom-1)*2,sizeof(int *));                             // 2015-02-16
    pairs1_4=(int **)calloc((numatom-1)*2,sizeof(int *));                             // 2015-02-16
    num1_4=(int *)calloc(numatom,sizeof(int ));                                       // 2015-02-16
    num1_5=(int *)calloc(numatom,sizeof(int ));                                       // 2015-02-16
    dihed=(double *)calloc((ECEPE_p.NUMVAR),sizeof(double));                          // 2015-02-16
    delta_dihed=(double *)calloc((ECEPE_p.NUMVAR),sizeof(double));                    // 2015-02-16
    mass=(double *)calloc(numatom,sizeof(double));                                    // 2015-02-16

    crd=(double *)calloc(numatom*3,sizeof(double));                                   // 2015-02-16
    crd_Term=(double *)calloc(numatom*3,sizeof(double));                              // 2015-02-16

    pick_ECEPP_data(preofilename,bd8filename,coofilename,coofilename_for_sflag,       // 2015-02-16
		    pairs1_5,pairs1_4,num1_4,num1_5,                                  // 2015-02-16
		    &ECEPE_p,&nb_p/*,delta_dihed*/,                                   // 2015-02-16
		    angfilename,flagang,flagpairs,pairsfilename,delta_dihed,dihed,crd,// 2015-02-16
		    inputfilename,parmfilename/*,oliflag*/);                              // 2015-02-16
    pick_mass_data(massfilename,ECEPE_p,numclut,mass);                                // 2015-02-16
  }                                                                                   // 2015-02-16

  clustfile=efopen(clustfilename,"r");
  //  clt=ABAp_clustscan(clustfile,&numclut); // 2014-06-20
  //  ABAp_clustscan(clustfile,&numclut,clt); // 2014-06-19
  // 2014-07-01
  //  int i,j,x;
  //  CLT *clt; // 2014-06-20

  fscanf(clustfile,"%d",&numclut);

  //  clt=(CLT *)gcemalloc(sizeof(CLT)*(numclut)); // 2014-07-22
  //  clt=(CLT *)emalloc(sizeof(CLT)*(numclut)); // 2014-07-22 // 2014-09-05
  clt=(CLT *)calloc(numclut,sizeof(CLT)); // 2014-09-05

  for(i=0;i<(numclut);++i) fscanf(clustfile,"%d",&clt[i].origin_atom_a);
  for(i=0;i<(numclut);++i) fscanf(clustfile,"%d",&clt[i].terminal);
  for(i=0;i<(numclut);++i) fscanf(clustfile,"%d",&clt[i].num_atom_clust);
  for(i=0;i<(numclut);++i) fscanf(clustfile,"%d",&clt[i].num_branch);
  for(i=0;i<(numclut);++i) fscanf(clustfile,"%d",&x);
  for(i=0;i<(numclut);++i) 
    for(j=0;j<clt[i].num_branch;++j)
      fscanf(clustfile,"%d",&clt[i].terminal_atom_a[j]);
  for (i=0;i<(numclut);++i) fscanf(clustfile, "%d", &clt[i].nNumClutOfParent);
  for (i=0;i<(numclut);++i)
    for(j=0;j<clt[i].num_branch;++j)
      fscanf(clustfile, "%d", &clt[i].nNumClutOfChild[j]);
  //  IndexOfABICycle=(int *)gcemalloc(sizeof(int)*(numclut)); // 2014-07-22
  //  IndexOfABICycle=(int *)emalloc(sizeof(int)*(numclut)); // 2014-07-22 // 2014-09-05
  IndexOfABICycle=(int *)calloc(numclut,sizeof(int)); // 2014-07-22 // 2014-09-05
  for (i=0;i<(numclut);++i) fscanf(clustfile, "%d", &IndexOfABICycle[i]);

  for(i=0;i<(numclut);++i) {
    //    clt[i].frc_e=(double *)gcemalloc(sizeof(double)*clt[i].num_atom_clust*3); // 2014-07-22
    //    clt[i].frc_e=(double *)emalloc(sizeof(double)*clt[i].num_atom_clust*3); // 2014-07-22 // 2014-09-05
    clt[i].frc_e=(double *)calloc(clt[i].num_atom_clust*3,sizeof(double)); // 2014-07-22 // 2014-09-05
    //    clt[i].frc_LJ=(double *)gcemalloc(sizeof(double)*clt[i].num_atom_clust*3); // 2014-07-22
    //    clt[i].frc_LJ=(double *)emalloc(sizeof(double)*clt[i].num_atom_clust*3); // 2014-07-22 // 2014-09-05
    clt[i].frc_LJ=(double *)calloc(clt[i].num_atom_clust*3,sizeof(double)); // 2014-07-22 // 2014-09-05
    //    clt[i].frc_1_4_e=(double *)gcemalloc(sizeof(double)*clt[i].num_atom_clust*3); // 2014-07-22
    //    clt[i].frc_1_4_e=(double *)emalloc(sizeof(double)*clt[i].num_atom_clust*3); // 2014-07-22 // 2014-09-05
    clt[i].frc_1_4_e=(double *)calloc(clt[i].num_atom_clust*3,sizeof(double)); // 2014-07-22 // 2014-09-05
    //    clt[i].frc_1_4_LJ=(double *)gcemalloc(sizeof(double)*clt[i].num_atom_clust*3); // 2014-07-22
    //    clt[i].frc_1_4_LJ=(double *)emalloc(sizeof(double)*clt[i].num_atom_clust*3); // 2014-07-22 // 2014-09-05
    clt[i].frc_1_4_LJ=(double *)calloc(clt[i].num_atom_clust*3,sizeof(double)); // 2014-07-22 // 2014-09-05
  }
  // 2014-07-01
  fclose(clustfile);

  //  numclutparent=(int *)gcemalloc(sizeof(int)*numclut); // 2014-07-22
  //  numclutparent=(int *)emalloc(sizeof(int)*numclut); // 2014-07-22 // 2014-09-05
  numclutparent=(int *)calloc(numclut,sizeof(int)); // 2014-07-22 // 2014-09-05
  //  terminal=(int *)gcemalloc(sizeof(int)*numclut); // 2014-07-22
  //  terminal=(int *)emalloc(sizeof(int)*numclut); // 2014-07-22 // 2014-09-05
  terminal=(int *)calloc(numclut,sizeof(int)); // 2014-07-22 // 2014-09-05
  //  origin=(int *)gcemalloc(sizeof(int)*numclut); // 2014-07-22
  //  origin=(int *)emalloc(sizeof(int)*numclut); // 2014-07-22 // 2014-09-05
  origin=(int *)calloc(numclut,sizeof(int)); // 2014-07-22 // 2014-09-05
  for (i=0;i<numclut;++i) {
    numclutparent[i]=clt[i].nNumClutOfParent;
    terminal[i]=clt[i].terminal_atom_a[0];
    origin[i]=clt[i].origin_atom_a;
  }

  /********************************************************************/
  /* ABAs_local_reference(clt,numclut,numatom,crd);     // 2014-08-13 */
  /* ABAs_trans_Matrix(clt,numclut,numatom,crd);        // 2014-08-13 */
  /* ABAs_inertia_matrix(clt,numclut,numatom,crd,mass); // 2014-08-13 */
  /********************************************************************/

  //  joinflag=0;                                            // 2014-08-13 // 2014-09-16
  //  for(i=0; i<numclut; ++i) {                             // 2014-08-13 // 2014-09-16
  //    joinflag=ABA_setJoin(clt,i,joinflag);                // 2014-08-13 // 2014-09-16
  //  }                                                      // 2014-08-13 // 2014-09-16
  ABAs_local_reference(clt,numclut,numatom,crd);         // 2014-08-13
  ABAs_trans_Matrix(clt,numclut,numatom,crd);            // 2014-08-13
  ABAs_inertia_matrix(clt,numclut,numatom,crd,mass);     // 2014-08-13

  //  Q=(double *)gcemalloc(sizeof(double)*numclut); // 2014-07-22
  //  Q=(double *)emalloc(sizeof(double)*numclut); // 2014-07-22 // 2014-09-05
  Q=(double *)calloc(numclut,sizeof(double)); // 2014-07-22 // 2014-09-05
  //  frc=(double *)gcemalloc(sizeof(double)*numatom*3); // 2014-07-22
  //  frc=(double *)emalloc(sizeof(double)*numatom*3); // 2014-07-22 // 2014-09-05
  frc=(double *)calloc(numatom*3,sizeof(double)); // 2014-07-22 // 2014-09-05

  if (AmberON==ON) {                                                 // 2015-02-16
    ffL_set_calcffandforce(&e,&f);
  }                                                                  // 2015-02-16
  else {                                                             // 2015-02-16
    charge=(double *)calloc(numatom,sizeof(double));                 // 2015-02-16
    nbtype=(int *)calloc(numatom,sizeof(int));                       // 2015-02-16
    for (i=0;i<numatom;++i) {                                        // 2015-02-16
      charge[i]=ECEPE_p.atom[i].charge;                              // 2015-02-16
      nbtype[i]=ECEPE_p.atom[i].nbtype;                              // 2015-02-16
    }                                                                // 2015-02-16
  }                                                                  // 2015-02-16

  //  q=(double *)gcemalloc(sizeof(double)*numclut); // 2014-07-22
  //  q=(double *)emalloc(sizeof(double)*numclut); // 2014-07-22
  q=(double *)calloc(numclut,sizeof(double)); // 2014-09-05
  //  qacc=(double *)gcemalloc(sizeof(double)*numclut); // 2014-07-22
  //  qacc=(double *)emalloc(sizeof(double)*numclut); // 2014-07-22 // 2014-09-05
  qacc=(double *)calloc(numclut,sizeof(double)); // 2014-07-22 // 2014-09-05
  //  qvel=(double *)gcemalloc(sizeof(double)*numclut); // 2014-07-22
  //  qvel=(double *)emalloc(sizeof(double)*numclut); // 2014-07-22 // 2014-09-05
  qvel=(double *)calloc(numclut,sizeof(double)); // 2014-07-22 // 2014-09-05
  //  qrot=(double *)gcemalloc(sizeof(double)*numclut); // 2014-07-22
  //  qrot=(double *)emalloc(sizeof(double)*numclut); // 2014-07-22 // 2014-09-05
  qrot=(double *)calloc(numclut,sizeof(double)); // 2014-07-22 // 2014-09-05
  //  predict=(double *)gcemalloc(sizeof(double)*numclut*6/**6*/); // 2014-07-22
  //  predict=(double *)emalloc(sizeof(double)*numclut*6/**6*/); // 2014-07-22 // 2014-09-05
  predict=(double *)calloc(numclut*6/**6*/,sizeof(double)); // 2014-07-22 // 2014-09-05
  //  correct=(double *)gcemalloc(sizeof(double)*numclut*6/**6*/); // 2014-07-22
  //  correct=(double *)emalloc(sizeof(double)*numclut*6/**6*/); // 2014-07-22 // 2014-09-05
  correct=(double *)calloc(numclut*6/**6*/,sizeof(double)); // 2014-07-22 // 2014-09-05

  qvel_b1=(double *)calloc(numclut*6,sizeof(double));       // 2015-01-19
  qvel_b2=(double *)calloc(numclut*6,sizeof(double));       // 2015-01-19
  qvel_b1_temp=(double *)calloc(numclut*6,sizeof(double));  // 2015-01-19

  //  delta_Term=(double *)gcemalloc(sizeof(double)*6); // 2014-07-22
  //  delta_Term=(double *)emalloc(sizeof(double)*6); // 2014-07-22 // 2014-08-13 // 2014-09-05
  delta_Term=(double *)calloc(6,sizeof(double)); // 2014-07-22 // 2014-08-13 // 2014-09-05
  //  vel_Term=(double *)gcemalloc(sizeof(double)*6); // 2014-07-22
  //  vel_Term=(double *)emalloc(sizeof(double)*6); // 2014-07-22 // 2014-09-05
  vel_Term=(double *)calloc(6,sizeof(double)); // 2014-07-22 // 2014-09-05
  //  acc_Term=(double *)gcemalloc(sizeof(double)*6); // 2014-07-22
  //  acc_Term=(double *)emalloc(sizeof(double)*6); // 2014-07-22 // 2014-09-05
  acc_Term=(double *)calloc(6,sizeof(double)); // 2014-07-22 // 2014-09-05
  //  acc_Term2=(double *)gcemalloc(sizeof(double)*6); // 2014-07-22
  //  acc_Term2=(double *)emalloc(sizeof(double)*6); // 2014-07-22 // 2014-09-05
  acc_Term2=(double *)calloc(6,sizeof(double)); // 2014-07-22 // 2014-09-05
  //  predict_Term=(double **)gcemalloc(sizeof(double *)*6); // 2014-07-22
  //  predict_Term=(double **)emalloc(sizeof(double *)*6); // 2014-07-22 // 2014-09-05
  predict_Term=(double **)calloc(6,sizeof(double *)); // 2014-07-22 // 2014-09-05
  //  predict_Term2=(double **)gcemalloc(sizeof(double *)*6); // 2014-07-22
  //  predict_Term2=(double **)emalloc(sizeof(double *)*6); // 2014-07-22 // 2014-09-05
  predict_Term2=(double **)calloc(6,sizeof(double *)); // 2014-07-22 // 2014-09-05
  //  correct_Term=(double **)gcemalloc(sizeof(double *)*6); // 2014-07-22
  //  correct_Term=(double **)emalloc(sizeof(double *)*6); // 2014-07-22 // 2014-08-13
  //  correct_Term2=(double **)gcemalloc(sizeof(double *)*6); // 2014-07-22
  //  correct_Term2=(double **)emalloc(sizeof(double *)*6); // 2014-07-22 // 2014-09-05
  correct_Term2=(double **)calloc(6,sizeof(double *)); // 2014-07-22 // 2014-09-05
  for (i=0;i<6;++i) {
    //    predict_Term[i]=(double *)gcemalloc(sizeof(double)*6); // 2014-07-22
    //    predict_Term[i]=(double *)emalloc(sizeof(double)*6); // 2014-07-22 // 2014-09-05
    predict_Term[i]=(double *)calloc(6,sizeof(double)); // 2014-07-22 // 2014-09-05
    //    predict_Term2[i]=(double *)gcemalloc(sizeof(double)*6); // 2014-07-22
    //    predict_Term2[i]=(double *)emalloc(sizeof(double)*6); // 2014-07-22 // 2014-09-05
    predict_Term2[i]=(double *)calloc(6,sizeof(double)); // 2014-07-22 // 2014-09-05
    //    correct_Term[i]=(double *)gcemalloc(sizeof(double)*6); // 2014-07-22
    //    correct_Term[i]=(double *)emalloc(sizeof(double)*6); // 2014-07-22 // 2014-08-13
    //    correct_Term2[i]=(double *)gcemalloc(sizeof(double)*6); // 2014-07-22
    //    correct_Term2[i]=(double *)emalloc(sizeof(double)*6); // 2014-07-22 // 2014-09-05
    correct_Term2[i]=(double *)calloc(6,sizeof(double)); // 2014-07-22 // 2014-09-05
  }

  vel_Term_b1=(double *)calloc(6,sizeof(double));        // 2015-01-19
  vel_Term_b2=(double *)calloc(6,sizeof(double));        // 2015-01-19
  vel_Term_b1_temp=(double *)calloc(6,sizeof(double));   // 2015-01-19

  for (i=0;i<6;++i) {   // 2015-01-19
    vel_Term_b1[i]=0.0; // 2015-01-19
    vel_Term_b2[i]=0.0; // 2015-01-19
  }                     // 2015-01-19
  for (i=0;i<numclut;++i) { // 2015-01-19
    qvel_b1[i]=0.00;        // 2015-01-19
    qvel_b2[i]=0.00;        // 2015-01-19
  }                         // 2015-01-19

  for (i=0;i<6;++i) vel_Term[i]=0.0;
  for (i=0;i<numclut;++i) qvel[i]=0.00/*1.00*/;

  DOF=(numclut-1);
  if (TERMMODE==ON) DOF+=6;
  KEobj=0.5*DOF*k_B*Tobj;


  //  KEobj=/*0.5**/DOF*k_B*Tobj;
  if (leapfrog_flag==OFF) { // 2015-01-19
    ABANH_set_new(s,s_vel,gzi,predict_gzi,correct_gzi,predict_s,correct_s,tau,&tau2,&Q_NH,KEobj,dt);
  }                         // 2015-01-19
  else {                    // 2015-01-19
    ABANH_set_new_mvV(s,s_vel,zeta,tau,&tau2,&Q_NH,KEobj,dt); // 2015-01-19
  }                         // 2015-01-19
  //  trans_A_to_CN_terminal=(double *)gcemalloc(sizeof(double)*3*3); // 2014-06-29 // 2014-07-22
  //  trans_A_to_CN_terminal=(double *)emalloc(sizeof(double)*3*3); // 2014-06-29 // 2014-07-22 // 2014-09-05
  trans_A_to_CN_terminal=(double *)calloc(3*3,sizeof(double)); // 2014-06-29 // 2014-07-22 // 2014-09-05
  //  l_Term=(double *)gcemalloc(sizeof(double)*3); // 2014-06-29 // 2014-07-22
  //  l_Term=(double *)emalloc(sizeof(double)*3); // 2014-06-29 // 2014-07-22 // 2014-09-05
  l_Term=(double *)calloc(3,sizeof(double)); // 2014-06-29 // 2014-07-22 // 2014-09-05
  for (i=0;i<3;++i) for (j=0;j<3;++j) trans_A_to_CN_terminal[i*3+j/*i][j*/]=0.0; // 2014-06-29
  for (i=0;i<3;++i) trans_A_to_CN_terminal[i*3+i/*i][i*/]=1.0; // 2014-06-29
  for (i=0;i<3;++i) l_Term[i]=0.0; // 2014-06-29

  if (leapfrog_flag==OFF) { // 2015-01-19
    ABA_integ_set(q,qvel,predict,correct,numclut,dt);
  }

  if (MODEV==ON) {
    if (leapfrog_flag==OFF) { // 2015-01-19
      ABAs_restat_read_new(inputvelofilename,numclut,correct,/*correct_Term 2014-08-13*/correct_Term_2014_08_13,correct_Term2,correct_gzi,MODE,TERMMODE);
    }                        // 2015-01-19
    else {                   // 2015-01-19
      ABAs_restat_read_new_mvV(inputvelofilename,numclut,qvel_b1,qvel_b2,vel_Term_b1,vel_Term_b2,zeta,MODE,TERMMODE);                            // 2015-01-19
    }                        // 2015-01-19
  }

  myncL_create_def_AMBER(trjfilename,numatom,&nc_id);
  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");

  //   i=2; // 2014-08-13
  //  if (i==1){ // 2014-08-13

  for (i=0;i<numstep;++i) {
    if (leapfrog_flag==OFF) { // 2015-01-19    
      ABA_integ_pret(qrot,qvel,q,predict,correct,dt,numclut);
      if (TERMMODE==ON) ABA_integ_pret_Term(predict_Term,predict_Term2,/*correct_Term 2014-08-13*/correct_Term_2014_08_13,correct_Term2,vel_Term,delta_Term,dt);

      if (MODE==NVT) ABANH_update_pret_new(&gzi,predict_gzi,correct_gzi,&s,&s_vel,predict_s,correct_s,dt);
      if (TERMMODE==ON) ABA_update_Term(crd,delta_Term,numatom,clt,numclut,trans_A_to_CN_terminal,l_Term); // 2014-06-29 // 2014-08-13 CHECK

      //    cos(qrot[0]); // 2014-08-13
      ABA_update(clt,crd,qrot,numclut,numatom);
      //    cos(qrot[0]); // 2014-08-13
      //    if (TERMMODE==ON) ABA_update_Term(crd,delta_Term,numatom,clt,numclut); // 2014-06-19

      for (j=0;j<numclut;++j) Q[j]=0.0;
      intervalflag=i%interval;
      if (AmberON==ON) {
	if (dfalg==ON) {
	  ffL_calcTorque(Q,crd,numclut,numclutparent,terminal,origin);
	}
	ffL_calcffandforce(crd,numatom,&e,&f);
	pot=0.0;
	if (dfalg==ON)
	  pot=e.p_d_t;
	for (j=0;j<numatom;++j) for (k=0;k<3;++k) frc[j*3+k]=0.0;
	for (j=0;j<numatom;++j) {
	  if (efalg   == ON) {
	    for (k=0;k<3;++k) frc[j*3+k]+=/*coff**/f.f_e[j*3+k];
	  }
	  if (LJfalg  == ON) {
	    for (k=0;k<3;++k) frc[j*3+k]+=/*coff**/f.f_LJ[j*3+k];
	  }
	  if (e14falg == ON) {
	    for (k=0;k<3;++k) frc[j*3+k]+=/*coff**/f.f_e_14[j*3+k];
	  }
	  if (LJ14falg== ON) {
	    for (k=0;k<3;++k) frc[j*3+k]+=/*coff**/f.f_LJ_14[j*3+k];
	  }
	}
	if (efalg   == ON)  pot+=0.5*e.p_e_t;
	if (LJfalg  == ON)  pot+=0.5*e.p_LJ_t;
	if (e14falg == ON)  pot+=0.5*e.p_e_14_t;
	if (LJ14falg== ON)  pot+=0.5*e.p_LJ_14_t;
      }
      else {                                                               // 2015-02-16
	calc_NB_ECEPE_for_db(&p_LJ_t,f_LJ,&p_e_t,f_e,                      // 2015-02-16
			     crd,                                          // 2015-02-16
			     pairs1_5,pairs1_4,num1_5,num1_4,              // 2015-02-16
			     ECEPE_p.NUMATM,charge,nbtype,nb_p);           // 2015-02-16
	//	calc_TORS_ECEPE2(crd,ECEPE_p,delta_dihed,&p_d_t,Q,numclut,origin); // 2015-02-16

	pot=/*p_d_t+*/p_LJ_t/*+p_e_t*/;                                            // 2015-02-16

	for (j=0;j<numatom;++j)                                          // 2015-02-16
	  for (k=0;k<3;++k)                                              // 2015-02-16
	    frc[j*3+k]=/*f_e[j*3+k]+*/f_LJ[j*3+k];                           // 2015-02-16
      }                                                                  // 2015-02-16

      ABA_calcKineE_TermOn(&KE,&KEv,&PEv,KEobj,clt,crd,qvel,s,s_vel,Q_NH,vel_Term,numclut,numatom,MODE,numclut+6);
      T=KE/(DOF*k_B)*2.0;
    //    T=KE/(DOF*k_B)/*/2.0*/;
    //    for (j=0;j<6;++j) cos(clt[0].Spacc[j]); // 2014-08-13
    //    for (j=0;j<6;++j) cos(acc_Term2[j]); // 2014-08-13
    //    for (j=0;j<6;++j) cos(acc_Term[j]); // 2014-08-13
    //    cos(qacc[1]); // 2014-08-13
    //    cos(predict[1*6+2]); // 2014-08-13
      if (TERMMODE==ON) 
	solverABA_TermOn_NH_new(qacc,qvel,clt,Q,frc,crd,numclut,numatom,q,gzi,&gzi_vel,s,&s_vel,tau2,acc_Term,acc_Term2,vel_Term,T,Tobj);
      else solverABA_NH_new(qacc,qvel,clt,Q,frc,crd,numclut,numatom,q,gzi,&gzi_vel,s,&s_vel,tau2,T,Tobj);
    //    cos(qacc[1]); // 2014-08-13
    //    cos(predict[1*6+2]); // 2014-08-13
    //    cos(qrot[0]); // 2014-08-13
    //    cos(qrot[1]); // 2014-08-13
      ABA_integ_cort(qrot,qvel,q,qacc,predict,correct,dt,numclut);
    //    cos(qrot[0]); // 2014-08-13
    //    cos(qrot[1]); // 2014-08-13
      if (TERMMODE==ON)
	ABA_integ_cort_Term(predict_Term,predict_Term2,/*correct_Term 2014-08-13*/correct_Term_2014_08_13,correct_Term2,acc_Term,acc_Term2,vel_Term,delta_Term,dt);
      if (MODE==NVT) ABANH_update_cort_new(&gzi,gzi_vel,&s,&s_vel,predict_gzi,correct_gzi,predict_s,correct_s,dt);

      if (TERMMODE==ON) ABA_update_Term(crd,delta_Term,numatom,clt,numclut,trans_A_to_CN_terminal,l_Term); // 2014-06-29 // 2014-08-13 CHECK
    //    cos(qrot[0]); // 2014-08-13
      ABA_update(clt,crd,qrot,numclut,numatom); // 2014-08-13 CHECK

    //    if (TERMMODE==ON) ABA_update_Term(crd,delta_Term,numatom,clt,numclut); // 2014-06-19
    //    cos(qrot[0]); // 2014-08-13
      if (TERMMODE==ON)
	ABA_calcKineE_TermOn(&KE,&KEv,&PEv,KEobj,clt,crd,qvel,s,s_vel,Q_NH,vel_Term,numclut,numatom,MODE,numclut+6);
      else
	ABA_calcKineE(&KE,&KEv,&PEv,KEobj,clt,crd,qvel,s,s_vel,Q_NH,numclut,numatom,MODE);
      if (MODE==NVT)
	ABANH_calcKE_new(gzi,s,s_vel,Q_NH,KEobj,&PEv,&KEv);
    }
    else { // 2015-01-19
      for (j=0;j<numclut;++j) qvel[j]=1.5*qvel_b1[j]-0.5*qvel_b2[j];
      if (TERMMODE==ON) for (j=0;j<6;++j) vel_Term[j]=1.5*vel_Term_b1[j]-0.5*vel_Term_b2[j];

      for (j=0;j<numclut;++j) Q[j]=0.0;
      intervalflag=i%interval;

      if (AmberON==ON) {
	if (dfalg==ON) {
	  ffL_calcTorque(Q,crd,numclut,numclutparent,terminal,origin);
	}
	ffL_calcffandforce(crd,numatom,&e,&f);
	pot=0.0;
	if (dfalg==ON)
	  pot=e.p_d_t;
	for (j=0;j<numatom;++j) for (k=0;k<3;++k) frc[j*3+k]=0.0;
	for (j=0;j<numatom;++j) {
	  if (efalg   == ON) {
	    for (k=0;k<3;++k) frc[j*3+k]+=/*coff**/f.f_e[j*3+k];
	  }
	  if (LJfalg  == ON) {
	    for (k=0;k<3;++k) frc[j*3+k]+=/*coff**/f.f_LJ[j*3+k];
	  }
	  if (e14falg == ON) {
	    for (k=0;k<3;++k) frc[j*3+k]+=/*coff**/f.f_e_14[j*3+k];
	  }
	  if (LJ14falg== ON) {
	    for (k=0;k<3;++k) frc[j*3+k]+=/*coff**/f.f_LJ_14[j*3+k];
	  }
	}
	if (efalg   == ON)  pot+=0.5*e.p_e_t;
	if (LJfalg  == ON)  pot+=0.5*e.p_LJ_t;
	if (e14falg == ON)  pot+=0.5*e.p_e_14_t;
	if (LJ14falg== ON)  pot+=0.5*e.p_LJ_14_t;
      }
      else {                                                             // 2015-02-16

	calc_NB_ECEPE_for_db(&p_LJ_t,f_LJ,&p_e_t,f_e,                      // 2015-02-16
			     crd,                                          // 2015-02-16
			     pairs1_5,pairs1_4,num1_5,num1_4,              // 2015-02-16
			     ECEPE_p.NUMATM,charge,nbtype,nb_p);           // 2015-02-16
	calc_TORS_ECEPE2(crd,ECEPE_p,delta_dihed,&p_d_t,Q,numclut,origin); // 2015-02-16

	pot=p_d_t+p_LJ_t+p_e_t;                                          // 2015-02-16

	for (j=0;j<numatom;++j)                                          // 2015-02-16
	  for (k=0;k<3;++k)                                              // 2015-02-16
	    frc[j*3+k]=f_e[j*3+k]+f_LJ[j*3+k];                           // 2015-02-16

      }                                                                  // 2015-02-16


      if (TERMMODE==ON) 
	solverABA_TermOn_NH_new_mvV(qacc,qvel,clt,Q,frc,crd,numclut,numatom,q,zeta,acc_Term,acc_Term2,vel_Term);
      else solverABA_NH_new_mvV(qacc,qvel,clt,Q,frc,crd,numclut,numatom,q,zeta);

      for (j=0;j<numclut;++j) {
	qvel_b1_temp[j]=qvel_b1[j]+dt*qacc[j];
	qvel[j]=0.5*qvel_b1[j]+0.5*qvel_b1_temp[j];
      }

      if (TERMMODE==ON) {
	for (j=0;j<6;++j) {
	  vel_Term_b1_temp[j]=vel_Term_b1[j]+dt*acc_Term[j];
	  vel_Term[j]=0.5*vel_Term_b1[j]+0.5*vel_Term_b1_temp[j];
	}
      }

      for (j=0;j<numclut;++j) {
	//	qrot[j]=/*dt*qvel_b1_temp[j]*/dt*qvel[j]+0.5*dt*dt*qacc[j];
	qrot[j]=dt*qvel_b1_temp[j]/*dt*qvel[j]+0.5*dt*dt*qacc[j]*/;
	qvel_b2[j]=qvel_b1[j];
	qvel_b1[j]=qvel_b1_temp[j];
      }
      if (TERMMODE==ON) {
	for (j=0;j<6;++j) {
	  //	  delta_Term[j]=/*dt*vel_Term_b1_temp[j]*/dt*vel_Term[j]+0.5*dt*dt*acc_Term[j];
	  delta_Term[j]=dt*vel_Term_b1_temp[j]/*dt*vel_Term[j]+0.5*dt*dt*acc_Term[j]*/;
	  vel_Term_b2[j]=vel_Term_b1[j];
	  vel_Term_b1[j]=vel_Term_b1_temp[j];
	}
      }
      
      if (MODE==NVT) { 
	ABA_calcKineE_TermOn(&KE_1_2,&KEv_1_2,&PEv_1_2,KEobj,clt,crd,qvel_b1,s,s_vel,Q_NH,vel_Term_b1,numclut,numatom,MODE,numclut+6);
	T_1_2=KE_1_2/(DOF*k_B)*2.0;
	zeta=zeta+dt*1.0/tau2*(T_1_2/Tobj-1.0);
	//////////////////////////////////////
	s_vel=zeta*s;
	s=s+dt*s_vel;
	//////////////////////////////////////
      }

      ABA_update(clt,crd,qrot,numclut,numatom);
      if (TERMMODE==ON) ABA_update_Term(crd,delta_Term,numatom,clt,numclut,trans_A_to_CN_terminal,l_Term);
    
      if (TERMMODE==ON)
	ABA_calcKineE_TermOn(&KE,&KEv,&PEv,KEobj,clt,crd,qvel,s,s_vel,Q_NH,vel_Term,numclut,numatom,MODE,numclut+6);
      else
	ABA_calcKineE(&KE,&KEv,&PEv,KEobj,clt,crd,qvel,s,s_vel,Q_NH,numclut,numatom,MODE);
      if (MODE==NVT)
	ABANH_calcKE_new(zeta,s,s_vel,Q_NH,KEobj,&PEv,&KEv);
      T=KE/(DOF*k_B)*2.0;
    } // 2015-01-19

    if (i%interval==0)                                                            // 2014-08-13
      fprintf(outputfile,"%d %24.20e %24.20e %24.20e %24.20e %24.20e %24.20e %24.20e\n", // 2014-08-13
	      i,pot,KE,KEv,PEv,pot+KE+PEv+KEv,s,T);                                      // 2014-08-13
    //    if (i%interval==0) {                                                        // 2014-08-13
    //      fprintf(outputfile,"%d %e %e %e %e\n",                                         // 2014-08-13
    //    	      i,s,PEv,KEv,PEv+KEv/*,pot,KEv,PEv,pot+KE+PEv+KEv,s,T*/);                      // 2014-08-13
    //    }                                                                           // 2014-08-13

      /*****************************************************************/
      /* pot_KE_KEv_PEv=pot+KE+PEv+KEv;				       */
      /* fprintf(outputfile,"%d ",i); // 2014-08-13		       */
      /* fprintf(outputfile,"%24.20lf ",pot); // 2014-08-13	       */
      /* fprintf(outputfile,"%24.20lf ",KE); // 2014-08-13	       */
      /* fprintf(outputfile,"%24.20lf ",KEv); // 2014-08-13	       */
      /* fprintf(outputfile,"%24.20lf ",PEv); // 2014-08-13	       */
      /* fprintf(outputfile,"%24.20lf ",pot_KE_KEv_PEv); // 2014-08-13 */
      /* fprintf(outputfile,"%24.20lf ",s); // 2014-08-13	       */
      /* fprintf(outputfile,"%24.20lf \n",T); // 2014-08-13	       */
      /*****************************************************************/

    // 2014-08-13
    if (i%intervalout==0) { 
      //      ffL_out_formated(outputfile2,e,KE,KEv,PEv,T,i,dt);
      fprintf(outputfile2,"/***********************************************/\n");
      fprintf(outputfile2,"steps            = %d  \n",i);
      fprintf(outputfile2,"total time       = %10.3lf ps  \n" ,dt*(double)i);
      fprintf(outputfile2,"T_kelvin         = %24.20e K  \n",T);
      fprintf(outputfile2,"toal_energy      = %24.20e kcal/mol  \n",e.p_t+KE);
      fprintf(outputfile2,"toal_vertial_energy      = %24.20e kcal/mol  \n",e.p_t+KE+KEv+PEv);
      fprintf(outputfile2,"kinetic_energy   = %24.20e kcal/mol  \n",KE);
      fprintf(outputfile2,"kinetic_energy_vertial   = %24.20e kcal/mol  \n",KEv);
      fprintf(outputfile2,"potential_energy_real = %24.20e kcal/mol  \n",e.p_t);
      fprintf(outputfile2,"potential_energy_vertial   = %24.20e kcal/mol  \n",PEv);
      fprintf(outputfile2,"dihedral_energy  = %24.20e kcal/mol  \n",e.p_d_t);
      fprintf(outputfile2,"elect_energy     = %24.20e kcal/mol  \n",e.p_e_t);
      fprintf(outputfile2,"VDW_energy       = %24.20e kcal/mol  \n",e.p_LJ_t);
      fprintf(outputfile2,"1_4_elect_energy = %24.20e kcal/mol  \n",e.p_e_14_t);
      fprintf(outputfile2,"1_4_VDW_energy   = %24.20e kcal/mol  \n",e.p_LJ_14_t);
    }
    // 2014-08-13

    // 2014-08-13
    if (i%intervalnc==0) {
      if (TERMMODE==ON) {                                                       // 2014-06-30
	ABA_update_Term2(crd,crd_Term,numatom,trans_A_to_CN_terminal,l_Term);   // 2014-06-30 CHECK
	for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crd_Term[j*3+k]; // 2014-06-30 CHECK
      }                                                                         // 2014-06-30
      else {                                                                    // 2014-06-30
	for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crd[j*3+k]; // CHECK
      }                                                                         // 2014-06-30 CHECK
      myncL_put_crd_AMBER(nc_id,l,crd_nc);
      ++l;
    }
    // 2014-08-13

    //    i=2; // 2014-08-13
    //    if (i==1){ // 2014-08-13
    //      ; // 2014-08-13
    //    } // 2014-08-13
  }
  fclose(outputfile);
  fclose(outputfile2);

  nc_close((nc_id.ncid));

  // 2014-08-13
  rstfile=efopen(rstfilename,"w");
  if (AmberON==ON) {             // 2015-02-16
    fprintf(rstfile,"ACE\n ");
  }                               // 2015-02-16
  fprintf(rstfile,"%d \n",numatom);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      fprintf(rstfile,"%10.8lf ",crd[i*3+j]);
    }
    fprintf(rstfile,"\n");
  }
  fclose(rstfile);
  // 2014-08-13

  if (leapfrog_flag==OFF) { // 2015-01-19
    ABAs_restat_write_vel_new(rstvelfilename,numclut,correct,correct_Term_2014_08_13,correct_Term2,correct_gzi,MODE,TERMMODE); // 2014-08-13 
  }
  else { // 2015-01-19
    ABAs_restat_write_vel_new_mvV(rstvelfilename,numclut,qvel_b1,qvel_b2,vel_Term_b1,vel_Term_b2,zeta,MODE,TERMMODE);        // 2015-01-19
  }      // 2015-01-19

  // 2014-07-22
  free(AP.AMASS);
  free(AP.IAC);
  free(AP.NUMEX);
  free(AP.ICO);
  free(AP.IPRES);
  free(AP.RK);
  free(AP.REQ);
  free(AP.TK);
  free(AP.TEQ);
  free(AP.PK);
  free(AP.PN);
  free(AP.PHASE);
  free(AP.SOLTY);
  free(AP.CN1);
  free(AP.CN2);
  free(AP.BH);
  free(AP.BA);
  free (AP.TH);
  free(AP.TA);
  free(AP.PH);
  free(AP.PA);
  free(AP.NATEX);
  free(AP.ASOL);
  free(AP.BSOL);
  free(AP.HBCUT);
  free(AP.JOIN);
  free(AP.IROTAT);
  free(AP.NSP);
  free(AP.BPER);
  free(AP.ICBPER);
  free(AP.TPER);
  free(AP.ICTPER);
  free(AP.PPER);
  free(AP.ICTPER);
  free(AP.IAPER);
  free(AP.IACPER);
  free(AP.CGPER);
  free(AP.LES_TYPE);
  free(AP.LES_FAC);
  free(AP.LES_CNUM);
  free(AP.LES_ID);

  free(mass);
  free(crd);
  free(crd_Term);
  free(clt);
  free(IndexOfABICycle);
  free(numclutparent);
  free(terminal);
  free(origin);

  free(Q); 
  free(frc);

  free(q);
  free(qacc);
  free(qvel);
  free(qrot);
  free(predict);
  free(correct);

  //  free(delta_Term); // 2014-08-13
  free(vel_Term);
  free(acc_Term);
  free(acc_Term2);
  free(predict_Term);
  free(predict_Term2);
  free(/*correct_Term 2014-08-13*/correct_Term_2014_08_13); // 2014-08-13
  free(correct_Term2);
  free(trans_A_to_CN_terminal);
  free(l_Term);

  // 2014-07-22

  free(qvel_b1);       // 2015-01-19
  free(qvel_b2);       // 2015-01-19
  free(qvel_b1_temp);  // 2015-01-19

  free(vel_Term_b1);        // 2015-01-19
  free(vel_Term_b2);        // 2015-01-19
  free(vel_Term_b1_temp);   // 2015-01-19

  free(charge);             // 2015-02-16
  free(nbtype);             // 2015-02-16

  free(pairs1_5);    // 2015-02-16
  free(pairs1_4);    // 2015-02-16
  free(num1_4);      // 2015-02-16
  free(num1_5);      // 2015-02-16
  free(delta_dihed); // 2015-02-16

  return 0;
}

int USAGE(char *progname) {                                                             
  printf("USAGE:\n");                                                                   
  printf("[-f]\n");
  printf("[--nve]\n");
  printf("[--termon]\n");
  printf("[--dih]\n");
  printf("[--els]\n");
  printf("[--lj]\n");
  printf("[--e14]\n");
  printf("[--l14]\n");
  printf("[-h]\n");
  printf("[--kc]\n");
  printf("[--cof]\n");
  printf("[--temp]\n");
  printf("[--dt]\n");
  printf("[--omg]\n");
  printf("[--nums]\n");
  printf("[--int]\n");
  printf("[--intnc]\n");
  printf("[--intout]\n");
  printf("[--vel]\n");
  printf("[--tau]\n");
  printf("[--rst]\n");
  printf("[--rstv]\n");
  printf("[-h] help \n");                                                               
  printf("%s [-h] inputfilename clustfilename parmfilename outputfilename\n",progname); 
  return 0; // 2014-07-04
}                                                                                       
// 2014-06-23
