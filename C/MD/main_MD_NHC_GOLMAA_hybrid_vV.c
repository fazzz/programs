#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "MD.h"
#include "NHC.h"
#include "GOLMAA_hybrid_set.h"
#include "GOLMAA_hybrid.h"

#include "PTL.h"
#include "FFL.h"
#include "EF.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,m,d;
  int numatom,numres,numstep=100000;
  int interval=1000,intervalout=1000,intervalnc=1000;
  double dt=0.001;
  double *frc,pot;
  double *a,*v;
  double KE,KEv,PEv;

  double Tobj=300,KEobj,KBT;
  double k_B=1.98723e-3;

  double *zeta,*zeta_vel,*zeta_acc,**predict_zeta,**correct_zeta;
  double *Q_NH,tau=0.1,tau2;
  double T;
  int DOF,M=4;

  int bflag=ON,aflag=ON,dflag=ON,eflag=ON,LJflag=ON,e14flag=ON,LJ14flag=ON,natflag=ON,nnatflag=ON;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_AMBER nc_id;

  double *crd,*refcrd,*mass;
  struct potential e;
  struct force f;
  struct potential_GOLMAA_hybrid e_GOLM;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*refcrdfilename,*inputvelofilename;
  char *outputfilename,*outputfilename2,*parmfilename,*trjfilename;
  char *rstfilename="rstfile",*rstvelfilename="rstvelfile";

  FILE *inputfile,*reffile,*inputvelofile;
  FILE *outputfile,*outputfile2,*parmfile;
  FILE *rstfile,*rstvelfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"f",0,NULL,'f'},
    {"nve",0,NULL,'*'},
    {"bnd",0,NULL,'b'},
    {"ang",0,NULL,'a'},
    {"dih",0,NULL,'d'},
    {"e14",0,NULL,'1'},
    {"l14",0,NULL,'4'},
    {"nat",0,NULL,'c'},
    {"repul",0,NULL,'N'},
    {"h",0,NULL,'h'},
    {"temp",1,NULL,'t'},
    {"dt",1,NULL,'@'},
    {"nums",1,NULL,'s'},
    {"int",1,NULL,'i'},
    {"intnc",1,NULL,'j'},
    {"intout",1,NULL,'k'},
    {"vel",1,NULL,'v'},
    {"tau",1,NULL,'?'},
    {"nchain",1,NULL,'M'},
    {"rst",1,NULL,'{'},
    {"rstv",1,NULL,'}'},
    {0,0,0,0}
  };

  MODE=NVT;
  MODEV=OFF;

  while((c=getopt_long(argc,argv,"*+bad14cNh@:t:k:c:t:o:s:i:j:k:v:?:M:{:}:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case '*':
      MODE=NVE;
      break;
    case 'b':
      bflag=OFF;
      break;
    case 'a':
      aflag=OFF;
      break;
    case 'd':
      dflag=OFF;
      break;
    case '1':
      e14flag=OFF;
      break;
    case '4':
      LJ14flag=OFF;
      break;
    case 'c':
      natflag=OFF;
      break;
    case 'N':
      nnatflag=OFF;
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
    case 'v':
      inputvelofilename=optarg;
      MODEV=ON;
      break;
    case '?':
      tau=atof(optarg);
      break;
    case 'M':
      M=atoi(optarg);
      break;
    case '{':
      rstfilename=optarg;
      break;
    case '}':
      rstvelfilename=optarg;
      break;
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

  if (argc < 6) {
    USAGE(progname);
    exit(1);
  }
  inputfilename  = *argv;
  refcrdfilename = *++argv;
  parmfilename   = *++argv;
  outputfilename = *++argv;
  outputfilename2= *++argv;
  trjfilename    = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  numres=AP.NRES;
  mass=(double *)gcemalloc(sizeof(double)*numatom);
  for (i=0;i<numatom;++i) mass[i]=AP.AMASS[i];

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      fscanf(inputfile,"%lf",&crd[i*3+j]);
    }
  }
  fclose(inputfile);

  refcrd=(double *)gcemalloc(sizeof(double)*numatom*3);
  reffile=efopen(refcrdfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(reffile,"%d",&d);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      fscanf(reffile,"%lf",&refcrd[i*3+j]);
    }
  }
  fclose(reffile);

  a=(double *)gcemalloc(sizeof(double)*numatom*3);
  v=(double *)gcemalloc(sizeof(double)*numatom*3);

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  ffL_set_calcffandforce(&e,&f);

  GOLMAA_hybrid_ff_set_calcff(&e_GOLM,refcrd,numatom,numres,e.parm.indexnb,e.parm.numnb);

  zeta=(double *)gcemalloc(sizeof(double)*M);
  zeta_vel=(double *)gcemalloc(sizeof(double)*M);
  zeta_acc=(double *)gcemalloc(sizeof(double)*M);
  predict_zeta=(double **)gcemalloc(sizeof(double *)*M);
  correct_zeta=(double **)gcemalloc(sizeof(double *)*M);
  for (i=0;i<M;++i) {
    predict_zeta[i]=(double *)gcemalloc(sizeof(double)*6);
    correct_zeta[i]=(double *)gcemalloc(sizeof(double)*6);
  }
  Q_NH=(double *)gcemalloc(sizeof(double)*M);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j)
      v[i*3+j]=0.000/*1.00*/;

  DOF=(numatom*3);
  KEobj=0.5*DOF*k_B*Tobj;
  KBT=k_B*Tobj;

  if (MODEV==ON) {
    inputvelofile=efopen(inputvelofilename,"r");
    fclose(inputvelofile);
  }
  
  myncL_create_def_AMBER(trjfilename,numatom,&nc_id);
  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");

  for (i=0;i<numstep;++i) {

    KE=MD_Propagetor_vV_NVE(crd,v,mass,numatom,dt,e,f);

    if (MODE==NVT) 
      NHC_update_pret(zeta,zeta_vel,predict_zeta,correct_zeta,M,dt);


    //    ffL_calcffandforce_14D_woH(crd,numatom,&e,&f);
    /*******************************************************************************/
    /* ffL_calcffandforce(crd,numatom,&e,&f);					   */
    /* for (j=0;j<numatom;++j) {						   */
    /*   //      if (strncmp(AP.IGRAPH[j],"H",1)!=0) {				   */
    /*   if(bflag==ON) for (k=0;k<3;++k) frc[j*3+k]-=f.f_b[j*3+k];		   */
    /*   if(aflag==ON) for (k=0;k<3;++k) frc[j*3+k]+=f.f_a[j*3+k];		   */
    /*   if(dflag==ON) for (k=0;k<3;++k) frc[j*3+k]+=f.f_d[j*3+k];		   */
    /*   if (e14flag == ON) for (k=0;k<3;++k) frc[j*3+k]+=f.f_e_14[j*3+k];	   */
    /*   if (LJ14flag== ON) for (k=0;k<3;++k) frc[j*3+k]+=f.f_LJ_14[j*3+k];	   */
    /* 	//      }								   */
    /* }									   */
    /*******************************************************************************/

    /********************************************************************/
    /* if (natflag == ON || nnatflag == ON ) 			        */
    /*   GOLMAA_hyb_ff_calcff(crd,numatom,&e_GOLM);		        */
    /* for (j=0;j<numatom;++j) {				        */
    /*    if (natflag == ON)					        */
    /* 	 for (k=0;k<3;++k) frc[j*3+k]+=e_GOLM.f_natatt[j][k];	        */
    /*    if (nnatflag == ON)					        */
    /* 	 for (k=0;k<3;++k) frc[j*3+k]+=e_GOLM.f_repul[j][k];	        */
    /* }							        */
    /********************************************************************/

    KE=0.0;
    for (j=0;j<numatom;++j) 
      for (k=0;k<3;++k) 
	KE+=0.5*mass[j]*v[j*3+k]*v[j*3+k]/(4.18407*100.0);

    T=KE/(DOF*k_B)*2.0;

    for (j=0;j<numatom;++j)
      for (k=0;k<3;++k)
	  a[j*3+k]=frc[j*3+k]/mass[j];
    if (MODE==NVT) {
      for (j=0;j<numatom;++j)
	for (k=0;k<3;++k)
	  a[j*3+k]-=v[j*3+k]*zeta_vel[0];
      NHC_solve(zeta_vel,zeta_acc,Q_NH,M,1,KBT,tau2,T,Tobj);
    }

    for (j=0;j<numatom;++j)
      for (k=0;k<3;++k)
	for (l=0;l<6;++l) 
	  correct[j][k][l] = predict[j][k][l]+GearsConstant[l]*(0.5*dt*dt*a[j*3+k]-predict[j][k][2]);

    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k) {
	crd[j*3+k]=correct[j][k][0];
	v[j*3+k]=correct[j][k][1]/dt;
      }
    }

    if (MODE==NVT) NHC_update_cort(zeta,zeta_vel,zeta_acc,predict_zeta,correct_zeta,M,dt);

    KE=0.0;
    for (j=0;j<numatom;++j) 
      for (k=0;k<3;++k) 
	KE+=0.5*mass[j]*v[j*3+k]*v[j*3+k]/(4.18407*100.0);

    T=KE/(DOF*k_B)*2.0;

    if (MODE==NVT)
      NHC_calcKE(zeta,zeta_vel,Q_NH,M,DOF,KBT,&PEv,&KEv);

    if (i%interval==0) {
      pot=0.0;
      if (bflag==ON) {
	//      for (j=AP.NBONH;j<AP.MBONA;++j) {    
	for (j=0;j<AP.NBONH+AP.MBONA;++j) {    
	  //	pot+=e.p_b[j];
	}
      }
      if (aflag==ON) {
	//      for (j=AP.NTHETH;j<AP.MTHETA;++j) {
	for (j=0;j<AP.NTHETH+AP.MTHETA;++j) {
	  pot+=e.p_a[j];
	}
      }
      if (dflag==ON) {
	//	for (j=AP.NPHIH;j<AP.MPHIA;++j) {
	for (j=0;j<AP.NPHIH+AP.MPHIA;++j) {
	  pot+=e.p_d[j];
	}
      }
      for (j=0;j<numatom;++j) {
	//	if (strncmp(AP.IGRAPH[j],"H",1)!=0) {
	  if (e14flag == ON)  pot+=0.5*e.p_e_14[j];
	  if (LJ14flag== ON)  pot+=0.5*e.p_LJ_14[j];
	  //	}
      }
      //      if (natflag == ON)  pot+=0.5*e_GOLM.p_natatt_t;
      //      if (nnatflag == ON) pot+=0.5*e_GOLM.p_repul_t;
      fprintf(outputfile,"%d %e %e %e %e %e %e %e\n",i,pot,KE,KEv,PEv,pot+KE+PEv+KEv,zeta[0],T);
    }

    if (i%intervalout==0) {
      pot=0.0;
      if (dflag==ON) pot=e.p_d_t;
      if (e14flag == ON)  pot+=0.5*e.p_e_14_t;
      if (LJ14flag== ON)  pot+=0.5*e.p_LJ_14_t;
      if (natflag == ON)  pot+=e_GOLM.p_natatt_t;
      if (nnatflag == ON) pot+=e_GOLM.p_repul_t;
      fprintf(outputfile2,"E_t    = %e \n",pot+KE+PEv+KEv);
      fprintf(outputfile2,"KE     = %e \n",KE);
      fprintf(outputfile2,"KEv    = %e \n",KEv);
      fprintf(outputfile2,"PEv    = %e \n",PEv);
      fprintf(outputfile2,"p_t    = %e \n",pot);
      fprintf(outputfile2,"p_NC   = %e \n",e_GOLM.p_natatt_t);
      fprintf(outputfile2,"p_NNC  = %e \n",e_GOLM.p_repul_t);
      fprintf(outputfile2,"p_14es = %e \n",e.p_e_14_t);
      fprintf(outputfile2,"p_14LJ = %e \n",e.p_LJ_14_t);
      fprintf(outputfile2,"p_dihe = %e \n",e.p_d_t);
    }

    if (i%intervalnc==0) {
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crd[j*3+k];
      myncL_put_crd_AMBER(nc_id,l,crd_nc);
      ++l;
    }
  }
  fclose(outputfile);
  fclose(outputfile2);
  nc_close((nc_id.ncid));

  rstfile=efopen(rstfilename,"w");
  fprintf(rstfile,"ACE\n ");
  fprintf(rstfile,"%d \n",numatom);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      fprintf(rstfile,"%10.8lf ",crd[i*3+j]);
    }
    fprintf(rstfile,"\n");
  }
  fclose(rstfile);

  if (MODEV==ON) {
    rstvelfile=efopen(rstvelfilename,"w");
    l=0;
    for (i=0;i<numatom;++i) {
      for (j=0;j<3;++j) {
	for (k=0;k<6;++k) {
	  ++l;
	  fprintf(inputvelofile,"%e ",correct[i][j][k]);
	  if (l%10==0)
	    fscanf(inputvelofile,"\n");
	}
      }
    }
    if (MODE==NVT)  {
      l=0;
      for (i=0;i<M;++i) {
	for (j=0;j<6;++j) {
	  ++l;
	  fprintf(inputvelofile,"%e",correct_zeta[i][j]);
	  if (l%10==0)
	    fscanf(inputvelofile,"\n");
	}
      }
    }
    fclose(inputvelofile);
  }

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename refcrdfilename clustfilename parmfilename outputfilename\n",progname);
}


