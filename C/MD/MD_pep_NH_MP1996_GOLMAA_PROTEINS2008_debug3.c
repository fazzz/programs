#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"

#include "PTL.h"
#include "EF.h"
#include "RAND.h"
#include "BOXMULL.h"
#include "MD.h"
#include "MD_NHC_MP1996.h"

// debug
#include "TOPO.h"
#include "LA.h"
#include "mymath.h"
// debug

#include "netcdf_mineL.h"

#define NVT 1
#define NVE 0

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,d;
  int numatom,numheavyatom,numres,numstep=10000,interval=100;
  double dt=0.001,dt2,wdt2[3],wdt4[3];
  //  double *frc,PE;
  double pi;

  double ep=ep_natatt_hybrid;

  int vMode=OFF,MODE=NVT,NCmode=3,nibnum=3,criteria=6.5;
  // debug
  int ii,jj,kk,ll;
  int Mode14=OFF;
  int Moded=OFF,bflag=ON,aflag=ON,dflag=ON,nflag=ON;
  int Moded2=OFF;
  double m[3],n[3],m_n[3],n_n[3],lm,ln;
  double vij[3],vkj[3],vkl[3];
  double lkj;
  double vijvkj,vklvkj;

  double atom[4][3];
  double dihed;
  // debug

  int nc=1;
  double T0=300,T,K0,KE;
  double k_B=1.98723e-3;
  double UNITT=418.4070;
  double NfKT,KT;
  double zeta=0.0,V_zeta=0.0,Q,tau=0.01,tau2;
  double PEv,KEv;

  double avePE=0.0,varPE=0.0,aveKE=0.0,varKE=0.0,aveT=0.0,varT=0.0;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;

  double *crd,*refcrd,*refcrdAA,*mass,*vel;

  double summass,COM[3];

  struct potential e;
  struct force f;
  struct potential_GOLMAA_PROTEINS2008 e_GOLM;
  double p_t=0.0,E_t;

  double x[3];

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*refcrdfilename,*velfilename,*parmfilename;
  char *trjfilename,*outputfilename,*outputfilename2,*rstfilename="rstcrd",*rstvelfilename="rstvel";

  char *logfilename="MD_pep_NH_MP1996_GOLMAA_PROTEINS2008.log";

  FILE *inputfile,*refcrdfile,*velfile,*parmfile;
  FILE *outputfile,*outputfile2,*rstfile,*rstvelfile;

  FILE *logfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"nve",0,NULL,'*'},
    {"vMode",1,NULL,'v'},
    {"ep",1,NULL,'e'},
    {"nums",1,NULL,'s'},
    {"temp",1,NULL,'t'},
    {"tau",1,NULL,'a'},
    {"int",1,NULL,'i'},
    {"rst",1,NULL,'{'},
    {"rstvel",1,NULL,'}'},
    {"dt",1,NULL,'x'},
    {"cutoff",1,NULL,'c'},
    {"14",0,NULL,'1'}, // debug
    {"bond",0,NULL,'2'}, // debug
    {"angl",0,NULL,'3'}, // debug
    {"dihe",0,NULL,'4'}, // debug
    {"nonb",0,NULL,'0'}, // debug
    {"debg",0,NULL,'d'}, // debug
    {"debg2",0,NULL,'k'}, // debug
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"*h12340dks:v:e:t:a:i:x:N:b:c:{:}:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case '*':
      MODE=NVE;
      break;
    case 's':
      numstep=atoi(optarg);
      break;
    case 'v':
      vMode=ON;
      velfilename=optarg;
      break;
    case 'e':
      ep=atof(optarg);
      break;
    case 't':
      T0=atof(optarg);
      break;
    case 'a':
      tau=atof(optarg);
      break;
    case 'i':
      interval=atoi(optarg);
      break;
    case '{':
      rstfilename=optarg;
      break;
    case '}':
      rstvelfilename=optarg;
      break;
    case 'x':
      dt=atof(optarg);
      break;
    case 'c':
      criteria=atof(optarg);
      break;
    case '1':
      Mode14=ON;
      break;
    case '2':
      bflag=OFF;
      break;
    case '3':
      aflag=OFF;
      break;
    case '4':
      dflag=OFF;
      break;
    case '0':
      nflag=OFF;
      break;
    case 'd':
      Moded=ON;
      break;
    case 'k':
      Moded2=ON;
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
  inputfilename     = *argv;
  refcrdfilename    = *++argv;
  parmfilename      = *++argv;
  outputfilename    = *++argv;
  outputfilename2   = *++argv;
  trjfilename       = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  j=0;
  for (i=0;i<numatom;++i) {
    if (strncmp(AP.IGRAPH[i],"H",1)==0) {
      ++j;
    }
  }
  numheavyatom=numatom-j;
  numres=AP.NRES;
  mass=(double *)gcemalloc(sizeof(double)*numatom);
  for (i=0;i<numatom;++i) mass[i]=AP.AMASS[i];
  
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrd=(double *)gcemalloc(sizeof(double)*numatom*3);
  vel=(double *)gcemalloc(sizeof(double)*numatom*3);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);

  refcrdfile=efopen(refcrdfilename,"r");
  getline(&line,&len,refcrdfile);
  fscanf(refcrdfile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(refcrdfile,"%lf",&refcrd[i*3+j]);
  fclose(refcrdfile);

  if ( vMode==OFF ) {
    MD_Generate_inivelo(vel,mass,numatom,k_B*T0*UNITT);
    for (i=0;i<numatom;++i) {
      if (strncmp(AP.IGRAPH[i],"H",1)==0) {
	for (j=0;j<3;++j) {
	  vel[i*3+j]=0.0;
	}
      }
    }
    zeta=0.0;
    V_zeta=0.0;
  }
  else {
    velfile=efopen(velfilename,"r");
    for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(velfile,"%lf",&vel[i*3+j]);
    fclose(velfile);
  }
  K0=0.0;
  for (i=0;i<numatom;++i) 
    for (j=0;j<3;++j) 
      K0+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];
  T=K0/((3*numheavyatom)*k_B)*2.0/UNITT;

  pi=acos(-1.0);
  tau=tau/2.0/pi;         
  tau2=tau*tau;          
  KT=k_B*T0;
  NfKT=(3.0*numheavyatom+1)*KT*UNITT;
  Q=tau2*KT*UNITT*(3.0*numheavyatom);

  summass=0.0;
  for (i=0;i<numatom;++i) summass+=mass[i];
  for (i=0;i<3;++i) COM[i]=0.0;
  for (i=0;i<numatom;++i) 
    for (j=0;j<3;++j) 
      COM[j]+=mass[i]*crd[i*3+j]/summass;
  for (i=0;i<numatom;++i) 
    for (j=0;j<3;++j) 
      crd[i*3+j]-=COM[j];

  //  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  ffL_set_calcffandforce(&e,&f);
  // debug
  if (Mode14==OFF) {
  // debug
    GOLMAA_PROTEINS2008_ff_set_calcff(&e_GOLM,refcrd,numatom,numres,e.parm.indexnb,e.parm.numnb,ep,nibnum,criteria);

    if (Moded2==ON) 
      GOLMAA_PROTEINS2008_ff_calcff_debug(crd,numatom,&e_GOLM,0);
    else
      GOLMAA_PROTEINS2008_ff_calcff(crd,numatom,&e_GOLM);
  // debug
  }
  else
    ffL_calcffandforce_14DAB_woH(crd,numatom,&e,&f);
  // debug

  MD_Propagetor_NH_Single_set_MP1996(nc,dt,&dt2,wdt2,wdt4);
  myncL_create_def_MCD(trjfilename,numatom,&nc_id_MCD);
  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");
  for (i=0;i<numstep;++i) {
    if (MODE==NVT) {
      // debug
      if (Mode14==OFF) {
	// debug
	if (Moded==ON) {
	  KE=MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008_debug(crd,vel,mass,&zeta,&V_zeta,Q,NfKT,numatom,&KEv,&PEv,dt,dt2,nc,wdt4,wdt2,&e_GOLM,bflag,aflag,dflag,nflag);
	}
	else if (Moded2==ON) {
	  KE=MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008_debug2(crd,vel,mass,&zeta,&V_zeta,Q,NfKT,numatom,&KEv,&PEv,dt,dt2,nc,wdt4,wdt2,&e_GOLM,i);
	}
	else {
	  KE=MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008(crd,vel,mass,&zeta,&V_zeta,Q,NfKT,numatom,&KEv,&PEv,dt,dt2,nc,wdt4,wdt2,&e_GOLM);
	}
      }
      // debug
      else {
      KE=MD_Propagetor_NH_MP1998_14LJdab(crd,vel,mass,&zeta,&V_zeta,Q,NfKT,numatom,&KEv,&PEv,dt,dt2,nc,wdt4,wdt2,&e,&f);
      }
    // debug
    }
    else {
      // debug
      if (Mode14==OFF) {
      // debug
      KE=MD_Propagetor_vV_NVE_GOLMAA_PROTEINS2008(crd,vel,mass,numatom,dt,&e_GOLM);
      // debug
      }
      else {
	KE=MD_Propagetor_vV_NVE_14LJdab(crd,vel,mass,numatom,dt,&e,&f);
      }
      // debug
    }

    // debug
    if (i>/*91001*/91084) {
      for (j=0;j<(e_GOLM).num_bond;++j) {
	printf("%10d %d-%d (%10.5lf %10.5lf %10.5lf), ( %10.5lf %10.5lf %10.5lf) p_b_%4d=%10.5lf \n",
	       i,
	       (e_GOLM).pairs_bond[j][0],
	       (e_GOLM).pairs_bond[j][1],
	       crd[(e_GOLM).pairs_bond[j][0]*3+0],
	       crd[(e_GOLM).pairs_bond[j][0]*3+1],
	       crd[(e_GOLM).pairs_bond[j][0]*3+2],
	       crd[(e_GOLM).pairs_bond[j][1]*3+0],
	       crd[(e_GOLM).pairs_bond[j][1]*3+1],
	       crd[(e_GOLM).pairs_bond[j][1]*3+2],
	       j,
	       (e_GOLM).p_b[j]);
      }
      for (j=0;j<(e_GOLM).num_angl;++j) {
	printf("%10d %d-%d-%d (%10.5lf %10.5lf %10.5lf),  (%10.5lf %10.5lf %10.5lf),  (%10.5lf %10.5lf %10.5lf), p_a_%4d=%10.5lf \n",
	       i,
	       (e_GOLM).pairs_angl[j][0],
	       (e_GOLM).pairs_angl[j][1],
	       (e_GOLM).pairs_angl[j][2],
	       crd[(e_GOLM).pairs_angl[j][0]*3+0],
	       crd[(e_GOLM).pairs_angl[j][0]*3+1],
	       crd[(e_GOLM).pairs_angl[j][0]*3+2],
	       crd[(e_GOLM).pairs_angl[j][1]*3+0],
	       crd[(e_GOLM).pairs_angl[j][1]*3+1],
	       crd[(e_GOLM).pairs_angl[j][1]*3+2],
	       crd[(e_GOLM).pairs_angl[j][2]*3+0],
	       crd[(e_GOLM).pairs_angl[j][2]*3+1],
	       crd[(e_GOLM).pairs_angl[j][2]*3+2],
	       j,
	       (e_GOLM).p_a[j]);
      }
      //      for (j=0;j<(e_GOLM).num_dihe;++j) {
      for (j=1050;j<1052;++j) {
	
	///////////////////////////////////////
	ii=(e_GOLM).pairs_dihe[j][0];
	jj=(e_GOLM).pairs_dihe[j][1];
	kk=(e_GOLM).pairs_dihe[j][2];
	ll=(e_GOLM).pairs_dihe[j][3];
	for (k=0;k<3;++k) {
	  atom[0][k]=crd[ii*3+k];
	  atom[1][k]=crd[jj*3+k];
	  atom[2][k]=crd[kk*3+k];
	  atom[3][k]=crd[ll*3+k];
	}
	
	for (k=0;k<3;++k) {
	  vij[k] = atom[1][k]-atom[0][k];
	  vkj[k] = atom[1][k]-atom[2][k];
	  vkl[k] = atom[3][k]-atom[2][k];
	}
	lkj=sqrt(inprod(vkj,vkj,3));
	
  	outprod(vij,vkj,m);
	outprod(vkj,vkl,n);
	lm=sqrt(inprod(m,m,3));
	ln=sqrt(inprod(n,n,3));
	for (k=0;k<3;++k) {
	  m_n[k]=m[k]/lm;
	  n_n[k]=n[k]/ln;
	}
	
	dihed=acos(inprod(m_n,n_n,3));
	/////////////////////////////////



	printf("%10d %d-%d-%d-%d (%10.5lf %10.5lf %10.5lf),  (%10.5lf %10.5lf %10.5lf),  (%10.5lf %10.5lf %10.5lf),  (%10.5lf %10.5lf %10.5lf), p_d_%4d=%10.5lf Kd1=%10.5lf Kd2=%10.5lf equ=%10.5lf dihe=%10.5lf\n",
	       i,
	       (e_GOLM).pairs_dihe[j][0],
	       (e_GOLM).pairs_dihe[j][1],
	       (e_GOLM).pairs_dihe[j][2],
	       (e_GOLM).pairs_dihe[j][3],
	       crd[(e_GOLM).pairs_dihe[j][0]*3+0],
	       crd[(e_GOLM).pairs_dihe[j][0]*3+1],
	       crd[(e_GOLM).pairs_dihe[j][0]*3+2],
	       crd[(e_GOLM).pairs_dihe[j][1]*3+0],
	       crd[(e_GOLM).pairs_dihe[j][1]*3+1],
	       crd[(e_GOLM).pairs_dihe[j][1]*3+2],
	       crd[(e_GOLM).pairs_dihe[j][2]*3+0],
	       crd[(e_GOLM).pairs_dihe[j][2]*3+1],
	       crd[(e_GOLM).pairs_dihe[j][2]*3+2],
	       crd[(e_GOLM).pairs_dihe[j][3]*3+0],
	       crd[(e_GOLM).pairs_dihe[j][3]*3+1],
	       crd[(e_GOLM).pairs_dihe[j][3]*3+2],
	       j,
	       (e_GOLM).p_d[j],
	       (e_GOLM).Kd1,
	       (e_GOLM).Kd2,
	       (e_GOLM).dih_equ[j],
	       dihed
	       );
	
      }
    }
    // debug

    if (i%interval==0) {
      KE=KE/UNITT;
      T=KE/((3*numheavyatom)*k_B)*2.0;
      PEv=PEv/UNITT;
      KEv=KEv/UNITT;

      E_t=e_GOLM.p_t+KE+KEv+PEv;
      fprintf(outputfile,"%d %e %e %e %e %e %e %e\n",i+1,e_GOLM.p_t,KE,KEv,PEv,E_t,T,e_GOLM.p_natatt_t);
      fprintf(outputfile2,"E_t    = %e \n",E_t);
      fprintf(outputfile2,"KE     = %e \n",KE);
      fprintf(outputfile2,"KEv    = %e \n",KEv);
      fprintf(outputfile2,"PEv    = %e \n",PEv);
      fprintf(outputfile2,"p_tot  = %e \n",e_GOLM.p_t);
      fprintf(outputfile2,"p_nat  = %e \n",e_GOLM.p_natatt_t);
      fprintf(outputfile2,"p_rep  = %e \n",e_GOLM.p_repul_t);
      fprintf(outputfile2,"p_dih  = %e \n",e_GOLM.p_d_t);
      fprintf(outputfile2,"p_ang  = %e \n",e_GOLM.p_a_t);
      fprintf(outputfile2,"p_bon  = %e \n",e_GOLM.p_b_t);
      fprintf(outputfile2,"T      = %e \n",T);

      for (j=0;j<3;++j) COM[j]=0.0;
      for (j=0;j<numatom;++j) 
	for (k=0;k<3;++k) 
	  COM[k]+=mass[j]*crd[j*3+k]/summass;
      for (j=0;j<numatom;++j) 
	for (k=0;k<3;++k) 
	  crd[j*3+k]-=COM[k];

      for (j=0;j<numatom;++j)
	for (k=0;k<3;++k) crd_nc[j][k]=crd[j*3+k];

      avePE=(i*avePE+e_GOLM.p_t)/(i+1);
      varPE=(i*varPE+e_GOLM.p_t*e_GOLM.p_t)/(i+1);

      aveKE=(i*aveKE+KE)/(i+1);
      varKE=(i*varKE+KE*KE)/(i+1);

      aveT=(i*aveT+T)/(i+1);
      varT=(i*varT+T*T)/(i+1);

      myncL_put_crd_ene_MCD(nc_id_MCD,l,crd_nc,e,0.0);
      ++l;
    }
  }
  fclose(outputfile);
  fclose(outputfile2);
  nc_close((nc_id_MCD.ncid));

  rstfile=efopen(rstfilename,"w");
  fprintf(rstfile,"ACE\n");
  fprintf(rstfile,"%d\n",numatom);
  for (i=0;i<numatom;++i) {
    for (k=0;k<3;++k) fprintf(rstfile,"%e ",crd[i*3+k]);
    fprintf(rstfile,"\n");
  }
  fclose(rstfile);

  rstvelfile=efopen(rstvelfilename,"w");
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) fprintf(rstvelfile,"%e ",vel[i*3+j]);
    if (MODE==NVT) fprintf(rstvelfile,"%lf %lf",zeta,V_zeta);
    fprintf(rstvelfile,"\n");
  }
  fclose(rstvelfile);

  logfile=efopen(logfilename,"w");
  fprintf(logfile,"initial   coordinate file = %s\n",inputfilename);
  fprintf(logfile,"reference coordinate file = %s\n",refcrdfilename);
  fprintf(logfile,"topology             file = %s\n",refcrdfilename);
  fprintf(logfile,"information          file = %s\n",outputfilename);
  fprintf(logfile,"output    data       file = %s\n",outputfilename2);
  fprintf(logfile,"trajectory           file = %s\n",trjfilename);
  varPE=sqrt(varPE);
  fprintf(logfile,"Potential energy = %10.5lf kcal/mol +- %10.5lf kcal/mol\n",avePE,varPE);
  varKE=sqrt(varKE);
  fprintf(logfile,"Kinetic   energy = %10.5lf kcal/mol +- %10.5lf kcal/mol\n",aveKE,varKE);
  varT=sqrt(varT);
  fprintf(logfile,"Temperature      = %10.5lf K        +- %10.5lf K\n",aveT,varT);
  fclose(logfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename refcrdfilename parmfilename outputfilename outputfilename2 trjfilename\n",progname);
}
