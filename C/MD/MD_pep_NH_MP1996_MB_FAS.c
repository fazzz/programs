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
#include "TOPO.h"
#include "LA.h"

#include "netcdf_mineL.h"

#define NVT 1
#define NVE 0

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,a;
  int numatom,numheavyatom,numres,numstep=10000,interval=100;
  double dt=0.001,dt2,wdt2[3],wdt4[3];
  double pi;

  double ep=ep_natatt_hybrid;
  double de=1.0,d=1.0,d2;

  int vMode=OFF,MODE=NVT,NCmode=3,nibnum=3,criteria=6.5;

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
  int numnb,num14;

  double A,B,C,D;
  double p_d1_t,p_d2_t,*p_d1,*p_d2,**f_d1,**f_d2;
  double *dih_equ1,*dih_equ2;
  double e1,e2,p_MB,**f1,**f2,**f_MB,c1_c2,kai;
  double numa,numb;

  int ii,jj,kk,ll;

  double m[3],n[3],m_n[3],n_n[3],lm,ln;
  double vij[3],vkj[3],vkl[3];
  double lkj;
  double vijvkj,vklvkj;

  double atom[4][3];
  double dihed;


  double x[3];

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*refcrdfilename,*velfilename,*parmfilename;
  char *trjfilename,*outputfilename,*outputfilename2,*rstfilename="rstcrd",*rstvelfilename="rstvel";

  char *logfilename="MD_pep_NH_MP1996_GOLMAA_MB_PROTEINS2008.log";

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
    {"de",1,NULL,'d'},
    {"dh",1,NULL,'2'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"*hs:v:e:t:a:i:x:N:b:c:d:2:{:}:",long_opt,&opt_idx))!=-1) {
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
    case 'd':
      de=atof(optarg);
      break;
    case '2':
      d=atof(optarg);
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
  parmfilename     = *++argv;
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
  fscanf(inputfile,"%d",&a);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);

  refcrdfile=efopen(refcrdfilename,"r");
  getline(&line,&len,refcrdfile);
  fscanf(refcrdfile,"%d",&a);
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

  de=de*KT;
  d=d*KT;
  
  d2=d*d;

  summass=0.0;
  for (i=0;i<numatom;++i) summass+=mass[i];
  for (i=0;i<3;++i) COM[i]=0.0;
  for (i=0;i<numatom;++i) 
    for (j=0;j<3;++j) 
      COM[j]+=mass[i]*crd[i*3+j]/summass;
  for (i=0;i<numatom;++i) 
    for (j=0;j<3;++j) 
      crd[i*3+j]-=COM[j];

  //  ffL_set_calcffandforce(&e,&f);
  //  GOLMAA_PROTEINS2008_ff_set_calcff(&e_GOLM,refcrd,numatom,numres,e.parm.indexnb,e.parm.numnb,ep,nibnum,criteria);
  ffL_set_non_bonding_index_1(&numnb,&num14);
  e.parm.numnb=numnb;
  e.parm.num14=num14;
  e.parm.indexnb=(int *)gcemalloc(sizeof(int)*numnb*2);
  e.parm.index14=(int *)gcemalloc(sizeof(int)*num14*2);
  ffL_set_non_bonding_index_2(e.parm.indexnb,e.parm.index14);

  GOLMAA_PROTEINS2008_ff_set_calcff_b(&e_GOLM,refcrd,numatom,numres,e.parm.indexnb,e.parm.numnb,ep,nibnum,criteria);

  d2=d*d;
  d2=d2*k_B;
  de=de*k_B;
  dih_equ1=(double *)gcemalloc(sizeof(double)*2);
  dih_equ2=(double *)gcemalloc(sizeof(double)*2);
  dih_equ1[0]=-0.5*pi;
  dih_equ1[0]=0.5*pi;
  dih_equ2[1]=0.5*pi;
  dih_equ2[1]=-0.5*pi;
  p_d1=(double *)gcemalloc(sizeof(double)*(e_GOLM).num_dihe);
  p_d2=(double *)gcemalloc(sizeof(double)*(e_GOLM).num_dihe);
  f1=(double **)gcemalloc(sizeof(double *)*numatom);
  f2=(double **)gcemalloc(sizeof(double *)*numatom);
  f_d1=(double **)gcemalloc(sizeof(double *)*numatom);
  f_d2=(double **)gcemalloc(sizeof(double *)*numatom);
  f_MB=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      f1[i]=(double *)gcemalloc(sizeof(double)*3);
      f2[i]=(double *)gcemalloc(sizeof(double)*3);
      f_d1[i]=(double *)gcemalloc(sizeof(double)*3);
      f_d2[i]=(double *)gcemalloc(sizeof(double)*3);
      f_MB[i]=(double *)gcemalloc(sizeof(double)*3);
    }
  }

  (e_GOLM).p_b_t=GOLMAA_PROTEINS2008_ff_calcBOND(crd,numatom,(e_GOLM).p_b,(e_GOLM).f_b,(e_GOLM).Kb,(e_GOLM).bon_equ,(e_GOLM).pairs_bond,(e_GOLM).num_bond);

  (e_GOLM).p_a_t=GOLMAA_PROTEINS2008_ff_calcANGLE(crd,numatom,(e_GOLM).p_a,(e_GOLM).f_a,(e_GOLM).Ka,(e_GOLM).ang_equ,(e_GOLM).pairs_angl,(e_GOLM).num_angl);


  p_d1_t=GOLMAA_PROTEINS2008_ff_calcDIHE(crd,numatom,p_d1,f_d1,
					 (e_GOLM).Kd1,0.0,0.0,dih_equ1,(e_GOLM).pairs_dihe,
					 (e_GOLM).num_dihe,(e_GOLM).impindex);
  
  p_d2_t=GOLMAA_PROTEINS2008_ff_calcDIHE(crd,numatom,p_d2,f_d2,
					 (e_GOLM).Kd1,0.0,0.0,dih_equ2,(e_GOLM).pairs_dihe,
					 (e_GOLM).num_dihe,(e_GOLM).impindex);

  e1=p_d1_t+(e_GOLM).p_a_t+(e_GOLM).p_b_t;
  e2=p_d2_t+(e_GOLM).p_a_t+(e_GOLM).p_b_t;

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      f1[i][j]=f_d1[i][j]+(e_GOLM).f_a[i][j]+(e_GOLM).f_b[i][j];
      f2[i][j]=f_d2[i][j]+(e_GOLM).f_a[i][j]+(e_GOLM).f_b[i][j];
    }
  }

  A=0.5*(e1+e2+de);
  B=0.5*(e1-e2-de);
  C=sqrt(B*B+d2);
  D=e1-e2-de;
  
  p_MB=A-C;
  //  p_MB=e1;

  for (i=0;i<numatom;++i) 
    for (j=0;j<3;++j) 
      f_MB[i][j]=0.5*(f1[i][j]+f2[i][j])-0.25*D/C*(f1[i][j]-f2[i][j]);
      //      f_MB[i][j]=f1[i][j];

  c1_c2=(e1-p_MB)/d;
  kai=log(c1_c2);
  c1_c2=d/(e2+de-p_MB);
  kai=log(c1_c2);

  numa=(e1-p_MB)*(e2+de-p_MB);
  numb=d2;

  MD_Propagetor_NH_Single_set_MP1996(nc,dt,&dt2,wdt2,wdt4);
  myncL_create_def_MCD(trjfilename,numatom,&nc_id_MCD);
  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");
  for (i=0;i<numstep;++i) {
    if (MODE==NVT)
      KE=MD_Propagetor_NH_MP1998_MB_FASYS(crd,vel,mass,&zeta,&V_zeta,Q,NfKT,numatom,&KEv,&PEv,dt,dt2,nc,wdt4,wdt2,de,d2,&e_GOLM,dih_equ1,dih_equ2,&e1,&e2,&p_MB,f_MB);
    /********************************************************************************************/
    /* else 										        */
    /*   KE=MD_Propagetor_vV_NVE_GOLMAA_MB_PROTEINS2008(crd,vel,mass,numatom,dt,de,d2,&e_GOLM); */
    /********************************************************************************************/
    
    if (i%interval==0) {
      KE=KE/UNITT;
      T=KE/((3*numheavyatom)*k_B)*2.0;
      PEv=PEv/UNITT;
      KEv=KEv/UNITT;

      c1_c2=(e1-p_MB)/d;
      kai=log((e1-p_MB)/d);
      /********************************/
      /* c1_c2=d/(e2+de-p_MB);	      */
      /* if (c1_c2<0.0) c1_c2=-c1_c2; */
      /* kai=log(c1_c2);	      */
      /********************************/
      E_t=p_MB+KE+KEv+PEv;
      fprintf(outputfile,"%d %e %e %e %e %e %e %e %e\n",i+1,p_MB,KE,KEv,PEv,E_t,T,c1_c2,kai);

      for (j=0;j<(e_GOLM).num_dihe;++j) {
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

	dihed=inprod(m_n,n_n,3);
	if (dihed>=1.0) dihed=0.0;
	else if (dihed<=-1.0)dihed=pi;
	else dihed=acos(dihed);
	if (inprod(vij,n,3)>0) dihed=-dihed;
	if (dihed<0.0) dihed=2.0*pi+dihed;

	if (dihed>1.0*pi) {
	  dihed-=2.0*pi;
	}
	else if (dihed<-1.0*pi) {
	  dihed+=2.0*pi;
	}

	fprintf(outputfile2,"%8.3e ",dihed*180.0/pi);
      }
      fprintf(outputfile2,"\n");
       
      for (j=0;j<3;++j) COM[j]=0.0;
      for (j=0;j<numatom;++j) 
	for (k=0;k<3;++k) 
	  COM[k]+=mass[j]*crd[j*3+k]/summass;
      for (j=0;j<numatom;++j) 
	for (k=0;k<3;++k) 
	  crd[j*3+k]-=COM[k];
       
      for (j=0;j<numatom;++j)
	for (k=0;k<3;++k) crd_nc[j][k]=crd[j*3+k];
       
      avePE=(i*avePE+p_MB)/(i+1);
      varPE=(i*varPE+p_MB*p_MB)/(i+1);
       
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
