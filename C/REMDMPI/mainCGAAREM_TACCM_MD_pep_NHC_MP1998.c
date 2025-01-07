
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "REMDCGAA_TACCM_MPI_2.h"

#define ON  0
#define OFF 1

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int numstep=1000,interval=100;
  int numEX=1,numRE=2;
  int vMode=OFF;

  struct AADataforREMD AAdata;
  struct CGDataforREMD CGdata;
  struct AACGCommonDataforREMD Cdata;
  struct TACCMDataforREMD Zdata;

  int nc=1;                          
  double T0AA,T0CG,T0Z;
  double k_B=1.98723e-3;             
  double UNITT=418.4070;             
  double KTAA,KTCG,KBTZ,tau=0.01,tau2,pi;                    

  double dt,dt2,wdt2[3],wdt4[3];

  double *refcrd;
  double ep=0.3;
  int NCmode=3,nibnum=3;
  double criteria=6.5;
  int numnb,num14;

  double **theta;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*refcrdfilename,*parmfilename,*TACCMfilename;

  char *outputfilenameAAbase,*trjfilenameAAbase,outputfilenameAA[2000],trjfilenameAA[2000];
  char *outputfilenameCGbase,*trjfilenameCGbase,outputfilenameCG[2000],trjfilenameCG[2000];
  char *trjfileZbase,*trjfileThetaAAbase,*trjfileThetaCGbase;
  char trjfilenameZ[2000],trjfilenameThetaAA[2000],trjfilenameThetaCG[2000];
  FILE *inputfile,*refcrdfile,*parmfile,*TACCMfile;

  char *progname;
  int opt_idx=1;

  int my_rank,num_procs,tag = 0;  // MPI
  MPI_Status status;              // MPI

  // MPI
  MPI_Init(&argc, &argv);
  
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  // MPI

  pi=acos(-1.0);

  struct option long_opt[] = {
    {"vMode",1,NULL,'v'},
    {"nums",1,NULL,'s'},
    {"numRE",1,NULL,'n'},
    {"numEX",1,NULL,'e'},
    {"tau",1,NULL,'a'},
    {"int",1,NULL,'i'},
    {"dt",1,NULL,'x'},
    {"TAA",1,NULL,'T'},
    {"TCG",1,NULL,'t'},
    {"TZ",1,NULL,'B'},
    {"mZ",1,NULL,'m'},
    {"massX",1,NULL,'X'},
    {"ep",1,NULL,'p'},
    {"cutoff",1,NULL,'c'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hv:n:E:s:e:a:i:x:T:t:B:m:X:p:c:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 's':
      numstep=atoi(optarg); break;
    case 'i':
      interval=atoi(optarg);  break;
    case 'e':
      numEX=atoi(optarg); break;
    case 'n':
      numRE=atoi(optarg); break;
    case 'T':
      T0AA=atof(optarg);  break;
    case 't':
      T0CG=atof(optarg);  break;
    case 'B':
      T0Z=atof(optarg);break;
    case 'v':
      vMode=ON;        break;
    case 'a':
      tau=atof(optarg); break;
    case 'x':
      dt=atof(optarg);   break;
    case 'm':
      Zdata.massZ=atof(optarg);  break;
      //    case 'X':
      //      massflag=ON;  massX=atof(optarg);   break;
    case 'p':
      ep=atof(optarg); break;
    case 'c':
      //      printf("cutoff=%s\n",optarg);
      criteria=atof(optarg); break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);   exit(1);
    }
  }
  //  printf("cutoff=%8.4lf\n",criteria);
  progname=*argv;  argc-=optind;  argv+=optind;

  if (argc < 11) {
    USAGE(progname);
    exit(1);
  }
  inputfilename        = *argv;
  refcrdfilename       = *++argv;			   
  TACCMfilename        = *++argv;			 
  parmfilename         = *++argv;
  outputfilenameAAbase = *++argv;
  trjfilenameAAbase    = *++argv;
  outputfilenameCGbase = *++argv;
  trjfilenameCGbase    = *++argv;
  trjfileZbase         = *++argv;
  trjfileThetaAAbase   = *++argv;
  trjfileThetaCGbase   = *++argv;
  //    printf("yes 146 in main\n");  
  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  Cdata.numatom=AP.NATOM;
  j=0;  for (i=0;i<Cdata.numatom;++i)  if (strncmp(AP.IGRAPH[i],"H",1)==0)  ++j;
  Cdata.numheavyatom=Cdata.numatom-j;  Cdata.numres=AP.NRES;

  ///////////////// TACCM //////////////////////
  TACCMfile=efopen(TACCMfilename,"r");
  fscanf(TACCMfile,"%d",&Zdata.numZ);
  Zdata.pairs=(int **)gcemalloc(sizeof(int *)*Zdata.numZ);
  for (i=0;i<Zdata.numZ;++i) Zdata.pairs[i]=(int *)gcemalloc(sizeof(int)*5);
  for (i=0;i<Zdata.numZ;++i) {
    for (j=0;j<4;++j) fscanf(TACCMfile,"%d",&(Zdata.pairs[i][j]));
    fscanf(TACCMfile,"%d",&(Zdata.pairs[i][j]));
  }
  fclose(TACCMfile);

  //    printf("yes 157 in main\n");  
  Cdata.mass=(double *)gcemalloc(sizeof(double)*Cdata.numatom);
  for (i=0;i<Cdata.numatom;++i) Cdata.mass[i]=AP.AMASS[i];
  //    printf("yes 158 in main\n");  
  refcrd=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  //    printf("yes 159 in main\n");  
  AAdata.crd=(double **)gcemalloc(sizeof(double *)*numRE);  
  AAdata.vel=(double **)gcemalloc(sizeof(double *)*numRE);
  CGdata.crd=(double **)gcemalloc(sizeof(double *)*numRE);
  CGdata.vel=(double **)gcemalloc(sizeof(double *)*numRE);
  //    printf("yes 163 in main\n");  
  for (i=0;i<numRE;++i) {  
    AAdata.crd[i]=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
    AAdata.vel[i]=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
    CGdata.crd[i]=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
    CGdata.vel[i]=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  }
  //    printf("yes 170 in main\n");  
  MPI_Barrier(MPI_COMM_WORLD);
  Zdata.Z=(double **)gcemalloc(sizeof(double)*numRE);
  //    printf("yes 171 in main\n");  
  Zdata.velZ=(double **)gcemalloc(sizeof(double *)*numRE);
  //    printf("yes 172 in main\n");  
  MPI_Barrier(MPI_COMM_WORLD);
  for (i=0;i<numRE;++i) {
    //    theta[i]=(double *)gcemalloc(sizeof(double)*Zdata.numZ);
    Zdata.Z[i]=(double *)gcemalloc(sizeof(double)*Zdata.numZ);
    Zdata.velZ[i]=(double *)gcemalloc(sizeof(double)*Zdata.numZ);
  }
  //    printf("yes 178 in main\n");  
  Zdata.KZAA=(double *)gcemalloc(sizeof(double)*numRE);
  Zdata.KZCG=(double *)gcemalloc(sizeof(double)*numRE);
  //    printf("yes 181 in main\n");  
  AAdata.e=(struct potential *)gcemalloc(sizeof(struct potential)*numRE);
  AAdata.f=(struct force *)gcemalloc(sizeof(struct force)*numRE);
  CGdata.e_GOLM=(struct potential_GOLMAA_PROTEINS2008 *)gcemalloc(sizeof(struct potential_GOLMAA_PROTEINS2008)*numRE);
  AAdata.T=(double *)gcemalloc(sizeof(double)*numRE);
  CGdata.T=(double *)gcemalloc(sizeof(double)*numRE);
  Zdata.T=(double *)gcemalloc(sizeof(double)*numRE);
  //    printf("yes 179 in main\n");  
  inputfile=efopen(inputfilename,"r");
  CGAAREMDreadInputs(inputfile,Cdata.numatom,numRE,
		     AAdata.crd,AAdata.vel,
		     CGdata.crd,CGdata.vel,
		     Zdata.KZAA,Zdata.KZCG);
  fclose(inputfile);
  MPI_Barrier(MPI_COMM_WORLD);
  //    printf("yes 186 in main\n");  
  refcrdfile=efopen(refcrdfilename,"r");
  getline(&line,&len,refcrdfile);  fscanf(refcrdfile,"%d",&d);
  for (i=0;i<Cdata.numatom;++i) for (j=0;j<3;++j) fscanf(refcrdfile,"%lf",&refcrd[i*3+j]);
  fclose(refcrdfile);
  //    printf("yes 192 in main\n");  
  AAdata.zeta=(double *)gcemalloc(sizeof(double)*numRE);
  AAdata.V_zeta=(double *)gcemalloc(sizeof(double)*numRE);
  //    printf("yes 195 in main\n");  
  CGdata.zeta=(double *)gcemalloc(sizeof(double)*numRE);
  CGdata.V_zeta=(double *)gcemalloc(sizeof(double)*numRE);
  //    printf("yes 198 in main\n");  
  Zdata.zetaZ=(double *)gcemalloc(sizeof(double)*numRE);
  Zdata.V_zetaZ=(double *)gcemalloc(sizeof(double)*numRE);
  //    printf("yes 201 in main\n");  
  MPI_Barrier(MPI_COMM_WORLD);
  if ( vMode==OFF ) {
    for (i=0;i<numRE;++i) {
      MD_Generate_inivelo(AAdata.vel[i],Cdata.mass,Cdata.numatom,k_B*T0AA*UNITT);
      MPI_Barrier(MPI_COMM_WORLD);
      //      printf("yes 206 in main\n");  
      MD_Generate_inivelo(CGdata.vel[i],Cdata.mass,Cdata.numatom,k_B*T0CG*UNITT);
      MPI_Barrier(MPI_COMM_WORLD);
      //      printf("yes 208 in main\n");  
      TACCM_MD_Generate_inivelo(Zdata.velZ[i],Zdata.massZ,Zdata.numZ,k_B*T0Z*UNITT);
      MPI_Barrier(MPI_COMM_WORLD);
      //      printf("yes 210 in main\n");  
      for (j=0;j<Cdata.numatom;++j) 
	if (strncmp(AP.IGRAPH[j],"H",1)==0)  
	  for (k=0;k<3;++k) 
	    CGdata.vel[i][j*3+k]=0.0;

      //      printf("yes 216 in main\n");  
      AAdata.zeta[i]=0.0;      AAdata.V_zeta[i]=0.0;
      CGdata.zeta[i]=0.0;      CGdata.V_zeta[i]=0.0;
      Zdata.zetaZ[i]=0.0;      Zdata.V_zetaZ[i]=0.0;
    }
  }  
  //    printf("yes 216 in main\n");  

  //  theta=(double **)gcemalloc(sizeof(double *)*numRE);
  MPI_Barrier(MPI_COMM_WORLD);

  for (i=0;i<numRE;++i) {
    TACCM_CTheta(AAdata.crd[i],Cdata.numatom,Zdata.Z[i],Zdata.numZ,Zdata.pairs,pi);
    //    TACCM_CTheta(CGdata.crd[i],Cdata.numatom,theta[i],Zdata.numZ,Zdata.pairs,pi);
  }
  //  for (i=0;i<numRE;++i) for (j=0;j<Zdata.numZ;++j) Zdata.Z[i][j]=theta[i][j];

  AAdata.avePE=(double *)gcemalloc(sizeof(double)*numRE);  AAdata.varPE=(double *)gcemalloc(sizeof(double)*numRE);
  AAdata.aveKE=(double *)gcemalloc(sizeof(double)*numRE);  AAdata.varKE=(double *)gcemalloc(sizeof(double)*numRE);
  AAdata.aveT=(double *)gcemalloc(sizeof(double)*numRE);   AAdata.varT=(double *)gcemalloc(sizeof(double)*numRE);

  CGdata.avePE=(double *)gcemalloc(sizeof(double)*numRE);  CGdata.varPE=(double *)gcemalloc(sizeof(double)*numRE);
  CGdata.aveKE=(double *)gcemalloc(sizeof(double)*numRE);  CGdata.varKE=(double *)gcemalloc(sizeof(double)*numRE);
  CGdata.aveT=(double *)gcemalloc(sizeof(double)*numRE);   CGdata.varT=(double *)gcemalloc(sizeof(double)*numRE);

  Zdata.avePEZ=(double *)gcemalloc(sizeof(double)*numRE);  Zdata.varPEZ=(double *)gcemalloc(sizeof(double)*numRE);
  Zdata.aveKEZ=(double *)gcemalloc(sizeof(double)*numRE);  Zdata.varKEZ=(double *)gcemalloc(sizeof(double)*numRE);
  Zdata.aveTZ=(double *)gcemalloc(sizeof(double)*numRE);   Zdata.varTZ=(double *)gcemalloc(sizeof(double)*numRE);
  //    printf("yes 275 in main\n");    

  tau=tau/2.0/pi;  tau2=tau*tau;
  
  KBTZ=k_B*T0Z;  Zdata.NfKTZ=(Zdata.numZ+1)*KBTZ*UNITT;  Zdata.QZ=tau2*KBTZ*UNITT*Zdata.numZ;

  KTAA=k_B*T0AA;  AAdata.NfKT=(3.0*Cdata.numatom+1)*KTAA*UNITT;  AAdata.Q=tau2*KTAA*UNITT*(3.0*Cdata.numatom);

  KTCG=k_B*T0CG;  CGdata.NfKT=(3.0*Cdata.numheavyatom+1)*KTCG*UNITT;  
  CGdata.Q=tau2*KTCG*UNITT*(3.0*Cdata.numheavyatom);
  //    printf("yes 285 in main\n");    
  MPI_Barrier(MPI_COMM_WORLD);
  //    printf("yes 287 in main\n");    
  for (i=0;i<numRE;++i)  {
    ffL_set_calcffandforce(&(AAdata.e[i]),&(AAdata.f[i]));
  }
  MPI_Barrier(MPI_COMM_WORLD);
  //    printf("yes 290 in main\n");    
  for (i=0;i<numRE;++i)  {
    ffL_set_non_bonding_index_1(&numnb,&num14);
    AAdata.e[i].parm.numnb=numnb;    AAdata.e[i].parm.num14=num14;
    AAdata.e[i].parm.indexnb=(int *)gcemalloc(sizeof(int)*numnb*2);
    AAdata.e[i].parm.index14=(int *)gcemalloc(sizeof(int)*num14*2);
    ffL_set_non_bonding_index_2(AAdata.e[i].parm.indexnb,AAdata.e[i].parm.index14);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  //  printf("yes 300 in main\n");    
  for (i=0;i<numRE;++i)  {
    //    printf("numatom=%4d numres=%4d refcrd[10]=%8.4lf criteria=%8.4lf\n",
    //	   Cdata.numatom,Cdata.numres,refcrd[0],criteria);

    GOLMAA_PROTEINS2008_ff_set_calcff_b(&(CGdata.e_GOLM[i]),refcrd,Cdata.numatom,Cdata.numres,
					AAdata.e[i].parm.indexnb,AAdata.e[i].parm.numnb,ep,nibnum,criteria);
    //    printf("yes 305 in main\n");    
  }
  for (i=0;i<numRE;++i)  ffL_set_calcffandforce(&(AAdata.e[i]),&(AAdata.f[i]));
  //  printf("yes 308 in main\n");    

  MD_Propagetor_NH_Single_set_MP1996(nc,dt,&dt2,wdt2,wdt4);
  //  printf("yes 323 in main\n");

  AAdata.nc_id_MCD=(struct my_netcdf_out_id_MCD *)gcemalloc(sizeof(struct my_netcdf_out_id_MCD)*numRE);
  AAdata.outputfile=(FILE **)gcemalloc(sizeof(FILE *)*numRE);

  CGdata.nc_id_MCD=(struct my_netcdf_out_id_MCD *)gcemalloc(sizeof(struct my_netcdf_out_id_MCD)*numRE);
  CGdata.outputfile=(FILE **)gcemalloc(sizeof(FILE *)*numRE);

  Zdata.trjfileZ=(FILE **)gcemalloc(sizeof(FILE *)*numRE);  
  Zdata.trjfilThetaAA=(FILE **)gcemalloc(sizeof(FILE *)*numRE);
  Zdata.trjfilThetaCG=(FILE **)gcemalloc(sizeof(FILE *)*numRE);
  //  printf("numRE=%d\n",numRE);
  sprintf(outputfilenameAA,"%s_%d",outputfilenameAAbase,my_rank+1);
  //  printf("%s\n",outputfilenameAA);
  sprintf(trjfilenameAA,"%s_%d",trjfilenameAAbase,my_rank+1);
  //  printf("%s\n",trjfilenameAA);

  sprintf(outputfilenameCG,"%s_%d",outputfilenameCGbase,my_rank+1);
  //  printf("%s\n",outputfilenameCG);
  sprintf(trjfilenameCG,"%s_%d",trjfilenameCGbase,my_rank+1);
  //  printf("%s\n",trjfilenameCG);

  sprintf(trjfilenameZ,"%s_%d",trjfileZbase,my_rank+1);
  //  printf("%s\n",trjfilenameZ);
  sprintf(trjfilenameThetaAA,"%s_%d",trjfileThetaAAbase,my_rank+1);
  //  printf("%s\n",trjfilenameThetaAA);
  sprintf(trjfilenameThetaCG,"%s_%d",trjfileThetaCGbase,my_rank+1);
  //  printf("%s\n",trjfilenameThetaCG);

  MPI_Barrier(MPI_COMM_WORLD);

  //  printf("yes 335 in main\n");
  myncL_create_def_MCD(trjfilenameAA,Cdata.numatom,&(AAdata.nc_id_MCD[my_rank]));
  AAdata.outputfile[my_rank]=efopen(outputfilenameAA,"w");

  myncL_create_def_MCD(trjfilenameCG,Cdata.numatom,&(CGdata.nc_id_MCD[my_rank]));
  CGdata.outputfile[my_rank]=efopen(outputfilenameCG,"w");
  MPI_Barrier(MPI_COMM_WORLD);
  //  printf("yes 340 in main\n");
  Zdata.trjfileZ[my_rank]=efopen(trjfilenameZ,"w");
  //  printf("yes 347 in main\n");
  Zdata.trjfilThetaAA[my_rank]=efopen(trjfilenameThetaAA,"w");
  //  printf("yes 349 in main\n");
  Zdata.trjfilThetaCG[my_rank]=efopen(trjfilenameThetaCG,"w");
  MPI_Barrier(MPI_COMM_WORLD);
  //  printf("yes 351 in main\n");
  MPI_Barrier(MPI_COMM_WORLD);
  if ( num_procs != numRE ) {    printf("condition error\n");    exit(1);  }
  //  printf("yes 310 in main\n");  
  //  exit(1);
  MPI_CGAATREM_TACCM_MD_pep_NHC_MP1998(my_rank,num_procs,tag,&status,numRE,numEX,
				       AAdata,CGdata,Zdata,Cdata,
				       T0AA,T0CG,T0Z,numstep,interval,
				       dt,dt2,wdt2,wdt4,nc,UNITT,k_B,tau,pi);
  
  fclose(AAdata.outputfile[my_rank]);  
  //  nc_close(AAdata.nc_id_MCD[my_rank]);
  fclose(CGdata.outputfile[my_rank]);  
  //  nc_close(CGdata.nc_id_MCD[my_rank]);
  fclose(Zdata.trjfileZ[my_rank]);  
  fclose(Zdata.trjfilThetaAA[my_rank]); 
  fclose(Zdata.trjfilThetaCG[my_rank]);
  
  MPI_Finalize();
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename refcrdfilename TACCMfilename parmfilename outputfilenameAAbase trjfilenameAAbase outputfilenameCGbase trjfilenameCGbase trjfileZbase trjfileThetaZ\n",progname);
}
