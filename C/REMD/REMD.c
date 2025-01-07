
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "REMD.h"
#include "MC.h"
#include "RAND.h"
#include "SHELL.h"
#include "PT.h"
#include "EF.h"
#include "FF.h"

#define kb 1.98723e-3

#define ON 1
#define OFF 0
#define ACEP 1

int Metropolisdb(double delta);
double CD_SDFF_sub(char *crdname,char *parmname,char *indexfilename,int numnb, int num14);

int HreplicaExchange(int numjudge,int numreplica,char *crd,char *vel, char *top,char *clt,char *pn,char *inpmd[MAXREP],char *rstin[MAXREP],char *pepca[MAXREP],char *dirbase,double temp,int waittime ) {
  int i,j,n,nr,num,judgeflag,dummy;
  int *npepca;
  double delta=0.0;
  char out[100],rst[100],rsttrj[100],rstvel[100],rstnew[100];
  char diri[100],dirj[100];
  char crdnamei[100],crdnamej[100];
  char qsubsh[100];

  char log[100];
  FILE *logf;

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  npepca=(int *)gcemalloc(sizeof(int)*numreplica);

  for (nr=0;nr<numreplica;++nr) {
    npepca[nr]=nr;
    sprintf(diri,"%s/%d",dirbase,nr+1);
    sprintf(qsubsh,"qsub_%d_0.sh",nr+1);
    sprintf(out,"%s/%s_remd_%d",diri,pn,nr+1);
    mkdir(diri);
    runMD(inpmd[nr],crd,vel,top,clt,out,out,pn,diri,pepca[npepca[nr]],qsubsh);
  }
  wait("remd",waittime);
  
  for (i=0;i<numjudge;++i) {
    if (i%2==0) judgeflag=ON;
    else judgeflag=OFF;
  
    for (nr=0;nr<numreplica;++nr) {
      num=nr+2;
      if (num>numreplica)
	num-=numreplica;
      sprintf(diri,"%s/%d",dirbase,nr+1);
      sprintf(dirj,"%s/%d",dirbase,num);
      sprintf(out,"%s/%s_remd_%d",diri,pn,nr+1);
      sprintf(qsubsh,"qsub_%d_%d.sh",nr+1,i+1);

      if (judgeflag==OFF) {
	sprintf(log,"%s/log_rex.txt",diri);
	logf=efopen(log,"a");
	fprintf(logf,"- %d\n",npepca[nr]+1);
	fclose(logf);
      }
      else {
	sprintf(crdnamei,"%s/%s_remd_%d.rst",diri,pn,nr+1);
	sprintf(crdnamej,"%s/%s_remd_%d.rst",dirj,pn,num);
	delta=calcdelta(crdnamei,crdnamej,pepca[npepca[nr]],pepca[npepca[(num-1)]],top,temp);
  	if (Metropolisdb(delta)==ACEP) {
          dummy=npepca[nr];
	  npepca[nr]=npepca[(num-1)];
	  npepca[(num-1)]=dummy;
	  sprintf(log,"%s/log_rex.txt",diri);
	  logf=efopen(log,"a");
	  fprintf(logf,"1 %d\n",npepca[nr]+1);
	  fclose(logf);
	}
  	else {
	  sprintf(log,"%s/log_rex.txt",diri);
	  logf=efopen(log,"a");
	  fprintf(logf,"0 %d\n",npepca[nr]+1);
	  fclose(logf);
	}
      }
  
      sprintf(rsttrj,"%s_remd_%d.rst",pn,nr+1);
      sprintf(rstvel,"%s_remd_%d.rve",pn,nr+1);
      runMD(rstin[nr],rsttrj,rstvel,top,clt,out,out,pn,diri,pepca[npepca[nr]],qsubsh);
      if (judgeflag==OFF) judgeflag=ON;
      else judgeflag=OFF;
    }
    wait("remd",waittime);
    printf("%d th simulations done\n",i);
  }
  return 0;
}

int  runMD(char *inpmd,char *crd, char *vel,char *top, char *clt,char *rst, char *out, char *pn, char *dir, char *pepca,char *qsubfilename) {
  char qsubsh[100];

  makeqsub(inpmd,crd,vel,top,clt,rst,out,pn,dir,pepca,qsubfilename);

  sprintf(qsubsh,"chmod +x %s;qsub %s;",qsubfilename,qsubfilename);
  chdir(dir);
  system(qsubsh);

}

double calcdelta(char *crdnamei,char *crdnamej,char *eignamei,char *eignamej,char *parmname,double temp) {
  int i,n;

  double *ele,*ALJ,*BLJ;
  double *p_ei,*p_ej,*p_LJi,*p_LJj,p_ii=0.0,p_ij=0.0,p_jj=0.0,p_ji=0.0,*f_e,*f_LJ;
  double *ui,*uj;
  int numnb,*indexnb;
  int numatom,numpara;
  double *crdi,*crdj;
  double facti,factj;
  int flagp,flagf;

  char *line;
  size_t len=0;

  FILE *fcrdi,*fcrdj,*parmtop,*eigi,*eigj;

  parmtop=efopen(parmname,"r");
  readParmtop(parmtop);
  numatom=AP.NATOM;
  numpara=AP.NTYPES*(AP.NTYPES+1)/2;
  numnb=ff_set_numnb();
  fclose(parmtop);

  indexnb=(int *)gcemalloc(sizeof(int)*numnb*2);
  ele=(double *)gcemalloc(sizeof(double)*numatom);
  ALJ=(double *)gcemalloc(sizeof(double)*numatom*numatom);
  BLJ=(double *)gcemalloc(sizeof(double)*numatom*numatom);

  crdi=(double *)gcemalloc(sizeof(double)*numatom*3);
  crdj=(double *)gcemalloc(sizeof(double)*numatom*3);

  p_ei=(double *)gcemalloc(sizeof(double)*numnb);
  p_LJi=(double *)gcemalloc(sizeof(double)*numnb);
  p_ej=(double *)gcemalloc(sizeof(double)*numnb);
  p_LJj=(double *)gcemalloc(sizeof(double)*numnb);

  ui=(double *)gcemalloc(sizeof(double)*numnb*2);
  uj=(double *)gcemalloc(sizeof(double)*numnb*2);

  ff_set_NB_PARM(ele,ALJ,BLJ,numatom);
  n=ff_set_NB_index(indexnb,numnb,numatom);

  fcrdi=efopen(crdnamei,"r");
  getline(&line,&len,fcrdi);
  getline(&line,&len,fcrdi);
  io_scanconf(fcrdi,numatom,crdi,'x');
  fclose(fcrdi);
  
  fcrdj=efopen(crdnamej,"r");
  getline(&line,&len,fcrdj);
  getline(&line,&len,fcrdj);
  io_scanconf(fcrdj,numatom,crdj,'x');
  fclose(fcrdj);
  
  eigi=efopen(eignamei,"r");
  for (i=0;i<numnb;++i)
    fscanf(eigi,"%lf",&ui[i*2]);
  for (i=0;i<numnb;++i)
    fscanf(eigi,"%lf",&ui[i*2+1]);
  fscanf(eigi,"%lf",&facti);
  fclose(eigi);
  eigj=efopen(eignamej,"r");
  for (i=0;i<numnb;++i)
    fscanf(eigj,"%lf",&uj[i*2]);
  for (i=0;i<numnb;++i)
    fscanf(eigj,"%lf",&uj[i*2+1]);
  fscanf(eigj,"%lf",&factj);
  fclose(eigj);

  ff_calcFFNB(ele,ALJ,BLJ,p_ei,p_LJi,f_e,f_LJ,numnb,indexnb,numatom,crdi,1,0);
  ff_calcFFNB(ele,ALJ,BLJ,p_ej,p_LJj,f_e,f_LJ,numnb,indexnb,numatom,crdj,1,0);

  for (i=0;i<numnb;++i) {
    p_ii+=ui[i*2]*p_ei[i]*facti;
    p_ii+=ui[i*2+1]*p_LJi[i]*facti;
  
    p_jj+=uj[i*2]*p_ej[i]*factj;
    p_jj+=uj[i*2+1]*p_LJj[i]*factj;
  
    p_ij+=uj[i*2]*p_ei[i]*factj;
    p_ij+=uj[i*2+1]*p_LJi[i]*factj;
  
    p_ji+=ui[i*2]*p_ej[i]*facti;
    p_ji+=ui[i*2+1]*p_LJj[i]*facti;
  }

  return 1.0/(kb*temp)*(p_ij+p_ji-p_ii-p_jj);
}

int Metropolisdb(double delta) {
  double expdelta; // exp(-delta)
  double u;

  char logn[20];
  FILE *log;

  expdelta=exp(-1.0*delta);

  u=genrand_real2();

  log=efopen("log_mc.txt","a");
  fprintf(log,"%lf\n",u);
  fclose(log);

  /*********************************************/
  /* sprintf(logn,"%s/log_remd.txt",filename); */
  /* log=efopen(logn,"a");		       */
  /* fprintf(log,"de=%lf u=%lf ",delta,u);     */
  /* fclose(log);			       */
  /*********************************************/

  if (expdelta >= 1.0)
    return 1; // ACEPT
  else if (/*(*/u/*=genrand_real2())*/ <= expdelta)
    return 1; // REJECT
  else
    return 0; // REJECT



}

int HREMD_vac_SDFF(int numexchange,int numreplica,char *minin,char *rexin, char *samin, char **crd,char **top,char *pn,char *dirbase,char *indexfilename,int numnb, int num14,double temp,int waittime,int numstep,int numatom) {
  int i,j,n,nr,num,judgeflag,nrex;
  char dummy[100];
  double delta=0.0;
  struct rep_HREMD_vac_SDFF *rep;

  char diri[100],dirmin[100],dirrex[100],dirsam[100];
  char minrst[100],minout[100];
  char rexrst[100],rexout[100],rextrj[100],rexvel[100];
  char samrst[100],samtrj[100],samvel[100],samout[100],samcrd[100];
  char *top1,*top2,rstnr[100],rstnrex[100];

  char log[100];
  FILE *logf;

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  rep=(struct rep_HREMD_vac_SDFF *)gcemalloc(sizeof(struct rep_HREMD_vac_SDFF)*numreplica);

  // make dirs
  for (nr=0;nr<numreplica;++nr) {
    sprintf(rep[nr].dir,"%s/rep%d",dirbase,nr+1);
    sprintf(rep[nr].dirmin,"%s/rep%d/min",dirbase,nr+1);
    sprintf(rep[nr].dirrex,"%s/rep%d/rex",dirbase,nr+1);
    sprintf(rep[nr].dirsam,"%s/rep%d/sam",dirbase,nr+1);
    rep[nr].parmtop=top[nr];
    mkdir(rep[nr].dir);
    mkdir(rep[nr].dirmin);
    mkdir(rep[nr].dirrex);
    mkdir(rep[nr].dirsam);
  }

  //minimization
  for (nr=0;nr<numreplica;++nr) {
    sprintf(minrst,"%s_min.rst",pn);
    sprintf(minout,"%s_min.out",pn);
    runMD_Amber(minin,crd[nr],crd[nr],"dummy",rep[nr].parmtop,minrst,"dummy",minout,pn,rep[nr].dirmin,"qsub_min.sh","min");
  }
  wait("min",waittime);
  printf("minimization is done!\n");

  //relax-calculation
  for (nr=0;nr<numreplica;++nr) {
    sprintf(minrst,"../min/%s_min.rst",pn);
    sprintf(rexrst,"%s_rex.rst",pn);
    sprintf(rexout,"%s_rex.out",pn);
    sprintf(rextrj,"%s_rex.trj",pn);
    sprintf(rexvel,"%s_rex.vel",pn);
    runMD_Amber(rexin,minrst,"dummy",rexvel,rep[nr].parmtop,rexrst,rextrj,rexout,pn,rep[nr].dirrex,"qsub_rex.sh","rex");
  }
  wait("rex",waittime);
  printf("relax-calculation is done!\n");

  //remd-1st-calculation
  i=0;
  for (nr=0;nr<numreplica;++nr) {
    sprintf(rexrst,"../rex/%s_rex.rst",pn);
    sprintf(samrst,"%s_sam_%d.rst",pn,1);
    sprintf(samout,"%s_sam_%d.out",pn,1);
    sprintf(samtrj,"%s_sam_%d.trj",pn,1);
    sprintf(samvel,"%s_sam_%d.vel",pn,1);
    runMD_Amber(samin,rexrst,"dummy",samvel,rep[nr].parmtop,samrst,samtrj,samout,pn,rep[nr].dirsam,"qsub_sam.sh","sam");
    printf("%d-th simulations of %d-th replica is submitted\n",i+1,nr+1);
  }
  wait("sam",waittime);
  printf("1st-sampling-calculation is done!\n");

  //remd-calculation
  for (i=1;i<numexchange;++i) {
    if (i%2==1) judgeflag=ON;
    else judgeflag=OFF;
    
    for (nr=0;nr<numreplica;++nr)
      sprintf(rep[nr].rstname,"%s/%s_sam_%d.rst",rep[nr].dirsam,pn,i);

    for (nr=0;nr<numreplica;++nr) {
      if (judgeflag==OFF || nr==numreplica-1) {
  	sprintf(log,"%s/log_remd.txt",rep[nr].dirsam);
  	logf=efopen(log,"a");
  	fprintf(logf,"- %s\n",rep[nr].parmtop);
  	fclose(logf);
      }
      else {
	/********************************/
        /* if (nr==numreplica-1)        */
	/*   nrex=0;		        */
	/* else			        */
        /********************************/
	nrex=nr+1;
	sprintf(rstnr,"%s/%s_sam_%d.rst",rep[nr].dirsam,pn,i);
	sprintf(rstnrex,"%s/%s_sam_%d.rst",rep[nrex].dirsam,pn,i);
  	top1=rep[nr].parmtop;
  	top2=rep[nrex].parmtop;
        delta=CD_SDFF(rstnr,rstnrex,top1,top2,indexfilename,numnb,num14,temp);
  	
  	if (Metropolisdb(delta)==ACEP) {
	  /*****************************************************/
          /* strcpy(dummy,rep[nr].rstname);		       */
	  /* strcpy(rep[nr].rstname,rep[nrex].rstname);	       */
  	  /* strcpy(rep[nrex].rstname,dummy);		       */
          /*****************************************************/
          strcpy(dummy,rep[nr].parmtop);
	  strcpy(rep[nr].parmtop,rep[nrex].parmtop);
  	  strcpy(rep[nrex].parmtop,dummy);
  	  sprintf(log,"%s/log_remd.txt",rep[nr].dirsam);
  	  logf=efopen(log,"a");
  	  fprintf(logf,"1 %s %lf %lf\n",rep[nr].parmtop,delta,exp(-1.0*delta));
  	  fclose(logf);
  	}
  	else {
  	  sprintf(log,"%s/log_remd.txt",rep[nr].dirsam);
  	  logf=efopen(log,"a");
  	  fprintf(logf,"0 %s %lf %lf\n",rep[nr].parmtop,delta,exp(-1.0*delta));
  	  fclose(logf);
        }
      }
  
      sprintf(samrst,"%s_sam_%d.rst",pn,i+1);
      sprintf(samout,"%s_sam_%d.out",pn,i+1);
      sprintf(samtrj,"%s_sam_%d.trj",pn,i+1);
      sprintf(samvel,"%s_sam_%d.vel",pn,i+1);
      runMD_Amber(samin,rep[nr].rstname,"dummy",samvel,rep[nr].parmtop,samrst,samtrj,samout,pn,rep[nr].dirsam,"qsub_sam.sh","sam");
      printf("%d-th simulations of %d-th replica is submitted\n",i+1,nr+1);
      if (judgeflag==OFF) judgeflag=ON;
      else judgeflag=OFF;
    }
    wait("sam",waittime);
    printf("%d-th simulations done\n",i+1);
  }

  //  gatherdata(numexchange,numreplica,numstep,numatom,dirbase,pn);


  return 0;
}

int  runMD_Amber(char *inpmd,char *crd,char *ref, char *vel,char *top, char *rst, char *trj, char *out, char *pn, char *dir, char *qsubfilename, char *simtype) {
  char qsubsh[100];

  makeqsub_serial_Amber(inpmd,crd,ref,vel,top,rst,trj,out,pn,dir,qsubfilename,simtype);
  
  sprintf(qsubsh,"chmod +x %s;qsub %s;",qsubfilename,qsubfilename);
  chdir(dir);
  system(qsubsh);

}

double CD_SDFF(char *crdnamei,char *crdnamej,char *parmname1,char *parmname2,char *indexfilename,int numnb, int num14,double temp) {
  double p1,p2;

  p1=CD_SDFF_sub(crdnamei,parmname1,indexfilename,numnb,num14);
  p2=CD_SDFF_sub(crdnamei,parmname2,indexfilename,numnb,num14);

  return 1.0/(kb*temp)*(p2-p1);
}

double CD_SDFF_sub(char *crdname,char *parmname,char *indexfilename,int numnb, int num14) {
  int i,j,k,f;
  double p_t;
  double *ele,*ALJ,*BLJ;
  double *p_e,*p_LJ,*p_d,*p_a,*p_b;
  double p_e_t,p_LJ_t,p_e_14_t,p_LJ_14_t,p_d_t,p_a_t,p_b_t;
  double *f_e,*f_LJ,*n_d;
  double *p_e_14,*p_LJ_14;
  double *f_e_14,*f_LJ_14;
  double *crd;
  char *line,dummy;
  size_t len=0;

  int *indexnb,*index14;
  int numatom,numpara;
  FILE *inputfile,*parmtop1,*parmtop2,*indexfile;

  parmtop1=efopen(parmname,"r");
  readParmtop(parmtop1);
  fclose(parmtop1);
  numatom=AP.NATOM;
  numpara=AP.NTYPES*(AP.NTYPES+1)/2;

  indexnb=(int *)emalloc(sizeof(int)*numnb*2);
  ele=(double *)emalloc(sizeof(double)*numatom);
  ALJ=(double *)emalloc(sizeof(double)*numatom*numatom);
  BLJ=(double *)emalloc(sizeof(double)*numatom*numatom);
  index14=(int *)emalloc(sizeof(int)*num14*2);

  ff_set_NB_PARM(ele,ALJ,BLJ,numatom);

  indexfile=efopen(indexfilename,"r");
  for (i=0;i<numnb;++i) {
    fscanf(indexfile,"%d",&f);indexnb[i*2]=f-1;
    fscanf(indexfile,"%s",&dummy);
    fscanf(indexfile,"%d",&f);indexnb[i*2+1]=f-1;
  }
  for (i=0;i<num14;++i) {
    fscanf(indexfile,"%d",&f);index14[i*2]=f-1;
    fscanf(indexfile,"%s",&dummy);
    fscanf(indexfile,"%d",&f);index14[i*2+1]=f-1;
  }
  fclose(indexfile);

  inputfile=efopen(crdname,"r");

  crd=(double *)ecalloc(sizeof(double),numatom*3);
  p_e=(double *)ecalloc(sizeof(double),numnb);
  p_LJ=(double *)ecalloc(sizeof(double),numnb);

  p_e_14=(double *)ecalloc(sizeof(double),numnb);
  p_LJ_14=(double *)ecalloc(sizeof(double),numnb);
  p_d=(double *)ecalloc(sizeof(double),(AP.NPHIH+AP.MPHIA));
  p_a=(double *)ecalloc(sizeof(double),(AP.NTHETH+AP.MTHETA));
  p_b=(double *)ecalloc(sizeof(double),(AP.NBONH+AP.MBONA));

  io_scanconf_Amber_rst(inputfile,crd);
  ff_calcFFNB(ele,ALJ,BLJ,p_e,p_LJ,f_e,f_LJ,numnb,indexnb,numatom,crd,2,0);
  ff_calcFFNB(ele,ALJ,BLJ,p_e_14,p_LJ_14,f_e_14,f_LJ_14,num14,index14,numatom,crd,2,0);
  ff_calcDIHE(p_d,n_d,crd,1,0,0);
  ff_calcANGLE(p_a,crd);
  ff_calcBOND(p_b,crd);

  p_t=0.0;
  p_e_t=0.0;
  p_LJ_t=0.0;
  p_e_14_t=0.0;
  p_LJ_14_t=0.0;
  p_d_t=0.0;
  p_a_t=0.0;
  p_b_t=0.0;
  
  for (j=0;j<numnb;++j) {
    p_t+=p_e[j]+p_LJ[j];
    p_e_t+=p_e[j];
    p_LJ_t+=p_LJ[j];
  }
  for (j=0;j<num14;++j) {
    p_t+=1.0/1.2*p_e_14[j]+0.5*p_LJ_14[j];
    p_e_14_t+=1.0/1.2*p_e_14[j];
    p_LJ_14_t+=0.5*p_LJ_14[j];
  }
  for (j=0;j<AP.NPHIH+AP.MPHIA;++j) {
    p_t+=p_d[j];
    p_d_t+=p_d[j];
  }
  for (j=0;j<AP.NTHETH+AP.MTHETA;++j) {
    p_t+=p_a[j];
    p_t+=p_b[j];
    p_a_t+=p_a[j];
    p_b_t+=p_b[j];
  }

  free(crd);
  free(p_e);
  free(p_LJ);
  free(p_e_14);
  free(p_LJ_14);
  free(p_d);
  free(p_a);
  free(p_b);

  fclose(inputfile);

  indexnb;
  free(ele);
  free(ALJ);
  free(BLJ);
  free(index14);

  return p_t;
}

void gatherdata(int numexchange,int numreplica,int numstep,int numatom,char *dirbase,char *pn) {
  int nr,i;
  double *trj,*vel;
  char dirdata[100],gtrjfilename[100],gvelfilename[100],trjfilename[100],velfilename[100];
  FILE *trjfile,*velfile,*gtrjfile,*gvelfile;

  for (nr=0;nr<numreplica;++nr) {
    sprintf(dirdata,"%s/rep%d/data",dirbase,nr+1);
    mkdir(dirdata);
    sprintf(gtrjfilename,"%s/rep%d/data/%s_sam.gtrj",dirbase,nr+1,pn);
    sprintf(gvelfilename,"%s/rep%d/data/%s_sam.gvel",dirbase,nr+1,pn);
    gtrjfile=efopen(gtrjfilename,"a");
    gvelfile=efopen(gvelfilename,"a");
    for (i=0;i<numexchange;++i) {
      sprintf(trjfilename,"%s/rep%d/sam/%s_sam_%d.trj",dirbase,nr+1,pn,i+1);
      sprintf(velfilename,"%s/rep%d/sam/%s_sam_%d.vel",dirbase,nr+1,pn,i+1);
      trjfile=efopen(trjfilename,"r");
      velfile=efopen(velfilename,"r");
      trj=(double *)gcemalloc(sizeof(double)*numstep*numatom*3);
      vel=(double *)gcemalloc(sizeof(double)*numstep*numatom*3);
      io_inputtrj_for_gather_Amberform(trjfile,trj,numstep,numatom);
      io_inputtrj_for_gather_Amberform(velfile,vel,numstep,numatom);
      io_addtrj_Amberform(gtrjfile,trj,numstep,numatom);
      io_addtrj_Amberform(gvelfile,vel,numstep,numatom);
      fclose(trjfile);
      fclose(velfile);
    }
    fclose(gtrjfile);
    fclose(gvelfile);
  }
}
