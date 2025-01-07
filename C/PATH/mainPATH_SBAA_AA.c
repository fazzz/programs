
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PATH.h"

int main(int argc, char *argv[]) {
  int i,j;
  int numpath;

  double beta,delta;
  int ns,interval;

  char *crd,*top,*ref1,*ref2,*inpf;
  char *ene,*dihed,*rmsdf1,*rmsdf2;
  char *inipath,*modpathbase,*crdforsam,*outbase,*pmf;
  char *dirbase;
  char *dirinipath,*dirinipathanal,*dirmodpath;
  char *dirsam,*dirsamanl,*dirsamcene,*dirsampmf;

  nodes=atoi(getenv("node"));
  node=getenv("node");
  dirbase=getenv("dirbase");
  //////////////////////////////////////////////////
  beta=atof(getenv("beta"));
  delta=atof(getenv("delta"));
  numstep=atoi(getenv("numstep"));
  interval=atof(getenv("interval"));
  //////////////////////////////////////////////////
  ref1=getenv("ref1");
  ref2=getenv("ref2");
  crd=getenv("crd");
  top=getenv("top");
  //////////////////////////////////////////////////
  dirinipath=getenv("dirinipath");
  inipathnc=getenv("inipathnc");
  //////////////////////////////////////////////////
  dirinipathanal=getenv("inipath");
  ene=getenv("ene");
  dihed=getenv("dihed");
  rmsdf1=getenv("rmsdf1");
  rmsdf2=getenv("rmsdf2");
  //////////////////////////////////////////////////
  dirmodpath=getenv("dirmodpath");
  inipathbase=getenv("inipathbase");
  criteria=atof(getenv("criteria"));
  condforginipath=getenv("condforginipath");
  modpathbase=getenv("modpathbase");
  numpoint=atoi(getenv("numpoint"));
  //////////////////////////////////////////////////
  betaforsam=getenv("betaforsam");
  numstepforsam=atoi(getenv("numstepforsam"));
  intervalforsam=atof(getenv("intervalforsam"));
  fconstforsam=atof(getenv("fconstforsam"));
  crdforsambase=getenv("crdforsambase");
  drestforsambase=atof(getenv("drestforsambase"));
  dirsam=getenv("dirsam");
  crdforsam=getenv("crdforsam");
  minrstforsam=getenv("minrstforsam");
  topforsam=getenv("topforsam");
  dirsam=getenv("dirsam");
  dirsaminp=getenv("dirsaminp");
  outbase=getenv("outbase");
  //////////////////////////////////////////////////
  dirsamanl=getenv("dirsamanl");
  dirsamcene=getenv("dirsamcene");
  dirsamfene=getenv("dirsamfene");
  dirsampmf=getenv("dirsampmf");
  samnc=getenv("samnc");
  eneforsam=getenv("samnc");
  dihedforsam=getenv("dihedforsam");
  //////////////////////////////////////////////////
  cene=getenv("cene");
  fene=getenv("fene");
  pmf=getenv("pmf");
  //////////////////////////////////////////////////

  execute_ini_path_wSBAAMC(beta,delta,
			   numstep,interval,
			   ref1,ref2,crd,top,
			   dirbase,dirinipath,inipathnc,
			   node,"INPH");
  waitjob("INPH");
  anl_ini_path_SBAA(dirbase,dirinipath,inipathnc,
		    dirinipathanal,
		    ene,dihed,rmsdf1,rmsdf2,ref1,ref2,
		    node,"AINPH");
  waitjob("AINPH");
  numpath=get_PATH_fSBAA(dirinipath,inipathnc,
			 criteria,condforginipath,
			 dirmodpath,inipathbase,
			 node,"GINPH");
  waitjob("GINPH");
  modify_ini_path(dirmodpath,inipathbase,
		  numpath,numpoint,
		  dirmodpath,modpathbase,
		  node,"MODPH");
  waitjob("MODPH");
  make_pdbfil_from_mpath(dirmodpath,modpathbase,
			 numpath,numpoint,
			 dirsaminp,pdbforsam,
			 node,"MKPDBSAM");
  waitjob("MKPDBSAM");
  make_cdfil_from_mpath(dirmodpath,modpathbase,
			numpath,numpoint,
			dirsaminp,crdforsam,
			node,"MKCDSAM");
  waitjob("MKCDSAM");
  execute_min_al_mpath_wAA(dirsaminp,crdforsam,
			   numpath,numpoint,
			   dirsaminp,minrstforsam,
			   node,"MINCDSAM");
  waitjob("MINCDSAM");
  execute_sam_al_mpath_wAA(betaforsam,
			   numstepforsam,intervalforsam,
			   minrstforsam,topforsam,
			   numpath,numpoint,
			   fconstforsam,drestforsambase,
			   dirsaminp,dirsam,samnc,
			   node,"SAMAMPATH");
  waitjob("SAMAMPATH");
  anl_sam_al_mpath_AA(dirsam,samnc,
		      numpath,numpoint,
		      dirsamanl,eneforsam,dihedforsam,
		      node,"ASAMAMP");
  waitjob("ASAMAMP");
  execute_cene(dirsam,samnc,
	       numpath,numpoint,
	       dirsamcene,cene,"CESAMAMP");
  waitjob("CESAMAMP");
  execute_MBAR(dirsam,samnc,
	       numpath,numpoint,
	       dirsamcene,cene,
	       dirsamanl,dihedforsam,
	       dirsamfene,fene,
	       dirsampmf,pmf,
	       node,"PMFSAMAM");

  return 0;
}
