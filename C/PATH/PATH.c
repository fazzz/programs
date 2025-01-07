
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <netcdf.h>

#include "PATH.h"
#include "TOPO.h"
#include "PT.h"
#include "rmsd.h"
#include "netcdf_mine.h"

#define INSIDE_ref1 0
#define INSIDE_ref2 1
#define OUT 2

#define INref1 0
#define INref2 1
#define PATH_1t2 2
#define PATH_2t1 3

int get_path_ftrj(char *trjname,
		  char *pathnamebase,
		  double *crdref1, 
		  double *crdref2,
		  double criteria, 
		  double *path, 
		  int *numpoints,
		  int MODE, 
		  int interval) {
  int i,j,k,l;
  int INOUT,state;
  int numatom,numstep,numpoint;
  int numpath;
  double *crd;
  double rmsd_from1,rmsd_from2;
  char pathname[100];

  FILE *pathfile;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_SBAAMCD nc_id_MCD;

  numpath=0;
  numstep=mync_get_present_step_SBAAMCD(trjname,&nc_id_MCD);
  numatom=AP.NATOM;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  state=-1;
  numpoint=0;
  for (i=0;i<numstep;++i) {
    if (i%interval==0) {
      mync_open_inq_get_sh_MCD(trjname,numatom,i,1,i+1,&nc_id_MCD,crd_nc);
      for (j=0;j<numatom;++j)  for (k=0;k<3;++k) crd[j*3+k]=crd_nc[j][k];

      rmsd_from1=rmsd_qcp(crd,crdref1,numatom,MODE);
      rmsd_from2=rmsd_qcp(crd,crdref2,numatom,MODE);

      if (rmsd_from1<criteria) INOUT=INSIDE_ref1;
      else if (rmsd_from2<criteria) INOUT=INSIDE_ref2;
      else INOUT=OUT;

      if (state==-1) {
	if (INOUT==INSIDE_ref1) state=INref1;
	else if (INOUT==INSIDE_ref2) state=INref2;
	else {
	  printf("error\n");
	  exit(1);
	}
      }
      else if (state==INref1) {
	if (INOUT==INSIDE_ref2) {
	  printf("error\n");
	  exit(1);
	}
	else if (INOUT== OUT) {
	  ++numpoint;
	  path=(double *)gcerealloc(path,sizeof(double)*numpoint*numatom*3);
	  for (j=0;j<numatom*3;++j) path[(numpoint-1)*numatom*3+j]=crd[j];
	  state=PATH_1t2;
	}
      }
      else if (state==INref2) {
	if (INOUT==INSIDE_ref1) {
	  printf("error\n");
	  exit(1);
	} 
	else if (INOUT==OUT) {
	  ++numpoint;
	  path=(double *)gcerealloc(path,sizeof(double)*numpoint*numatom*3);
	  for (j=0;j<numatom*3;++j) path[(numpoint-1)*numatom*3+j]=crd[j];
	  state=PATH_2t1;
	}
      }
      else if (state==PATH_1t2) {
	if (INOUT==INSIDE_ref1) {
	  numpoint=0;
	  path=(double *)gcerealloc(path,sizeof(double)*numatom*3);
	  state=INSIDE_ref1;
	}
	else if (INOUT==INSIDE_ref2) {
	  ++numpoint;
	  path=(double *)gcerealloc(path,sizeof(double)*numpoint*numatom*3);
	  for (j=0;j<numatom*3;++j) path[(numpoint-1)*numatom*3+j]=crd[j];
	  sprintf(pathname,"%s_path=%d",pathnamebase,numpath);
	  pathfile=efopen(pathname,"w");
	  for (j=0;j<numpoint;++j) { 
	    for (k=0;k<numatom;++k) { 
	      for (l=0;l<3;++l)
		fprintf(pathfile,"%e ",path[j*numatom*3+k*3+l]);
	      fprintf(pathfile,"\n");
	    }
	  }
	  state=INSIDE_ref2;
	  numpath++;
	  path=(double *)gcerealloc(path,sizeof(double)*numatom*3);
	  for (j=0;j<numatom*3;++j) path[j]=crd[j];
	}
      }
      else if (INOUT==OUT) {
	++numpoint;
	path=(double *)gcerealloc(path,sizeof(double)*numpoint*numatom*3);
	for (j=0;j<numatom*3;++j) path[(numpoint-1)*numatom*3+j]=crd[j];
      }
      else if (state==PATH_2t1) {
	if (INOUT==INSIDE_ref1) {
	  ++numpoint;
	  path=(double *)gcerealloc(path,sizeof(double)*numpoint*numatom*3);
	  for (j=0;j<numatom*3;++j) path[(numpoint-1)*numatom*3+j]=crd[j];
	  sprintf(pathname,"%s_path=%d",pathnamebase,numpath);
	  pathfile=efopen(pathname,"w");
	  for (j=0;j<numpoint;++j) { 
	    for (k=0;k<numatom;++k) { 
	      for (l=0;l<3;++l)
		fprintf(pathfile,"%e ",path[j*numatom*3+k*3+l]);
	      fprintf(pathfile,"\n");
	    }
	  }
	  state=INSIDE_ref1;
	  numpath++;
	  path=(double *)gcerealloc(path,sizeof(double)*numatom*3);
	  for (j=0;j<numatom*3;++j) path[j]=crd[j];
	} 
	else if (INOUT==INSIDE_ref2) {
	  numpoint=0;
	  path=(double *)gcerealloc(path,sizeof(double)*numatom*3);
	  for (j=0;j<numatom*3;++j) path[j]=crd[j];
	  state=INSIDE_ref2;
	}
	else if (INOUT==OUT) {
	  ++numpoint;
	  path=(double *)gcerealloc(path,sizeof(double)*numpoint*numatom*3);
	  for (j=0;j<numatom*3;++j) path[(numpoint-1)*numatom*3+j]=crd[j];
	}
      }
    }
  }
  
  return numpath+1;
}

int get_path_ftrj_bydha(double *trjname,
			char *pathnamebase,
			double *crdref1, 
			double *crdref2,
			double criteria,
			double *dha_spe,
			int *atomp,
			int numdha,
			double *path, 
			int *numpoints,
			int MODE, 
			int interval) {
  int i,j,k,l;
  int INOUT,state;
  int numatom,numstep,numpoint;
  int numpath;
  double *crd;
  double *dha;
  double pi,dang;
  double atom[4][3];
  char pathname[100];

  FILE *pathfile;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_SBAAMCD nc_id_MCD;

  pi=acos(-1.0);

  numpath=0;
  numstep=mync_get_present_step_SBAAMCD(trjname,&nc_id_MCD);
  numatom=AP.NATOM;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  dha=(double *)gcemalloc(sizeof(double)*numdha);

  state=-1;
  numpoint=0;
  for (i=0;i<numstep;++i) {
    if (i%interval==0) {
      mync_open_inq_get_sh_MCD(trjname,numatom,i,1,i+1,&nc_id_MCD,crd_nc);
      for (j=0;j<numatom;++j)  for (k=0;k<3;++k) crd[j*3+k]=crd_nc[j][k];

      for (j=0;j<numdha;++j) {
	for (k=0;k<4;++k) {
	  for (l=0;l<3;++l) {
	    atom[k][l]=crd[atomp[j*4+k]*3+l];
	  }
	  dha[j]=dih(atom[0],atom[1],atom[2],atom[3]);
	  if (dha[j] > pi) dha[j]-=2.0*pi;
	  else if (dha[j] <  -1.0*pi) dha[j]+=2.0*pi;
	}
      }

      INOUT=-1;
      for (j=0;j<numdha;++j) {
	dang=dha[j]-dha_spe[j*2];
	if (dang > pi) dang-=2.0*pi;
	else if (dang <  -1.0*pi) dang+=2.0*pi;
	if (fabs(dang)>criteria) {
	  INOUT=-1;
	  break;
	}
	else {
	  INOUT=INSIDE_ref1;
	}
      }
      if (INOUT==-1) {
	for (j=0;j<numdha;++j) {
	  dang=dha[j]-dha_spe[j*2+1];
	  if (dang > pi) dang-=2.0*pi;
	  else if (dang <  -1.0*pi) dang+=2.0*pi;
	  if (fabs(dang)>criteria) {
	    INOUT=-1;
	    break;
	  }
	  else {
	    INOUT=INSIDE_ref2;
	  }
	}
      }
      if (INOUT==-1) INOUT=OUT;

      if (state==-1) {
	if (INOUT==INSIDE_ref1) state=INref1;
	else if (INOUT==INSIDE_ref2) state=INref2;
	else {
	  printf("error\n");
	  exit(1);
	}
      }
      else if (state==INref1) {
	if (INOUT==INSIDE_ref2) {
	  printf("error\n");
	  exit(1);
	}
	else if (INOUT== OUT) {
	  ++numpoint;
	  path=(double *)gcerealloc(path,sizeof(double)*numpoint*numatom*3);
	  for (j=0;j<numatom*3;++j) path[(numpoint-1)*numatom*3+j]=crd[j];
	  state=PATH_1t2;
	}
      }
      else if (state==INref2) {
	if (INOUT==INSIDE_ref1) {
	  printf("error\n");
	  exit(1);
	} 
	else if (INOUT==OUT) {
	  ++numpoint;
	  path=(double *)gcerealloc(path,sizeof(double)*numpoint*numatom*3);
	  for (j=0;j<numatom*3;++j) path[(numpoint-1)*numatom*3+j]=crd[j];
	  state=PATH_2t1;
	}
      }
      else if (state==PATH_1t2) {
	if (INOUT==INSIDE_ref1) {
	  numpoint=0;
	  path=(double *)gcerealloc(path,sizeof(double)*numatom*3);
	  state=INSIDE_ref1;
	}
	else if (INOUT==INSIDE_ref2) {
	  ++numpoint;
	  path=(double *)gcerealloc(path,sizeof(double)*numpoint*numatom*3);
	  for (j=0;j<numatom*3;++j) path[(numpoint-1)*numatom*3+j]=crd[j];
	  sprintf(pathname,"%s_path=%d",pathnamebase,numpath);
	  pathfile=efopen(pathname,"w");
	  for (j=0;j<numpoint;++j) { 
	    for (k=0;k<numatom;++k) { 
	      for (l=0;l<3;++l)
		fprintf(pathfile,"%e ",path[j*numatom*3+k*3+l]);
	      fprintf(pathfile,"\n");
	    }
	  }
	  state=INSIDE_ref2;
	  numpath++;
	  path=(double *)gcerealloc(path,sizeof(double)*numatom*3);
	  for (j=0;j<numatom*3;++j) path[j]=crd[j];
	}
      }
      else if (INOUT==OUT) {
	++numpoint;
	path=(double *)gcerealloc(path,sizeof(double)*numpoint*numatom*3);
	for (j=0;j<numatom*3;++j) path[(numpoint-1)*numatom*3+j]=crd[j];
      }
      else if (state==PATH_2t1) {
	if (INOUT==INSIDE_ref1) {
	  ++numpoint;
	  path=(double *)gcerealloc(path,sizeof(double)*numpoint*numatom*3);
	  for (j=0;j<numatom*3;++j) path[(numpoint-1)*numatom*3+j]=crd[j];
	  sprintf(pathname,"%s_path=%d",pathnamebase,numpath);
	  pathfile=efopen(pathname,"w");
	  for (j=0;j<numpoint;++j) { 
	    for (k=0;k<numatom;++k) { 
	      for (l=0;l<3;++l)
		fprintf(pathfile,"%e ",path[j*numatom*3+k*3+l]);
	      fprintf(pathfile,"\n");
	    }
	  }
	  state=INSIDE_ref1;
	  numpath++;
	  path=(double *)gcerealloc(path,sizeof(double)*numatom*3);
	  for (j=0;j<numatom*3;++j) path[j]=crd[j];
	} 
	else if (INOUT==INSIDE_ref2) {
	  numpoint=0;
	  path=(double *)gcerealloc(path,sizeof(double)*numatom*3);
	  for (j=0;j<numatom*3;++j) path[j]=crd[j];
	  state=INSIDE_ref2;
	}
	else if (INOUT==OUT) {
	  ++numpoint;
	  path=(double *)gcerealloc(path,sizeof(double)*numpoint*numatom*3);
	  for (j=0;j<numatom*3;++j) path[(numpoint-1)*numatom*3+j]=crd[j];
	}
      }
    }
  }
  
  return numpath+1;
}

int execute_ini_path_wSBAAMC(double beta,double delta,
			     int ns,int interval,
			     char *ref1,char *ref2,char *crd,char *top,
			     char *dirbase,char *dirinipath,
			     char *inipathnc,
			     int numnode, char *name) {
  FILE *qsubsh;
  char qsubfilenamefullpath[100],dircommand[100],dirlog[100];

  mkdir(dirinipath);
  sprintf(dircommand,"%s/command",dirbase);
  sprintf(dirlog,"%s/log",dirbase);
  mkdir(dircommand);
  mkdir(dirlog);
  sprintf(qsubfilenamefullpath,"%s/qsub_%s.sh",dircommand,name);
  makeqsubbase(dircommand,dirlog,dir,name,numnode);  
  qsubsh=efopen(qsubfilenamefullpath,"a");
  fprintf(qsubsh,"~/mybin/MC_pep_SBAAMBFF_nc ");
  fprintf(qsubsh,"-b %e ",beta);
  fprintf(qsubsh,"-l %e ",delta);
  fprintf(qsubsh,"-n %d ",numstep);
  fprintf(qsubsh,"-i %d ",interval);
  fprintf(qsubsh,"-j %d ",interval);
  fprintf(qsubsh,"-k %d ",interval);
  fprintf(qsubsh,"-l %d ",interval);
  fprintf(qsubsh,"-d 0.001 ");
  fprintf(qsubsh,"%s %s %s %s ",ref1,ref2,crd,top);
  fprintf(qsubsh,"%s/%s \n",dirinipath,inipathnc);
  fclose(qsubsh);

  system("chmod +x %s",qsubfilenamefullpath);
  //  system(qsubfilenamefullpath);

  return 1;
}

int anl_ini_path_SBAA(char *dirbase,char *dirinipath,char *inipathnc,
		      char *dirinipathanal, 
		      char *ene, 
		      char *dihed ,
		      char *rmsdf1 , char *rmsdf2, 
		      char *ref1, char *ref2,
		      int numnode, char *name) {
  FILE *qsubsh;
  char qsubfilenamefullpath[100],dircommand[100],dirlog[100];

  sprintf(dircommand,"%s/command",dirbase);
  sprintf(dirlog,"%s/log",dirbase);
  sprintf(qsubfilenamefullpath,"%s/qsub_%s.sh",dircommand,name);
  makeqsubbase(dircommand,dirlog,dir,name,numnode);  
  qsubsh=efopen(qsubfilenamefullpath,"a");
  fprintf(qsubsh,"~/mybin/MC_get_ene_term_SBAA ");
  fprintf(qsubsh,"%s/%s %s/%s.ene \n",
	  dirinipath,inipathnc,
	  dirinipathanal,ene);
  fprintf(qsubsh,"~/mybin/CD_nc ");
  fprintf(qsubsh,"-P -p %s %s/%s %s/%s.dtrj \n",
	  top,dirinipath,inipathnc,dirinipathanal,dihed);
  fprintf(qsubsh,"~/mybin/rmsd_nc ");
  fprintf(qsubsh,"-c %s/%s %s %s %s/%s.rmsd_f1 \n",
	  dirinipath,inipathnc,ref1,top,dirinipathanal,rmsd_f1);
  fprintf(qsubsh,"~/mybin/rmsd_nc ");
  fprintf(qsubsh,"-c %s/%s %s %s %s/%s.rmsd_f2 \n",
	  dirinipath,inipathnc,ref2,top,dirinipathanal,rmsd_f2);
  fclose(qsubsh);

  system("chmod +x %s",qsubfilenamefullpath);
  //  system(qsubfilenamefullpath);

  return 1;
}

int get_PATH_fSBAA(char *dirinipath,char *inipathnc,double criteia,
		   char *cond, int numnode, char *name) {
  int numpath;
  FILE *qsubsh;
  char qsubfilenamefullpath[100],dircommand[100],dircond[100];

  sprintf(dircommand,"%s/command",dirbase);
  sprintf(dilog,"%s/log",dirbase);
  sprintf(qsubfilenamefullpath,"%s/qsub_%s.sh",dircommand,name);
  makeqsubbase(dircommand,dirlog,dir,name,numnode);
  qsubsh=efopen(qsubfilenamefullpath,"a");
  fprintf(qsubsh,"~/mybin/get_path_ftrj_bydha ");
  fprintf(qsubsh,"-c %e %s/%s %s/%s %s  %s/%s\n",
	  criteia,
	  dirinipathanal,inipathnc,
	  dircond,cond,top,
	  dirmodpath,inipathbase);
  fclose(qsubsh)
  system("chmod +x %s",qsubfilenamefullpath);
  //  system(qsubfilenamefullpath);

  return 0/*numpath*/;
}

int modify_ini_path(char *dirinp,char *inipathbase,
		    int numpath,int numpoint,
		    char *dirout,char *modpsthbase,
		    int numnode, char *name) {
  int i;
  FILE *qsubsh;
  char qsubfilenamefullpath[100],dircommand[100],dirlog[100],,inipath[100];
  char modpath[100],modpathene[100];

  //  for (i=0;i<numpath;++i) {
    sprintf(dircommand,"%s/command",dirbase);
    sprintf(dirlog,"%s/log",dirbase);
    sprintf(qsubfilenamefullpath,"%s/qsub_%s_%d.sh",dircommand,name,numpath);
    sprintf(inipath,"%s_%d",inipathbase,numpath);
    sprintf(modpath,"%s_%d",modpsthbase,numpath);
    sprintf(modpathene,"%s_%d_ene",modpsthbase,numpath);
    makeqsubbase(dircommand,dirlog.dir,name,numnode);
    qsubsh=efopen(qsubfilenamefullpath,"a");
    fprintf(qsubsh,"~/mybin/zstring_woH ");
    fprintf(qsubsh,"-t 0.0001 ");
    fprintf(qsubsh,"%d %s/%s %s %s/%s %s/%s\n",
	    numpoint,
	    dirini,inipathpath,
	    dircon,top,
	    dirout,modpath,dirmodpath,modpathene);
    fclose(qsubsh);
    system("chmod +x %s",qsubfilenamefullpath);
    //    system(qsubfilenamefullpath);
    //  }
}

int make_pdbfil_from_mpath(char *dirinp,char *modpsthbase,
			   int numpath, int numpoint,
			   char *dirout,char *pdbbase,
			   int numnode, char *name) {
  int i;
  FILE *qsubsh;
  char qsubfilenamefullpath[100],dircommand[100],dirlog[100],filinp[100];

  //  for (i=0;i<numpath;++i) {
    sprintf(dircommand,"%s/command",dirbase);
    sprintf(dirlog,"%s/log",dirbase);
    sprintf(qsubfilenamefullpath,"%s/qsub_%s_%d.sh",dircommand,name,numpath);
    makeqsubbase(dircommand,dirlog.dir,name,numnode);
    qsubsh=efopen(qsubfilenamefullpath,"a");
    fprintf(qsubsh,"~/mybin/get_pdbs_fpath -c ");
    sprintf(filinp,"%s/%s_%d",dirinp,modpsthbase,numpath);
    fprintf(qsubsh,"%s %s/%s.%d\n",top,dirout,pdbbase,numpath);
    fclose(qsubsh);
    system("chmod +x %s",qsubfilenamefullpath);
    //    system(qsubfilenamefullpath);
    //  }
}

int make_crdfil_from_mpath(char *dirmodpath,char *modpsthbase,
			   int numpath, int *numpoint,
			   char *dirsaminp,char *crdbase,
			   int *nodes,int numnode, char *name)  {
  int i,j;
  FILE *qsubsh,*cmd;
  char qsubfilenamefullpath[100],dircommand[100];
  char filinp[100],cmdname[100],crdname[100];

  //  for (i=0;i<numpath;++i) {
  for (j=0;j<numpoint[i];++j) {
    sprintf(cmdname,"%s/cmd_%d",dircommand,i);
    cmd=efopen(cmdname,"w");
    sprintf(filinp,"%s_%d.pdb",filinpbase,i+1);
    fprintf(qsubsh,"pn=loadPdb %s\nsaveAmberparm pn top .crd\nquit",filinp);
    fclose(cmd);
    sprintf(dircommand,"%s/command",dirbase);
    sprintf(qsubfilenamefullpath,"%s/qsub_mCDFINPH.sh",dircommand);
    makeqsubbase(dircommand,dir,"mCDFINPH",nodes,numnode,log,qsubfilename);  
    qsubsh=efopen(qsubfilenamefullpath,"a");
    fprintf(qsubsh,"/home/appl/amber10/exe/tleap ");
    fprintf(qsubsh,"-s -f //home/appl/amber10/dat/leap/cmd/oldff/leaprc.ff99 -f cmd\n");
    system("chmod +x %s",qsubfilenamefullpath);
    system(qsubfilenamefullpath);
  }
  //  }
}

int execute_sam_al_mpath_wAA(char *dirinp,char *filinp,
			     char *dirout,char *filout,
			     int *nodes,int numnode) {
  int i,j;
  FILE *qsubsh;
  char qsubfilenamefullpath[100],dircommand[100],filinp[100];

  //  for (i=0;i<numpath;++i) {
  for (j=0;j<numpath;++j) {
    sprintf(dircommand,"%s/command",dirbase);
    sprintf(qsubfilenamefullpath,"%s/qsub_MCSAFFAA.sh",dircommand);
    makeqsubbase(dircommand,dir,"MCSAFFAA",nodes,numnode,log,qsubfilename);  
    fprintf(qsubsh,"~/mybin/MC_pep_nc -a -b ");
    fprintf(qsubsh,"-b %e",beta);
    fprintf(qsubsh,"-n %d",numstep);
    fprintf(qsubsh,"-i %d -j %d -k %d",interval,interval,interval);
    fprintf(qsubsh,"-d 0.001 ");
    fprintf(qsubsh,"%s/%s ",dirinp,filinp);
    fprintf(qsubsh,"%s ",top);
    fprintf(qsubsh,"%s/%s ",dirout,filout);
    fprintf(qsubsh,"%s/%s ",dirinp,condfil);
    system("chmod +x %s",qsubfilenamefullpath);
    system(qsubfilenamefullpath);
  }
  //  }
}

int anl_sam_al_mpath_AA(char *dirinp,char *filinp,
			char *dirout,char *filout,
			int *nodes,int numnode) {
  int i,j;
  FILE *qsubsh;
  char qsubfilenamefullpath[100],dircommand[100],filinp[100];

  //  for (i=0;i<numpath;++i) {
  for (j=0;j<numpath;++j) {
    sprintf(dircommand,"%s/command",dirbase);
    sprintf(qsubfilenamefullpath,"%s/qsub_MCSAFFAA.sh",dircommand);
    makeqsubbase(dircommand,dir,"MCSAFFAA",nodes,numnode,log,qsubfilename);  
    fprintf(qsubsh,"~/mybin/MC_pep_nc -a -b ");
    fprintf(qsubsh,"-b %e",beta);
    fprintf(qsubsh,"-n %d",numstep);
    fprintf(qsubsh,"-i %d -j %d -k %d",interval,interval,interval);
    fprintf(qsubsh,"-d 0.001 ");
    fprintf(qsubsh,"%s/%s ",dirinp,filinp);
    fprintf(qsubsh,"%s ",top);
    fprintf(qsubsh,"%s/%s ",dirout,filout);
    fprintf(qsubsh,"%s/%s ",dirinp,condfil);
    system("chmod +x %s",qsubfilenamefullpath);
    system(qsubfilenamefullpath);
  }
  //  }
}

int execute_cene(char *dirinp,char *filinp,
		 char *dirout,char *filout,
		 int *nodes,int numnode) {

}

int execute_MBAR(char *dirinp,char *filinp,
		 char *dirout,char *filout,
		 int *nodes,int numnode) {

}
