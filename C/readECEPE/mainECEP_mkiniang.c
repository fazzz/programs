#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ECEPE.h"
#include "PROTOPO.h"
#include "PT.h"
#include "FF.h"
#include "TOPO.h"

#include "EF.h"

#define ON  1
#define OFF 0

#define MAXATOM 1000

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,m,n;
  int amberrstflag=OFF;
  int crdflag=OFF;
  int lastcrdflag=OFF;
  int outdflag=OFF;
  int deltadihedflag=OFF;
  int preoflag=OFF;
  int woddangflag=OFF;
  int nums;
  int d,numatom,inistep=0,numstep=-1;
  int num;
  int argnum;
  int *numdiofres;
  double fd;

  double theta,*old_theta;
  double crd_nc[MAXATOM][3],*crd;

  double *co;
  double *dihed,*delta_dihed,*dihedave,*dihed_dummy,*dihed_dummy2;
  int *index_dihed;
  double pi;

  int *bp;
  int **bp_f;
  int *numb;
  double atom[4][3];

  struct ECEPE_parms ECEPE_p;
  struct pnb nb_p;
  struct ECEPE_pote p;
  struct ECEPE_force f;

  char *name_atom_list,*name_res_list;

  char *progname;
  char *preofilename,*bd8filename,*coofilename;
  char  *angfilename,*trjfilename,*outputfilename,*outdihedfilename;
  char *deltadihedfilename;
  FILE *parmtopfile,*coofile,*angfile,*trjfile,*outputfile,*outdihedfile,*deltadihedfile;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double *dihed_dummy_dummy;
  int ntotaldih,ndihinres,numres;

  double *co_dummy;

  int PROflag=OFF;
  int *PROflaglist;

  progname=argv[0];
  while((c=getopt(argc,argv,"halocwpn:i:d:"))!=-1) {
    switch(c) {
    case 'a':
      amberrstflag=ON;
      break;
    case 'c':
      crdflag=ON;
      break;
    case 'p':
      preoflag=ON;
      break;
    case 'w':
      woddangflag=ON;
      break;
    case 'i':
      inistep=atoi(optarg);
      break;
    case 'n':
      numstep=atoi(optarg);
      break;
    case 'l':
      lastcrdflag=ON;
      break;
    case 'o':
      outdflag=ON;
      break;
    case 'd':
      deltadihedflag=ON;
      deltadihedfilename=optarg;
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  argc-=optind;
  argv+=optind;

  if (outdflag==ON) argnum=7;
  else argnum=6;
  if (argc < argnum) {
    USAGE(progname);
    exit(1);
  }
  preofilename   = *argv;
  bd8filename    = *++argv;
  coofilename    = *++argv;
  angfilename    = *++argv;
  trjfilename    = *++argv;
  outputfilename = *++argv;
  if (outdflag==ON) outdihedfilename = *++argv;

  pi=acos(-1.0);

  read_ECEPE_parm_wtransindex(preofilename,bd8filename,&ECEPE_p,&nb_p);

  PROflaglist=(int *)gcemalloc(sizeof(int)*ECEPE_p.NUMRES); //082211
  numdiofres=(int *)gcemalloc(sizeof(int)*ECEPE_p.NUMRES);
  for (i=0;i<ECEPE_p.NUMRES;++i){
    numdiofres[i]=0;
    PROflaglist[i]=OFF; //082211
  }
  num=1;
  for (i=0;i<ECEPE_p.NUMVAR;++i) {
    if (num==ECEPE_p.dihed[i].indexv1) {
      numdiofres[num-1]+=1;
    }
    else {
      ++num;
      numdiofres[num-1]+=1;
      if (ECEPE_p.dihed[i].indexv2==2)
	PROflaglist[num-1]=ON;
    }
  }

  bp=(int *)gcemalloc(sizeof(int)*(ECEPE_p.NUMATM-1)*2);
  bp_f=(int **)gcemalloc(sizeof(int *)*ECEPE_p.NUMATM);
  numb=(int *)gcemalloc(sizeof(int)*ECEPE_p.NUMATM);
  //  make_bd_pair_list(ECEPE_p,bp,bp_f,numb);
  name_atom_list=(char *)gcemalloc(sizeof(char)*ECEPE_p.NUMATM*4);
  name_res_list=(char *)gcemalloc(sizeof(char)*ECEPE_p.NUMATM*3);
  for (i=0;i<ECEPE_p.NUMATM;++i) {
    for (j=0;j<4;++j)
      name_atom_list[i*4+j]=ECEPE_p.atom[i].name_atom[j];
    for (j=0;j<3;++j)
      name_res_list[i*3+j]=ECEPE_p.atom[i].name_res[j];
  }
  make_bp_v2(ECEPE_p.NUMATM,bp,bp_f,numb,name_atom_list,name_atom_list,ON,ON);

  for (i=0;i<ECEPE_p.NUMATM;++i) {
    for (j=0;j<4;++j) printf("%c",name_atom_list[i*4+j]);
    printf("(%3d)-",i);
    for (j=0;j<numb[i];++j) {
      n=bp_f[i][j];
      for (k=0;k<4;++k) printf("%c",name_atom_list[n*4+k]);
      printf("(%3d) ",n);
    }
    printf("\n");
  }

  make_dihed_pairs_list_v2(ECEPE_p,bp_f,numb);
  //  chang_index_dihed_pairs(ECEPE_p.NUMVAR,ECEPE_p);

  co=(double *)gcemalloc(sizeof(double)*(ECEPE_p.NUMATM)*3);
  dihed=(double *)gcemalloc(sizeof(double)*(ECEPE_p.NUMVAR));
  dihedave=(double *)gcemalloc(sizeof(double)*(ECEPE_p.NUMVAR));
  delta_dihed=(double *)gcemalloc(sizeof(double)*(ECEPE_p.NUMVAR));
  dihed_dummy=(double *)gcemalloc(sizeof(double)*(ECEPE_p.NUMVAR));
  for (i=0;i<ECEPE_p.NUMVAR;++i) delta_dihed[i]=0.0;
  for (i=0;i<ECEPE_p.NUMVAR;++i) dihed_dummy[i]=0.0;

  if (deltadihedflag==OFF && woddangflag==OFF) {
    angfile=efopen(angfilename,"r");
    /*************************************************************************/
    /* fscanf(angfile,"%lf",&dihed_dummy[0]);				   */
    /* for (i=0;i<9;++i) fscanf(angfile,"%lf",&fd);			   */
    /* for (i=1;i<ECEPE_p.NUMVAR;++i) fscanf(angfile,"%lf",&dihed_dummy[i]); */
    /*************************************************************************/
    k=0;
    for (i=0;i<ECEPE_p.NUMRES;++i) {
      if (PROflaglist[i]==ON) {     // 082211
	fscanf(angfile,"%lf",&fd);  // 082211
	num=9;                      // 082211
      }                             // 082211
      else num=10;                  // 082211
      for (j=0;j<numdiofres[i];++j) {
	fscanf(angfile,"%lf",&dihed_dummy[k]);
	++k;
      }
      for (j=0;j<num/*10 082211*/-numdiofres[i];++j) {
	fscanf(angfile,"%lf",&fd);
      }
    }
    fclose(angfile);
  }
  ntotaldih=0;
  ndihinres=1;
  numres=1;
  dihed_dummy_dummy=(double *)gcemalloc(sizeof(double)*ECEPE_p.NUMVAR);
  index_dihed=(int *)gcemalloc(sizeof(int)*ECEPE_p.NUMVAR);
  for (i=0;i<ECEPE_p.NUMVAR;++i) {
    if (ECEPE_p.dihed[i].indexv1>numres) {
      numres=ECEPE_p.dihed[i].indexv1;
      ntotaldih+=ndihinres;
      ndihinres=1;
      if (ECEPE_p.dihed[i].indexv2==2) {  // 082211
	PROflag=ON;                        // 082211
      }                                    // 082211
      else {                               // 082211
	PROflag=OFF;                       // 082211
      }                                    // 082211
    }
    if (ECEPE_p.dihed[i].indexv2>ndihinres) {
      ndihinres=ECEPE_p.dihed[i].indexv2;
      if (PROflag==ON) ndihinres=2;        // 082211
    }
    if (PROflag==ON) {                                                          // 082211
      dihed_dummy_dummy[i]=dihed_dummy[ntotaldih+ECEPE_p.dihed[i].indexv2-1-1]; // 082211
      index_dihed[i]=ntotaldih+ECEPE_p.dihed[i].indexv2-1-1;                    // 082211
    }                                                                           // 082211
    else {                                                                      // 082211
      dihed_dummy_dummy[i]=dihed_dummy[ntotaldih+ECEPE_p.dihed[i].indexv2-1];
      index_dihed[i]=ntotaldih+ECEPE_p.dihed[i].indexv2-1;
    }                                                                           // 082211
  }
  for (i=0;i<ECEPE_p.NUMVAR;++i) {
    dihed_dummy[i]=dihed_dummy_dummy[i];
  }

  for (i=0;i<ECEPE_p.NUMVAR;++i) {
    dihed_dummy[i]=dihed_dummy[i]/180.0*pi;
    if (dihed_dummy[i]<-pi) dihed_dummy[i]+=2.0*pi;
    else if (dihed_dummy[i]>pi) dihed_dummy[i]-=2.0*pi;
  }

  if (deltadihedflag==OFF && woddangflag==OFF) {
    coofile=efopen(coofilename,"r");
    read_ECEPE_coo(coofile,co,dihed,ECEPE_p.NUMATM);
    fclose(coofile);
    co_dummy=(double *)gcemalloc(sizeof(double)*ECEPE_p.NUMATM*3);
    for (i=0;i<ECEPE_p.NUMATM;++i) 
      for (j=0;j<3;++j) 
	co_dummy[(ECEPE_p.atom[i].katom-1)*3+j]=co[i*3+j];
    //      co_dummy[i*3+j]=co[(ECEPE_p.atom[i].katom-1)*3+j];
    for (i=0;i<ECEPE_p.NUMATM;++i) 
      for (j=0;j<3;++j) 
	co[i*3+j]=co_dummy[i*3+j];
    
    calc_TORS_for_get_sabun(ECEPE_p.NUMVAR,co,ECEPE_p,dihed_dummy,delta_dihed);
  }

  if (deltadihedflag==ON) {
    deltadihedfile=efopen(deltadihedfilename,"r");
    get_sabun_from_delta_dihedfile(deltadihedfile,ECEPE_p.NUMVAR,ECEPE_p,delta_dihed);
    fclose(deltadihedfile);
  }

  trjfile=efopen(trjfilename,"r");

  for (i=0;i<ECEPE_p.NUMVAR;++i) dihedave[i]=0.0;

  old_theta=(double  *)gcemalloc(sizeof(double )*ECEPE_p.NUMVAR);
  for (i=0;i<ECEPE_p.NUMVAR;++i) old_theta[i]=0.0;

  crd=(double *)gcemalloc(sizeof(double)*(ECEPE_p.NUMATM)*3);
  //  trjfile=efopen(trjfilename,"r");
  if (outdflag==ON) outdihedfile=efopen(outdihedfilename,"w");
  numatom=ECEPE_p.NUMATM;
  nums=/*0*/-1;
  d=0;
  while (d!=-1 ) {
    if (amberrstflag==OFF && crdflag==OFF && preoflag==OFF)
      d=io_scanconfwj(trjfile,numatom,crd,'x');
    if (amberrstflag==ON && crdflag==OFF && preoflag==OFF)
      d=io_scanconf_Amber_rst(trjfile,crd);
    if (crdflag==ON && amberrstflag==OFF && preoflag==OFF) {
      fscanf(trjfile,"%d",&numatom);
      for (j=0;j<numatom;++j)
	for (k=0;k<3;++k)
	  fscanf(trjfile,"%lf",&crd[j*3+k]);
      d=-1;
    }      
    if (amberrstflag==OFF && crdflag==OFF && preoflag==OFF) {
      for (j=0;j<ECEPE_p.NUMATM;++j)
	for (k=0;k<3;++k)
	  //	  crd_nc[j][k]=crd[(ECEPE_p.atom[j].katom-1)*3+k];
	  crd_nc[j][k]=crd[j*3+k]; //0811
    }
    else if (preoflag==ON) {
      for (j=0;j<ECEPE_p.NUMATM;++j)
	for (k=0;k<3;++k)
	  crd_nc[j][k]=ECEPE_p.atom[(ECEPE_p.atom[j].katom-1)].refcoord[k];
    }
    else {
      for (j=0;j<ECEPE_p.NUMATM;++j)
	for (k=0;k<3;++k)
	  crd_nc[j][k]=crd[j*3+k];
    }      
    ++nums;
    if (nums== numstep)
      break;
    
    for (k=0;k<ECEPE_p.NUMVAR;++k) {
      for (m=0;m<4;++m)
    	for (l=0;l<3;++l)
    	  atom[m][l]=crd_nc[(ECEPE_p.dihed[k].dpairs[m])][l];
    
      theta=dih(atom[0],atom[1],atom[2],atom[3]);
      if (theta > pi)  theta-=2.0*pi;
      else if (theta < -pi)  theta+=2.0*pi;
      theta+=delta_dihed[k];
      if (outdflag==ON) fprintf(outdihedfile,"%10.8lf ",theta*180/pi);
      //      theta+=delta_dihed[k];
      //      if (theta > pi)  theta-=2.0*pi;
      //      else if (theta < -pi)  theta+=2.0*pi;
      //      if (fabs(theta-old_theta[k]) < fabs((theta-2.0*pi)-old_theta[k])) {
      //    	if (fabs(theta-old_theta[k]) < fabs((theta+2.0*pi)-old_theta[k])) {
      //    	  ;
      //    	}
      //    	else {
      //    	  theta=theta+2.0*pi;
      //    	}
      //      }
      //      else {
      //    	if (fabs((theta-2.0*pi)-old_theta[k]) < fabs((theta+2.0*pi)-old_theta[k])) {
      //    	  theta=theta-2.0*pi;
      //    	}
      //    	else {
      //    	  theta=theta+2.0*pi;
      //    	}
      //      }
      if (lastcrdflag==OFF && nums >= inistep ) dihedave[k]+=theta;
      else dihedave[k]=theta;      
    }
    if (outdflag==ON) fprintf(outdihedfile,"\n");
    if (amberrstflag==ON)
      break;
  }
  if (outdflag==ON)
    fclose(outdihedfile);
  fclose(trjfile);

  if (lastcrdflag==OFF)
    for (i=0;i<ECEPE_p.NUMVAR;++i)
      dihedave[i]/=nums;
  
  dihed_dummy2=(double *)gcemalloc(sizeof(double)*ECEPE_p.NUMVAR);
  for (i=0;i<ECEPE_p.NUMVAR;++i) {
    dihed_dummy2[index_dihed[i]]=dihedave[i];
  }

  outputfile=efopen(outputfilename,"w");
  k=0;
  for (i=0;i<ECEPE_p.NUMRES;++i) {
      if (PROflaglist[i]==ON) {     // 082211
	fprintf(outputfile," -75.000");  // 082211
	num=9;                      // 082211
      }                             // 082211
      else num=10;                  // 082211
    for (j=0;j<numdiofres[i];++j) {
      fprintf(outputfile,"%8.3lf",dihed_dummy2[k]*180.0/pi);
      ++k;
    }
    for (j=0;j<num/*10 082211*/-numdiofres[i];++j) {
      fprintf(outputfile,"%8.3lf",0.0);
    }
    fprintf(outputfile,"\n");
  }
  //  fprintf(outputfile,"%8.3lf",dihedave[0]*180/pi);
  //  for (i=0;i<9;++i) fprintf(outputfile,"%8.3lf",0.0);
  //  fprintf(outputfile,"\n");
  //  for (i=1;i<ECEPE_p.NUMVAR;++i) fprintf(outputfile,"%8.3lf",dihedave[i]*180/pi);
  //  for (i=0;i<10-ECEPE_p.NUMVAR+1;++i) fprintf(outputfile,"%8.3lf",0.0);
  //  fprintf(outputfile,"\n");
  //  for (i=0;i<2;++i) fprintf(outputfile,"%8.3lf",0.0);
  //  fprintf(outputfile,"%8.3lf",180.0);
  //  for (i=0;i<7;++i) fprintf(outputfile,"%8.3lf",0.0);
  //  fprintf(outputfile,"\n");
  fclose(outputfile);

  return 0;
}

void USAGE(char *progname) {
  printf("-a -- amberrstmode \n");
  printf("-c -- crdflag \n");
  printf("-p -- preoflag \n");
  printf("-n -- numstep \n");
  printf("-i -- inistep \n");
  printf("-l -- lastcrdflag \n");
  printf("-o -- outdihedangflag \n");
  printf("-d -- deltadihedflag \n");
  printf("-h -- help\n");
  printf("USAGE: %s profilename bd8filename coofilename angfilename trjfilename outputfilename (outdihedfilename)\n", progname);
}

