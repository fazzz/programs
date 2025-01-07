
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <netcdf.h>
#include <getopt.h>

#include "EM.h"
#include "EMalg.h"
#include "EMalg_non_square_range.h"

#include "Simpson_integ.h"

#include "K_means.h"
#include "K_means_non_square_range.h"

#include "Gaussian.h"

#include "STRING_CV.h"
#include "CSI.h"
#include "EF.h"

#define ON 0
#define OFF 1

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,n;
  double fx,fx_old,fy_old,fy,fp;

  int N,K;
  int Nx,Ny;
  int Nx_old,Ny_old;

  int dim;
  int flag,floatflag,flagnum,periodflag=OFF,flagd=OFF,flagminmaxcheck[4];

  int PathSearchflag=OFF;

  double **prob,*prob2,sum;

  double **prob_app;

  double KL;

  double *x;
  double *ax;
  
  // parameters of Mixed-Gaussian
  int *num_k;
  double **nyu_k,***Sigma_k, *pi_k,***gamma_nk,**ave_k;
  double Sigma_initial=1.0,Sigma_initial2=0.0;
  // parameters of Mixed-Gaussian
  int *vec;
  double ***dist,J=0.0,J_new,dJ;
  double lnrou,lnrou_new,dlnrou,threhold=1.0e-4;

  double pi,x_map[2],minx, maxx, dx=0.01, miny, maxy, dy=0.01,minpmf;

  int Nx_app,Ny_app;
  double x_app[2];
  double dx_app=0.01,dy_app=0.01;

  double minxD, maxxD, minyD, maxyD;

  double xi,xf,yi,yf,dx2;
  double a,b;

  int numpoint,numiteration=10000,outinterval=1;

  int flagparmout=OFF;

  double dt=0.01;
  double maxdelta;
  double threhold_path=1.0e-5;

  double **path,**path_evoluted;
  double **fe,fe_dummy,*pe,pe_dummy,*pe_old;

  double *CV,*CV_dummy;

  char *inputfilename,*outputfilename,*pathfilename,*pmffilename,*KLdfilename,*parmfilename;
  FILE *inputfile,*outputfile,*pathfile,*pmffile,*KLdfile,*parmfile;
  
  char *line;
  size_t len=0;
  
  int c,d;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *progname;
  int opt_idx=1;
  
  pi=acos(-1.0);
  
  struct option long_opt[] = {
    {"outi",1,NULL,'o'},
    {"d",0,NULL,'d'},
    {"minx",1,NULL,'x'},
    {"maxx",1,NULL,'a'},
    {"miny",1,NULL,'y'},
    {"maxy",1,NULL,'m'},
    {"pfile",1,NULL,'f'},
    {"sigi",1,NULL,'i'},
    {"sigi2",1,NULL,'j'},
    {"dx",1,NULL,'1'},
    {"dy",1,NULL,'2'},
    {"dt",1,NULL,'t'},
    {"p",0,NULL,'p'},
    {"P",0,NULL,'P'},
    {"K",1,NULL,'K'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"hPdpK:o:x:a:y:m:f:i:j:1:2:t:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'K':
      K = atoi(optarg);  break;
    case 'P':
      PathSearchflag=ON;  break;
    case 'd':
      flagd=ON;  break;
    case 'x':
      minxD=atof(optarg);  break;
    case 'a':
      maxxD=atof(optarg);  break;
    case 'y':
      minyD=atof(optarg);  break;
    case 'm':
      maxyD=atof(optarg);  break;
    case 'o':
      outinterval=atoi(optarg);  break;
    case 'i':
      Sigma_initial=atof(optarg);  break;
    case 'j':
      Sigma_initial2=atof(optarg);  break;
    case '1':
      dx_app=atof(optarg);  break;
    case '2':
      dy_app=atof(optarg);  break;
    case 't':
      dt=atof(optarg);  break;
    case 'p':
      periodflag=ON;  break;
    case 'f':
      parmfilename=optarg; 
      flagparmout=ON;
      break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);  exit(1);
    }
  }

  //  numpoint=atoi(*++argv);
  //  xi=atof(*++argv);
  //  yi=atof(*++argv);
  //  xf=atof(*++argv);
  //  yf=atof(*++argv);
  //  pathfilename   = *++argv;
  //  KLdfilename    = *++argv;
  
  progname=*argv;  argc-=optind;  argv+=optind;
  
  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  numsim = atoi(*argv);
  inputfilelistname = *++argv;
  outputfilename    = *++argv;
  pmffilename       = *++argv;

  x_ij=(double **)gcemalloc(sizeof(double *)*(num_sim));
  for (i=0;i<num_sim;++i) x_ij[i]=(double *)gcemalloc(sizeof(double));

  inputfilelist=efopen(inputfilelistname,"r");
  for (i=0;i<num_sim;++i) {
    getline(&line,&len,enelistfile);
    l=strlen(line);
    line[l-1]='\0';

    inputfile=efopen(line,"r");
    num=0;
    d = 1;
    while ( d != -1  )  {
      d=fscanf(enefile,"%lf",&f);
      x_ij[i]=(double *)gcerealloc(x_ij[j],sizeof(double)*(num+1));
      x_ij[i][num]=f;
      ++num;
    }
    fclose(inputfile);
  }
  fclose(inputfilelist);

  prob=(double **)gcemalloc(sizeof(double *)*Nx);
  for (i=0;i<Nx;++i) {
    prob[i]=(double *)gcemalloc(sizeof(double)*Ny);
  }

  prob_app=(double **)gcemalloc(sizeof(double *)*Nx);
  for (i=0;i<Nx;++i) {
    prob_app[i]=(double *)gcemalloc(sizeof(double)*Ny);
  }

  dist=(double ***)gcemalloc(sizeof(double **)*Nx);
  for (i=0;i<Nx;++i) {
    //    dist[i]=(double **)gcemalloc(sizeof(double *)*Nx);
    dist[i]=(double **)gcemalloc(sizeof(double *)*Ny);
    //    for (j=0;j<Ny;++j) dist[i][j]=(double *)gcemalloc(sizeof(double)*Ny);
    for (j=0;j<Ny;++j) dist[i][j]=(double *)gcemalloc(sizeof(double)*K); // 2013-06-07
  }

  n=0;
  for (i=0;i<Nx;++i) {
    for (j=0;j<Ny;++j) {
      prob[i][j]=prob2[n];
      ++n;
    }
  }

  x=(double *)gcemalloc(sizeof(double)*2);

  sum=0.0;
  for (i=0;i<Nx;++i) {
    for (j=0;j<Ny;++j) {
      if (prob[i][j]!=0.0) {
	sum+=prob[i][j];
      }
    }
  }
  for (i=0;i<Nx;++i) {
    for (j=0;j<Ny;++j) {
      if (prob[i][j]!=0.0) {
	prob[i][j]=prob[i][j]/sum;
      }
    }
  }

  nyu_k=(double **)gcemalloc(sizeof(double *)*K);
  for (i=0;i<K;++i) nyu_k[i]=(double *)gcemalloc(sizeof(double)*2);

  Sigma_k=(double ***)gcemalloc(sizeof(double **)*K);
  for (i=0;i<K;++i) {
    Sigma_k[i]=(double **)gcemalloc(sizeof(double *)*2);
    for (j=0;j<2;++j) Sigma_k[i][j]=(double *)gcemalloc(sizeof(double )*2);
  }

  ave_k=(double **)gcemalloc(sizeof(double *)*K);
  for (i=0;i<K;++i) {
    ave_k[i]=(double *)gcemalloc(sizeof(double)*2);
  }

  pi_k=(double *)gcemalloc(sizeof(double )*K);

  gamma_nk=(double ***)gcemalloc(sizeof(double **)*Nx);
  for (i=0;i<Nx;++i) gamma_nk[i]=(double **)gcemalloc(sizeof(double *)*Ny);
  for (i=0;i<Nx;++i) 
    for (j=0;j<Ny;++j)
      gamma_nk[i][j]=(double *)gcemalloc(sizeof(double)*K);

  for (i=0;i<K;++i) {
    nyu_k[i][0]=minx+dx*i;
    nyu_k[K-i-1][1]=miny+dy*i;
  }

  //  x=(double *)gcemalloc(sizeof(double)*2);

  for (i=0;i<Nx;++i) {
    x[0]=minx+i*dx;
    for (j=0;j<Ny;++j) {
      x[1]=miny+i*dy;
      if (prob[i][j]!=0.0) {
	for (k=0;k<K;++k) {
	  dist[i][j][k]=(x[0]-nyu_k[k][0])*(x[0]-nyu_k[k][0])+(x[1]-nyu_k[k][1])*(x[1]-nyu_k[k][1]);
	}
      }
    }
  }
  Kmeans_Estep_fprob_nsr(Nx,minx,dx,Ny,miny,dy,
			 K,prob,nyu_k,gamma_nk,dist);
  J=Kmeans_J_fprob_nsr(Nx,minx,dx,Ny,miny,dy,
		       K,prob,nyu_k,gamma_nk,dist);

  i=0;

  printf("# ofiteration = %3d J = %8.3lf\n",i,J);
  for (j=0;j<K;++j) {
    printf("nyu (%3d) = ( ",j);
    for (k=0;k<2;++k) {
      printf("%8.3lf ",nyu_k[j][k]);
    }
    printf(")\n");
  }

  while (dJ > threhold || i==0) {

    Kmeans_Estep_fprob_nsr(Nx,minx,dx,Ny,miny,dy,
			   K,prob,nyu_k,gamma_nk,dist);

    Kmeans_Mstep_fprob_nsr(Nx,minx,dx,Ny,miny,dy,
			   K,prob,nyu_k,gamma_nk);

    J_new=Kmeans_J_fprob_nsr(Nx,minx,dx,Ny,miny,dy,
			     K,prob,nyu_k,gamma_nk,dist);

    dJ=fabs(J_new-J);
    J=J_new;

    ++i;
    printf("# ofiteration = %3d J = %8.3lf\n",i,J);
    for (j=0;j<K;++j) {
      printf("nyu (%3d) = ( ",j);
      for (k=0;k<2;++k) {
	printf("%8.3lf ",nyu_k[j][k]);
      }
      printf(")\n");
    }
  }

  vec=(int *)gcemalloc(sizeof(int)*Nx*Ny);

  n=0;
  for (i=0;i<Nx;++i) {
    for (j=0;j<Ny;++j) {
      for (k=0;k<K;++k) {
	if (gamma_nk[i][j][k]==1.0)
	  vec[n]=k;
      }
      ++n;
      //      }
    }
  }

  for (i=0;i<K;++i) {
    pi_k[i]=0.0;
  }

  n=0;
  for (i=0;i<Nx;++i) {
    for (j=0;j<Ny;++j) {
      pi_k[vec[n]]+=1.0/(Nx*Ny);
      ++n;
    }
  }

  num_k=(int *)gcemalloc(sizeof(int)*K);

  for (i=0;i<K;++i) {
    Sigma_k[i][0][0]=0.0;
    Sigma_k[i][1][1]=0.0;    
    ave_k[i][0]=0.0;
    ave_k[i][1]=0.0;    
  }

  for (i=0;i<Nx;++i) {
    x[0]=minx+i*dx;
    for (j=0;j<Ny;++j) {
      x[1]=miny+j*dy;
      for (k=0;k<K;++k) {
	Sigma_k[k][0][0]+=gamma_nk[i][j][k]*x[0]*x[0]*prob[i][j];
	Sigma_k[k][1][1]+=gamma_nk[i][j][k]*x[1]*x[1]*prob[i][j];    
	ave_k[k][0]+=gamma_nk[i][j][k]*x[0]*prob[i][j];
	ave_k[k][1]+=gamma_nk[i][j][k]*x[1]*prob[i][j];    
      }
    }
  }

  for (i=0;i<K;++i) {
    pi_k[i]=0.1/*1.0*//*10.0*/;
  }

  for (i=0;i<K;++i) {
    Sigma_k[i][0][0]-=ave_k[i][0]*ave_k[i][0];
    Sigma_k[i][1][1]-=ave_k[i][1]*ave_k[i][1];    
  }

  for (i=0;i<K;++i) {
    Sigma_k[i][0][0]=Sigma_initial;
    Sigma_k[i][1][0]=Sigma_initial2;
    Sigma_k[i][0][1]=Sigma_initial2;
    Sigma_k[i][1][1]=Sigma_initial;
  }

  lnrou=EM_lnrou_fprob_nsr(Nx,minx,dx,Ny,miny,dy,
			   K,prob,nyu_k,Sigma_k,pi_k,pi);
  i=0;
  dlnrou=0.0;
  printf("# ofiteration = %3d ln(rou) = %8.3lf\n",i,lnrou);
  for (j=0;j<K;++j) {
    printf("nyu (%3d) = ( ",j);
    for (k=0;k<2;++k) {
      printf("%8.3lf ",nyu_k[j][k]);
    }
    printf("Sigma (%3d) = ( ",j);
    for (k=0;k<2;++k) {
      printf("%10.6lf ",Sigma_k[j][k][k]);
    }
    printf(")\n");
  }

  while (dlnrou > threhold || i==0) {

    E_step_fprob_nsr(Nx,minx,dx,Ny,miny,dy,K,prob,nyu_k,Sigma_k,pi_k,gamma_nk,pi);

    //    M_step_fprob_nsr(Nx,minx,dx,Ny,miny,dy,K,prob,nyu_k,Sigma_k,pi_k,gamma_nk);
    M_step_fprob_nsr_2013_06_07(Nx,minx,dx,Ny,miny,dy,K,prob,nyu_k,Sigma_k,pi_k,gamma_nk); // 2013-06-07
    
    //lnrou_new=EM_lnrou_fprob_nsr(Nx,minx,dx,Ny,miny,dy,K,prob,nyu_k,Sigma_k,pi_k,pi);
    lnrou_new=EM_lnrou_fprob_nsr2(Nx,minx,dx,Ny,miny,dy,K,prob,nyu_k,Sigma_k,pi_k,pi);     // 2013-06-07

    dlnrou=fabs(lnrou_new-lnrou);
    lnrou=lnrou_new;

    ++i;
    printf("# ofiteration = %3d ln(rou) = %8.3lf\n",i,lnrou);
    for (j=0;j<K;++j) {
      printf("nyu (%3d) = ( ",j);
      for (k=0;k<2;++k) {
	printf("%8.3lf ",nyu_k[j][k]);
      }
      printf(") ");
      printf("Sigma (%3d) = ( ",j);
      for (k=0;k<2;++k) {
	printf("%10.6lf ",Sigma_k[j][k][k]);
      }
      printf(") ");
      printf("pi=%10.6lf \n",pi_k[j]);
    }
  }

  if (flagparmout==ON) {
    parmfile=efopen(parmfilename,"w");
    for (i=0;i<K;++i) {
      fprintf(parmfile,"%10.8e %10.8e ",nyu_k[i][0],nyu_k[i][1]);
      fprintf(parmfile,"%10.8e %10.8e ",Sigma_k[i][0],Sigma_k[i][1]);
      fprintf(parmfile,"%10.8e \n",pi_k[i]);
    }    
    fclose(parmfile);
  }

  prob2=Create_mixed_twoD_GaussianMap(minx, maxx, dx_app,
				      miny, maxy, dy_app,
				      nyu_k,  Sigma_k, pi_k, K,  pi);

  n=0;
  x_app[0]=minx;
  for (i=0;x_app[0]<maxx;++i){
    x_app[0]=minx+dx_app*i;
    x_app[1]=miny;
    for (j=0;x_app[1]<maxy;++j){
      x_app[1]=miny+dy_app*j;
      prob2[n]=-1.0*log(prob2[n]);
      ++n;
    }
  }

  n=0;
  minpmf=prob2[0];
  x_app[0]=minx;
  for (i=0;x_app[0]<maxx;++i){
    x_app[0]=minx+dx_app*i;
    x_app[1]=miny;
    for (j=0;x_app[1]<maxy;++j){
      x_app[1]=miny+dy_app*j;
      if (minpmf>prob2[n])
	minpmf=prob2[n];
      ++n;
    }
  }
  
  n=0;
  x_app[0]=minx;
  for (i=0;x_app[0]<maxx;++i){
    x_app[0]=minx+dx_app*i;
    x_app[1]=miny;
    for (j=0;x_app[1]<maxy;++j){
      x_app[1]=miny+dy_app*j;
      prob2[n]=prob2[n]-minpmf;
      ++n;
    }
  }

  pmffile=efopen(pmffilename,"w");
  n=0;
  x_map[0]=minx;
  for (i=0;x_map[0]<maxx;++i) {
    x_map[0]=minx+dx_app*i;
    x_map[1]=miny;
    for (j=0;x_map[1]<maxy;++j){
      x_map[1]=miny+dy_app*j;	
      fprintf(pmffile,"%8.3lf %8.3lf %8.3e\n",x_map[0],x_map[1],prob2[n]);
      ++n;
    }
  }
  fclose(pmffile);

  a=(yi-yf)/(xi-xf);
  b=yi-(yi-yf)/(xi-xf)*xi;

  dx2=(xf-xi)/(numpoint-1);

  path=(double **)gcemalloc(sizeof(double *)*numpoint);
  path_evoluted=(double **)gcemalloc(sizeof(double *)*numpoint);
  fe=(double **)gcemalloc(sizeof(double *)*numpoint);
  pe=(double *)gcemalloc(sizeof(double)*numpoint);
  pe_old=(double *)gcemalloc(sizeof(double)*numpoint);
  for (i=0;i<numpoint;++i) {
    path[i]=(double *)gcemalloc(sizeof(double)*2);
    path_evoluted[i]=(double *)gcemalloc(sizeof(double)*2);
    fe[i]=(double *)gcemalloc(sizeof(double)*2);
  }

  CV=(double *)gcemalloc(sizeof(double)*2);
  CV_dummy=(double *)gcemalloc(sizeof(double)*2);

  for (i=0;i<numpoint;++i) {
    path[i][0]=xi+dx2*i;
    path[i][1]=a*path[i][0]+b;

  }

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numpoint;++i) {
    for (j=0;j<2;++j) {
      fprintf(outputfile,"%8.3lf ",path[i][j]);
    }
    fprintf(outputfile,"\n");
  }
  fprintf(outputfile,"\n");

  for (i=0;i<numiteration;++i) {
    for (j=0;j<numpoint;++j) {
      pe[j]=0.0;
      for (k=0;k<2;++k) fe[j][k]=0.0;
    }

    for (j=0;j<numpoint;++j) {

      for (k=0;k<2;++k) {
	CV[k]=path[j][k];
      }

      if (periodflag==ON) {
	for (k=0;k<2;++k) {	
	  while (CV[k] <= -1.0*pi ) {
	    CV[k]+=2.0*pi;
	  }
	  while (CV[k] > 1.0*pi ) {
	    CV[k]-=2.0*pi;
	  }
	}
      }

      CV_dummy[0]=CV[0]+0.0001;
      CV_dummy[1]=CV[1];

      pe[j]=mixed_twod_Gaussian(CV,nyu_k,Sigma_k,pi_k,K,pi);
      pe[j]=log(pe[j]);

      de_ln_twoD_Mixed_Gaussian(CV,fe[j],nyu_k,Sigma_k,pi_k,K,pi);
      for (k=0;k<2;++k) {
      	fe[j][k]=-1.0*fe[j][k];
      }
    }

    z_string_CV(path,path_evoluted,fe,numpoint,2,dt);

    if (i>0) {
      maxdelta=fabs(pe[0]-pe_old[0]);
      for (j=1;j<numpoint;++j) {
	if (maxdelta<fabs(pe[j]-pe_old[j])) maxdelta=fabs(pe[j]-pe_old[j]);
      }
      printf("%lf\n",fabs(maxdelta));
      if (fabs(maxdelta)<threhold_path) 
	break;
    }

    for (j=0;j<numpoint;++j) pe_old[j]=pe[j];

    if (i%outinterval==0) {
      for (j=0;j<numpoint;++j) {
      	for (k=0;k<2;++k) {
      	  fprintf(outputfile,"%8.3lf ",path[j][k]);
      	}
      	fprintf(outputfile,"\n");
      }
      fprintf(outputfile,"\n");
    }
  }
  fclose(outputfile);

  pathfile=efopen(pathfilename,"w");
  for (j=0;j<numpoint;++j) {
    for (k=0;k<2;++k) {
      if (periodflag==ON) {
	while (path[j][k] <= -1.0*pi ) {
 	  path[j][k]+=2.0*pi;
	}
	while (path[j][k] > 1.0*pi ) {
	  path[j][k]-=2.0*pi;
	}
      }
      fprintf(pathfile,"%8.3lf ",path[j][k]);
    }
    fprintf(pathfile,"\n");
  }
  fclose(pathfile);

  KL=0.0;
  for (i=0;i<Nx;++i) {
    for (j=0;j<Ny;++j) {
      if (prob[i][j]!=0.0) {
	x[0]=minx+dx*i;
	x[1]=miny+dy*j;
	prob_app[i][j]=0.0;
	for (k=0;k<K;++k){
	  prob_app[i][j]+=pi_k[k]*twoD_Gaussian(x,nyu_k[k],Sigma_k[k],pi);
	}
      }
    }
  }

  sum=0.0;
  for (i=0;i<Nx;++i) {
    for (j=0;j<Ny;++j) {
      if (prob[i][j]!=0.0) {
	sum+=prob_app[i][j];
      }
    }
  }
  for (i=0;i<Nx;++i) {
    for (j=0;j<Ny;++j) {
      if (prob[i][j]!=0.0) {
	prob_app[i][j]=prob_app[i][j]/sum;
      }
    }
  }

  KL=0.0;
  for (i=0;i<Nx;++i) {
    for (j=0;j<Ny;++j) {
      if (prob[i][j]!=0.0)
	KL+=prob[i][j]*log(prob[i][j]/prob_app[i][j]);
    }
  }
  
  KLdfile=efopen(KLdfilename,"w");
  fprintf(KLdfile,"%12.8lf\n",KL);
  fclose(KLdfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename \n",progname);
}
