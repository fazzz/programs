
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PDBedit.h"

#include "EF.h"

#define start 0
#define fin 1

int USAGE(char *progname);
int USAGE2(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,d;
  int w;
  int initialatom,finalatom;
  int len1,len2;
  int numatom;

  int flagsof,flagnuma,flagnum;

  int nums,n;
  int numf;

  int replaceflag=OFF,replaceallflag=OFF,deleteflag=OFF,pickupflag=OFF,putOcuflag=OFF;
  char flag;

  char string[1000];

  PDBF PDB,PDBundo;

  char ATOMNAME[4],ATOMNAME2[4];
  char ATOMNAMEreplace[4],ATOMNAMEreplaceall[4],ATOMNAMEdelete[4],ATOMNAMEpickup[4];
  double occupancy;
  
  char *pdbfilename1,*outputfilename1,pdbfilename[1000],outputfilename[1000];
  FILE *pdbfile,*outputfile;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"rep",0,NULL,'r'},
    {"repa",1,NULL,'a'},
    {"del",1,NULL,'d'},
    {"pick",1,NULL,'p'},
    {"putOcu",1,NULL,'O'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hr:a:d:p:O:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'r':
      if (replaceallflag==ON) {
	replaceflag=ON;
	replaceallflag=OFF;
	for (i=0;i<4;++i) ATOMNAMEreplace[i]=optarg[i];
      }
      else {
	printf("error:\n");
	exit(1);
      }
      break;
    case 'a':
      replaceallflag=ON;
      for (i=0;i<4;++i) ATOMNAMEreplaceall[i]=optarg[i];
      break;
    case 'd':
      deleteflag=ON;
      for (i=0;i<4;++i) ATOMNAMEdelete[i]=optarg[i];
      break;
    case 'p':
      pickupflag=ON;
      for (i=0;i<4;++i) ATOMNAMEpickup[i]=optarg[i];
      break;
    case 'O':
      putOcuflag=ON;
      occupancy=atof(optarg);
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

  if (argc < 1) {
    printf("PDBfilename:\n");
    scanf("%s",pdbfilename);
    i=strlen(pdbfilename);
    pdbfilename[i]='\0';
    pdbfile=efopen(pdbfilename,"r");
    readPDBatomnum(pdbfile,&numatom);
    pdbfile=efopen(pdbfilename1,"r");
    readPDB(pdbfile,PDB,numatom);
  }
  else {
    pdbfilename1=*argv;
    pdbfile=efopen(pdbfilename1,"r");
    readPDBatomnum(pdbfile,&numatom);

    PDB.PDBa=(struct PDBatom *)gcemalloc(sizeof(struct PDBatom)*numatom);
    for (i=0;i<numatom;++i) {
      PDB.PDBa[i].HETEROflag=0;
      PDB.PDBa[i].serial=0;
      PDB.PDBa[i].ChainID=0;
      PDB.PDBa[i].resSeq=0;
      for (j=0;j<3;++j) PDB.PDBa[i].coord[j]=0.0;
      PDB.PDBa[i].occupancy=0.0;
      PDB.PDBa[i].tempfact=0.0;
    }
    pdbfile=efopen(pdbfilename1,"r");
    readPDB(pdbfile,PDB,numatom);
    PDB.numatom=numatom;
  }

  initialatom=0;
  finalatom=numatom;

  if (argc > 1) {
    outputfilename1 = *++argv;
    if (deleteflag==ON) {
      for (i=0;i<finalatom;++i) {
	w=delete_ATOMNAME(PDB,ATOMNAMEdelete,i);
	i-=w;
	finalatom-=w;
	numatom-=w;
      }
    }
    if (pickupflag==ON) {
      for (i=0;i<numatom;++i) {
	w=pickup_ATOMNAME(PDB,ATOMNAMEpickup,i);
	i-=w;
	finalatom-=w;
	numatom-=w;
      }
    }
    if (replaceflag==ON) {
      for (i=0;i<numatom;++i) replace_ATOMNAME(PDB,ATOMNAMEreplaceall,ATOMNAMEreplace,i);
    }
    if (replaceallflag==ON) {
      for (i=0;i<numatom;++i) replaceall_ATOMNAME(PDB,ATOMNAMEreplaceall,i);
    }
    if (putOcuflag==ON) {
      for (i=0;i<numatom;++i) put_occupancy(PDB,occupancy,i);
    }

    outputfile=efopen(outputfilename1,"w");
    PDB.numatom=numatom;
    writPDB(outputfile,PDB);
    fclose(outputfile);

    exit(1);
  }

  printf(">>\n");
  while (scanf("%c",&flag)!=-1) {
    switch (flag) {
    case 'i':
      printf("initialatom=:\n");
      scanf("%d",&initialatom);
      
      printf(">>\n");
      break;
    case 'f':
      printf("finalatom=:\n");
      scanf("%d",&finalatom);

      printf(">>\n");
      break;
    case 'd':
      printf("ATOMNAME to delete:\n");
      scanf("%4s",ATOMNAME);

      memcpy(&PDBundo,&PDB,sizeof(PDBF));

      for (i=initialatom;i<finalatom;++i) {
	w=delete_ATOMNAME(PDB,ATOMNAME,i);
	i-=w;
	finalatom-=w;
	numatom-=w;
      }

      printf(">>\n");
      break;
    case 'r':
      printf("ATOMNAME1 :\n");
      scanf("%4s",ATOMNAME);

      printf("ATOMNAME2:\n");
      scanf("%4s",ATOMNAME2);

      memcpy(&PDBundo,&PDB,sizeof(PDBF));

      for (i=initialatom;i<finalatom;++i) replace_ATOMNAME(PDB,ATOMNAME,ATOMNAME2,i);

      printf(">>\n");
      break;
    case 'l':
      printf("ATOMNAME :\n");
      scanf("%4s",ATOMNAME);

      memcpy(&PDBundo,&PDB,sizeof(PDBF));

      for (i=initialatom;i<finalatom;++i) replaceall_ATOMNAME(PDB,ATOMNAME,i);
      printf(">>\n");
      break;
    case 'p':
      printf("ATOMNAME to pick up:\n");
      scanf("%4s",ATOMNAME);

      memcpy(&PDBundo,&PDB,sizeof(PDBF));

      for (i=initialatom;i<finalatom;++i) {
	w=pickup_ATOMNAME(PDB,ATOMNAME,i);
	i-=w;
	finalatom-=w;
	numatom-=w;
      }
      printf(">>\n");
      break;
    case 's':
      printf("number of serial=:\n");
      scanf("%s",string);

      flagsof=start;     
      nums=0;
      numf=0;
      w=0;
      while ((c=string[w])!=-'\0'){
	++w;
	if (c=='\n') {
	  break;
	}
	else if (c==' ') {
	  d=0;
	  n=-1;
	}
	else if (isdigit(c)) {
	  flagnum=ON;
	  d=(c-'0');
	  ++n;
	  if (flagsof==start) {
	    nums=nums*10;
	    nums+=d;
	  }
	  if (flagsof!=start) {
	    numf=numf*10;
	    numf+=d;
	  }
	}
	else if (c=='-')
	  flagsof=fin;
	else
	  break;
      }

      if (flagsof==start) {
	if ( i > 0 )
	  i=nums-1;
	show_line(PDB,i);
      }
      else {
	for (i=nums-1;i<numf;++i ) {
	  show_line(PDB,i);
	}
      }

      printf(">>\n");
      break;
    case 'a':
      printf("outputfilename:\n");
      scanf("%s",outputfilename);
      len=strlen(outputfilename);
      outputfilename[len]=0;

      outputfile=efopen(outputfilename,"w");
      PDB.numatom=numatom;
      writPDB(outputfile,PDB);
      fclose(outputfile);

      exit(1);
    case 'u':
      memcpy(&PDB,&PDBundo,sizeof(PDBF));

      break;
    case 'h':
      USAGE2(progname);
      printf(">>\n");
      break;
    default:
      //      USAGE2(progname);
      printf(">>\n");
      break;
    }
  }

  return 0.0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("[-rep] \n");
  printf("[-repa] \n");
  printf("[-del] \n");
  printf("[-pick] \n");
  printf("[-putOcu] \n");
  printf("%s outputfilename \n",progname);
}

int USAGE2(char *progname) {
  printf("USAGE:\n");
  printf("[-i] initialatom\n");
  printf("[-f] finalatom\n");
  printf("[-d] delete ATOM\n");
  printf("[-r] replace atom\n");
  printf("[-l] replace all atom\n");
  printf("[-p] pick up ATOM name\n");
  printf("[-s] show line\n");
  printf("[-u] undo\n");
  printf("[-a] save PDB file name\n\n");
}
