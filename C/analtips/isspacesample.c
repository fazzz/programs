
#include <stdio.h>
#include <stdlib.h>

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,ii,jj;
  int c1,c2;
  char num1[20],num2[20];
  double a,b,d;

  char *progname;
  char *inpfile1name,*inpfile2name,*outfilename;
  FILE *inpfile1,*inpfile2,*outfile;

  char *line1,*line2,*dummy;
  char ch;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  int space_num=0;
  char str[]="Aa Bb 1-2:3 Cc     ";

  progname=argv[0];
  while((c=getopt(argc,argv,"h"))!=-1) {
    switch(c) {
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

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  inpfile1name  = *argv;
  inpfile2name = *++argv;
  outfilename = *++argv;

  inpfile1=/*e*/fopen(inpfile1name,"r");
  inpfile2=/*e*/fopen(inpfile2name,"r");
  outfile=/*e*/fopen(outfilename,"w");

  while ( (c1=getline(&line1,&len,inpfile1))!=-1 && (c2=getline(&line2,&len,inpfile2))!=-1 ) {

    printf("kigou:");
    for(i=0;i<20;++i)
      if(ispunct(str[i])!=0) printf("%c",str[i]);

    printf("\nkuhaku:");
    for(i=0;i<20;++i)
      if(isspace(str[i])!=0) space_num++;
    printf("%d\n",space_num);
  }

  fclose(inpfile1);
  fclose(inpfile2);
  fclose(outfile);


  return 0;
}

void USAGE(char *progname) {
  printf("-h -- help\n");
  printf("USAGE: %s inpfile2name inpfile2name outfilename \n", progname);
}
