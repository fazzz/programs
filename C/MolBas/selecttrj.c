
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "IO.h"
#include "PT.h"
#include "EF.h"

int main(int argc, char *argv[]) {
  int i,j;
  int numstep,numatom,numdihed;
  int *num;
  int **atomdihedpairs;
  double pca,*data;
  char *line;
  size_t len=0;
  
  char *inputfilename1,*inputfilename2,*inputfilename3,*inputfilename4,*outputfilename,*outputfilename2;
  FILE *inputfile1,*inputfile2,*inputfile3,*inputfile4,*outputfile,*outputfile2;
  
  if (argc < 6) {
    printf("USAGE: %s inputfilename1(dihed) inputfilename2(pepca) inputfilename3(cond) inputfilename4(parm) outputfilename1(p) outputfilename2(m) \n",argv[0]);
    exit(1);
  }
  inputfilename1 = *++argv;
  inputfilename2 = *++argv;
  inputfilename3 = *++argv;
  inputfilename4 = *++argv;
  outputfilename = *++argv;
  outputfilename2 = *++argv;
  
  inputfile3=efopen(inputfilename3,"r");
  fscanf(inputfile3,"%d",&numstep);
  fscanf(inputfile3,"%d",&numdihed);
  fclose(inputfile3);

  inputfile4=efopen(inputfilename4,"r");
  readParmtop(inputfile4);
  fclose(inputfile4);
  /***************************************************************/
  /* num=(int *)gcemalloc(sizeof(int)*5);			 */
  /* atomdihedpairs[0]=(int *)emalloc(sizeof(int)*4); // PHI,PSI */
  /* atomdihedpairs[1]=(int *)emalloc(sizeof(int)*4); // OMEGA	 */
  /* atomdihedpairs[2]=(int *)emalloc(sizeof(int)*4); // KI	 */
  /* atomdihedpairs[3]=(int *)emalloc(sizeof(int)*8); // ACE	 */
  /* atomdihedpairs[4]=(int *)emalloc(sizeof(int)*8); // NME	 */
  /* numdihed=readdihedpairs(atomdihedpairs,num);		 */
  /***************************************************************/

  data=(double *)gcemalloc(sizeof(double)*(numdihed+1));
  
  inputfile1=efopen(inputfilename1,"r");
  getline(&line,&len,inputfile1);
  inputfile2=efopen(inputfilename2,"r");
  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");
  for (i=0;i<numstep;++i) {
    io_scan_data(inputfile1,numdihed+1,data,'x');
    fscanf(inputfile2,"%lf",&pca);

    if (pca>0)
      io_outputdata(outputfile,numdihed+1,data);
    else
      io_outputdata(outputfile2,numdihed+1,data);
  }
  
  fclose(inputfile1);
  fclose(inputfile2);
  fclose(outputfile);
  fclose(outputfile2);

  /****************************/
  /* free(atomdihedpairs[0]); */
  /* free(atomdihedpairs[1]); */
  /* free(atomdihedpairs[2]); */
  /* free(atomdihedpairs[3]); */
  /* free(atomdihedpairs[4]); */
  /****************************/
  
  return 0;
}

