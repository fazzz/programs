
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "IO.h"
#include "EF.h"
#include "PT.h"

int main(int argc, char *argv[]) {
  int i,j,k,n;
  int numstep,numatom,numcrd,interval;
  double *data;
  
  char *inputfilename,*inputfilename2,*inputfilename3;
  char outputfilename[100],*outputfilebasename;
  FILE *inputfile,*inputfile2,*inputfile3,*outputfile;
  
  if (argc < 4) {
    printf("USAGE: ./pickcoord.exe inputfilename(traj) inputfilename2(cond) inputfilename3(parm) outputfilebasename\n");
    exit(1);
  }
  inputfilename =  *++argv;
  inputfilename2 = *++argv;
  inputfilename3 = *++argv;
  outputfilebasename = *++argv;

  inputfile2=efopen(inputfilename2,"r");
  fscanf(inputfile2,"%d",&numstep);
  fscanf(inputfile2,"%d",&numcrd);
  fclose(inputfile2);

  inputfile3=efopen(inputfilename3,"r");
  readParmtop(inputfile3);
  numatom=AP.NATOM;
  fclose(inputfile3);

  data=(double *)emalloc(sizeof(double)*numatom*3);
  n=0;
  interval=(int)numstep/numcrd;  
  inputfile=efopen(inputfilename,"r");
  for (i = 0; i < numstep; ++i){
    io_scanconf(inputfile,numatom,data,'x');
    if (i%interval==0) {
      ++n;
      sprintf(outputfilename,"%s%d.crd",outputfilebasename,n);
      outputfile=efopen(outputfilename,"w");
      /*********************************************************/
      /* for (j=0;j<numatom;++j) {			       */
      /* 	for (k=0;k<3;++k)			       */
      /* 	  fprintf(outputfile,"%10.6lf",data[j*3+k]);   */
      /* 	fprintf(outputfile,"\n");		       */
      /* }						       */
      /* fprintf(outputfile,"\n");			       */
      /*********************************************************/
      io_outputconf_Amberform(outputfile,numatom,data);
  //      io_outputconf(outputfile,numatom,data,'x');
      fclose(outputfile);
    }
  }
  fclose(inputfile);
  
  return 0;
}

