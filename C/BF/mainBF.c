
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bestfit.h"
#include "PT.h"
#include "EF.h"

int main( int argc, char *argv[] ) 
{
  int numatom,numstep;
  double *coord,*coord_ref,*coord_ave;
  double *mass;

  char *inputfilename1,*inputfilename2,*inputfilename3,*outputfilename;
  FILE *inputfile1,*inputfile2, *outputfile;
  
  if (argc < 5) {
    printf("USAGE: BFC inputfilename1(traj) inputfilename2(cond) inputfilename3(parm) outputfilename(traj_bf)\n");
    exit(1);
  }
  inputfilename1 = *++argv;
  inputfilename2 = *++argv;
  inputfilename3 = *++argv;
  outputfilename = *++argv;
  
  inputfile2=efopen(inputfilename2,"r");
  fscanf("%d",&numstep);
  fclose(inputfile2);

  inputfile3=efopen(inputfilename3,"r");
  readParmTop(inputfilename3);
  fclose(inputfile3);
  numatom=AP.NATOM;
  for (i=0;i<numatom;++i)
    mass[i]=AP.AMASS[i];
  
  inputfile1=efopen(inputfilename1,"r");
  io_scantrj_Amber_form(traj,numatom,numstep);
  fclose(inputfile1);

  rmsd=0.0;
   do {
    rmsd_old=rmsd;
    ave_coord(coord_ref,traj,numatom,numstep);
    rmsd=bestfit_trajectry(numatom,numstep,mass,coord_ref,rmsd_trj,traj);
    printf("%d cyc: rmsd=%lf\n",cyc,rmsd);
    ++cyc;
  } while (abf(rmsd_old-rmsd) < LIMIT)
  
  outputfile=efopen(outputfilename,"w");
  IO_outputtraj(traj,numatom,numstep);
  fclose(outputfile);
  
  return 0;
}

