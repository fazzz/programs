
#ifndef INCLUDE_OUT
#define INCLUDE_OUT

double *readOUT(FILE *inputfile,char *type,int *numstep);

double *readOUTwnum(FILE *inputfile,char *type,int numstep);

double *readOUTwab(FILE *inputfile,char *type,int *numstep,int numabondon);

double *readOUTwabwnum(FILE *inputfile,char *type,int numstep,int numabondon);

#endif
