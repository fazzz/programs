#ifndef INCLUDE_IOV2
#define INCLUDE_IOV2

#define AA 0
#define CA 1
#define HV 2

int IOV2_scantrj(FILE* trj,double ***crd,int MODE);

int IOV2_scantrj_wrange(FILE* trj,int inistep,int numstep,double ***crd,int MODE);

#endif
