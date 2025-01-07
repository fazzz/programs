#include <stdio.h>
#include <stdlib.h>

#include "EF.h"
//#include "gc.h"

#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;}

FILE *efopen(char *filename,char *flag){
  FILE *file;

  if((file=fopen(filename,flag))==NULL) {
    printf("open error about %s\n",filename);
    exit(1);
  }
  return file;
}

void *emalloc(size_t n){
  void *p;

  if ((p=malloc(n))==NULL) {
    printf("open about allocation\n");
    exit(1);
  }

  return p;
}

void *gcemalloc(size_t n){
  void *p;

  if ((p=GC_MALLOC(n))==NULL) {
    printf("open about allocation\n");
    exit(1);
  }

  return p;
}

void *ecalloc(size_t t,size_t size){
  void *p;

  if ((p=calloc(t,size))==NULL) {
    printf("open about allocation\n");
    exit(1);
  }

  return p;
}

/******************************************/
/* void *gcecalloc(size_t t,size_t size){ */
/*   void *p;				  */
/* 					  */
/*   if ((p=GC_CALLOC(t,size))==NULL) {	  */
/*     printf("open about allocation\n"); */
/*     exit(1);				  */
/*   }					  */
/* 					  */
/*   return p;				  */
/* }					  */
/******************************************/

void *gcerealloc(void *ptr, size_t n){
  void *p;

  if ((p=GC_REALLOC(ptr, n))==NULL) {
    printf("open about allocation\n");
    exit(1);
  }

  return p;
}

/********************************************************************************/
/* char *egetenv(const char *name) {					        */
/*   char *env;								        */
/* 									        */
/*   if ((env=getenv(name))==NULL) {					        */
/*     printf("error; there is not eiviromental variable %s\n",name);	        */
/*     exit(1);								        */
/*   }									        */
/*   else {								        */
/*     return env;							        */
/*   }									        */
/* }									        */
/* 									        */
/* int enc_create(const char *path, int cmode, int *ncidp){		        */
/*   int c;								        */
/* 									        */
/*   if ((c = nc_create(path,cmode,ncidp)))				        */
/*     ERR(c);								        */
/* }									        */
/* 									        */
/* int enc_def_dim(int ncid, const char *name, size_t len, int *idp){	        */
/*   int c;								        */
/* 									        */
/*   if ((c = nc_def_dim(ncid,name,len,idp)))				        */
/*     ERR(c);								        */
/* }									        */
/* 									        */
/* int enc_def_var(int ncid, const char *name, nc_type xtype, int ndims,        */
/* 		const int *dimidsp, int *varidp){			        */
/*   int c;								        */
/* 									        */
/*   if ((c = nc_def_var(ncid,name,xtype,ndims,dimidsp,varidp)))	        */
/*     ERR(c);								        */
/* }									        */
/* 									        */
/* int enc_inq_varid(int ncid, const char *name, int *varidp){		        */
/*   int c;								        */
/* 									        */
/*   if ((c = nc_inq_varid(ncid,name,varidp)))				        */
/*     ERR(c);								        */
/* }									        */
/* 									        */
/* int enc_open(const char *path, int mode, int *ncidp){		        */
/*   int c;								        */
/* 									        */
/*   if ((c = nc_open(path,mode,ncidp)))				        */
/*     ERR(c);								        */
/* }									        */
/* 									        */
/* int encopen(const char* path, int mode) {				        */
/*   int c;								        */
/* 									        */
/*   if ((c = ncopen(path,mode)))					        */
/*     ERR(c);								        */
/* }									        */
/* 									        */
/* int encclose(int ncid){						        */
/*   int c;								        */
/* 									        */
/*   if ((c = ncclose(ncid)))						        */
/*     ERR(c);								        */
/* }									        */
/********************************************************************************/
