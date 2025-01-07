#include <stdio.h>
#include <stdlib.h>

void *emalloc(size_t n) {
  void *p;
  
  p = malloc(n);
  if (p == NULL)
    printf ("malloc of %u bytes failed:", n);
}

FILE *efopen(char *filename,char *flag){
  FILE *file;

  if((file=fopen(filename,flag))==NULL) {
    printf("open error about %s\n",filename);
    exit(1);
  }
  return file;
}
