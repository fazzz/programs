
#include <stdlib.h>
#include "mysort.h"

int compare_int(const void*a , const void*b){
  return (*(int*)a -*(int*)b);
}

int compare_double(const void*a , const void*b){
  return (*(double*)a -*(double*)b);
}

int compare_str(const void*a , const void*b){
  return (*(char**)a -*(char**)b);
}

myqsortint(int *array){
  return qsort(array,ARRAY_SIZE(array),sizeof(int),compare_int);
}

myqsortstr(char **array){
  return qsort(array,ARRAY_SIZE(array),sizeof(char *),compare_int);
}
