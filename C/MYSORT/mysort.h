#ifndef INCLUDE_MS
#define INCLUDE_MS

#define ARRAY_SIZE(array) (sizeof(array)/sizeof(array[0]))

int compare_int(const void*a , const void*b);
int compare_double(const void*a , const void*b);
int compare_str(const void*a , const void*b);

myqsortint(int *array);
myqsortstr(char **array);


#endif
