#include <stdio.h>
#include <stdlib.h>

#include "MD.h"

void getoption(int num, char *option)
{
	int i;

	switch(num) {
		case 'c': 
			InpfilCLUST = option;
		case 'r':
			InpfilCOORD = option;
		case 't':
			InpfilSEQ = option;
		case 's':
			InpfilTOP = option;
		default:
			printf("too many options !!! \n");
			exit(1);
		}
}
