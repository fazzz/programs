#!/bin/gawk -f

BEGIN{
  l=0
}

{
    if (l>0) {
	for (i=1;i<NF;++i) {
	    printf("%lf \n",$(i+1)) >> i".dihed"
	}
    }
    ++l
}


