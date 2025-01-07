#!/bin/awk

BEGIN{i=0}
{++i;if((i%intv)==ame){print $0}}
