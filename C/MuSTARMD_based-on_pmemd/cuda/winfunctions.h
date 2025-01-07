// Matthew Clark
// July 2011

// Pmemd may use native linux functions and libraries which
// are not present under Windows, so provide equivalent functions

int ffs ( int i )
{
	if ( i == 0 ) return 0;
	int j = 0;
	while ( 1 )
	{
		if ( (i>>j) & 1 == 1 ) break;
		j++;
	}
	return j + 1;
}
