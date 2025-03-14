/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header: /home/case/cvsroot/amber11/AmberTools/src/netcdf/src/ncgen/escapes.c,v 9.1 2007/11/15 21:44:48 jmongan Exp $
 *********************************************************************/

#include <stdlib.h>
#include <netcdf.h>
#include "generic.h"
#include "ncgen.h"
#include "genlib.h"

/*
 * "Expands" valid escape sequences in yystring (read by lex) into the
 * apropriate characters in termstring.  For example, the two character
 * sequence "\t" in yystring would be converted into a single tab character
 * in termstring.  On return, termstring is properly terminated.
 */

void
expand_escapes(
     char *termstring,		/* returned, with escapes expanded */
     char *yytext,
     int yyleng)
{
    char *s, *t, *endp;
    
    yytext[yyleng-1]='\0';	/* don't copy quotes */
    /* expand "\" escapes, e.g. "\t" to tab character  */
    s = termstring;
    t = yytext+1;
    while(*t) {
	if (*t == '\\') {
	    t++;
	    switch (*t) {
	      case 'a':
		*s++ = '\007'; t++; /* will use '\a' when STDC */
		break;
	      case 'b':
		*s++ = '\b'; t++;
		break;
	      case 'f':
		*s++ = '\f'; t++;
		break;
	      case 'n':
		*s++ = '\n'; t++;
		break;
	      case 'r':
		*s++ = '\r'; t++;
		break;
	      case 't':
		*s++ = '\t'; t++;
		break;
	      case 'v':
		*s++ = '\v'; t++;
		break;
	      case '\\':
		*s++ = '\\'; t++;
		break;
	      case '?':
		*s++ = '\177'; t++;
		break;
	      case '\'':
		*s++ = '\''; t++;
		break;
	      case '\"':
		*s++ = '\"'; t++;
		break;
	      case 'x':
		t++; /* now t points to one or more hex digits */
		*s++ = (char) strtol(t, &endp, 16);
		t = endp;
		break;
	      case '0':
	      case '1':
	      case '2':
	      case '3':
	      case '4':
	      case '5':
	      case '6':
	      case '7':
		/* t now points to octal digits */
		*s++ = (char) strtol(t, &endp, 8);
		t = endp;
		break;
	      default:
		*s++ = *t++;
		break;
	    }
	} else {
	    *s++ = *t++;
	}
    }
    *s = '\0';
    return;
}
