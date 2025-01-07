#ident	"@(#)Profile.sh	1.2"	6/20/88

#Copyright (C) 1988, by Ardent Computer Corp.
#        All Rights Reserved
#This program is a  trade secret of  Ardent Computer Corp. and it is not to be
#reproduced, published, disclosed  to others, copied, adapted, distributed,
#or displayed without the prior authorization of Ardent Computer Corp.  Licensee
#agrees  to  attach  or  embed  this  Notice on all copies  of the program,
#including partial copies or modified versions thereof.

# Try copying the following to .profile in your home directory.

stty erase '^h' kill '^u' intr '^c' echoe -echok tab0
PATH=$PATH:/usr/ucb:/usr/X11/bin
TERM=tigr40
MAIL=/usr/mail/$LOGNAME
export TERM PATH TZ MAIL
