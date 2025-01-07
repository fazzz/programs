#								| /
# This script is written by T. Tohdo in 1990			|/   *
#								..  )
#								******

set path=(. ~/bin /usr/local/bin /bin /usr/bin /usr/ucb /usr/X11/bin \
          ~/tools)

set prompt="% "
set history=100
set cdpath=(~ /tmp1/nisigaki /tmp1)

umask 022

alias	ls	"ls -F -x" 
alias	df	"df -b"
alias	lpr	"lpb"

alias   cd      'cd \!*;set cwd=`pwd`;d'
alias   d       echo '"'`hostname`' :$cwd"'

if      ( $?Setup ) then
        set gamichan=($Setup)
        unsetenv Setup
        $gamichan
        unset gamichan
endif

limit coredumpsize 0

