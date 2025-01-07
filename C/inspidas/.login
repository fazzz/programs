#								| /
# This script is written by T. Tohdo in 1990			|/   *
#								..  )
#								******

stty tabs kill '^u' intr '^c' erase '^h' echoe echok

set savehist=40
set notify
set history = 20 
#set noclobber
set term=tigr40

setenv TZ Japan
setenv shell /bin/csh

if ( `tty` == "/dev/console" ) then
	set term=tigr40
else if ( `tty` == "/dev/tty1" ) then
	set term=tigr40
	setenv DISPLAY titan:0.0
else if ( $term == "vt100" ) then
	setenv DISPLAY titan:0.0
	stty erase '^h'
else
	stty erase '^h'
endif













