SHELL=/bin/sh
VPATH=

all:
	cd src && $(MAKE) all

install:
	cd src && $(MAKE) install
	@echo "Installation of pmemd complete"

depends:
	cd src && ../makef90depends > f90depends

clean:
	cd src && rm -f *.f90 *.o *.mod pmemd *.d work.pc*
