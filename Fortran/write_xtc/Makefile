#
# Makefile for example
#

TARGET = woxtc

FC   = gfortran -g

FLAG = -O2

SRC = woxtc.f90 xdr.f90

OBJS = woxtc.o xdr.o

MOD_FILES = xdr.mod

.SUFFIXES: .f90

All:$(TARGET)

${TARGET}:${OBJS}
	$(FC) $(FLAG) -o $@ $(OBJS) \
	-L/Users/yamamoriyuu/onedrive/Research/work/programs/Fortran/xdrfile/lib -lxdrfile ;

woxtc.o: xdr.o
xdr.o:

.f90.o:
	$(FC) -c $<;

.PHONY:clean install
clean:
	rm $(OBJS) $(MOD_FILES) $(TARGET)

