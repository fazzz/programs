#
# Makefile for g_catxtc
#

TARGET = g_catxtc

FC   = gfortran

FLAG = -O2

SRC = g_catxtc.f90 xdr.f90

OBJS = g_catxtc.o xdr.o

MOD_FILES = xdr.mod

.SUFFIXES: .f90

All:$(TARGET)

${TARGET}:${OBJS}
	$(FC) $(FLAG) -o $@ $(OBJS) -L/Users/yamamoriyuu/onedrive/Research/work/programs/Fortran/xdrfile/lib -lxdrfile ;

g_catxtc.o: xdr.o
xdr.o:

.f90.o:
	$(FC) -c $<;

.PHONY:clean
clean:
	rm $(OBJS) $(MOD_FILES) $(TARGET)
