#
# Makefile for g_dih
#

TARGET = g_dih

FC   = gfortran

FLAG = -O2

SRC = g_dih-xtc.f90 xdr.f90

OBJS = g_dih-xtc.o xdr.o

MOD_FILES = xdr.mod

.SUFFIXES: .f90

All:$(TARGET)

${TARGET}:${OBJS}
	$(FC) $(FLAG) -o $@ $(OBJS) -L/Users/yamamoriyuu/onedrive/Research/work/programs/Fortran/xdrfile/lib -lxdrfile ;

g_dih-xtc.o: xdr.o
xdr.o:

.f90.o:
	$(FC) -c $<;

.PHONY:clean
clean:
	rm $(OBJS) $(MOD_FILES) $(TARGET)
