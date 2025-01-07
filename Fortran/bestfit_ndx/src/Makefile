#
# Makefile for bestfit_2pdb_wndx.f90
#

TARGET = bestfit_wndx

FC   = gfortran

FLAG = -O2

#FLAG = -fbacktrace -ffpe-trap=invalid,zero,overflow -O -Wuninitialized -g -pg -Wall -pedantic -std=f95

SRC = bestfit.f90 pdb.f90 gettop.f90 ndx.f90  bestfit_2pdb_wndx.f90

OBJS = bestfit.o pdb.o gettop.o ndx.o bestfit_2pdb_wndx.o

MOD_FILES = bestfit.mod quaternion.mod pdb.mod gettop.mod ndx.mod

.SUFFIXES: .f90

All:$(TARGET)

${TARGET}:${OBJS}
	$(FC) $(FLAG) -o $@ $(OBJS) -llapack -lblas -L${HOME}/lib ;

bestfit_2pdb_wndx.f90.o: bestfit.o pdb.o gettop.o ndx.o
bestfit.o :
ndx.o:
pdb.o:
gettop.o:

.f90.o:
	$(FC) $(FLAG) -c $<;

.PHONY:clean
clean:
	rm $(OBJS) $(MOD_FILES) $(TARGET)
install:
	cp $(TARGET) ../bin

