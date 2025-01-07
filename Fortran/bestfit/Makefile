#
# Makefile for bestfit_2pdb.f90
#

TARGET = bestfit

FC   = gfortran

#FLAG = -fbacktrace -ffpe-trap=invalid,zero,overflow \
-O -Wuninitialized -g -pg -Wall -pedantic -std=f95
FLAG = -O2

SRC = bestfit_2pdb.f90 bestfit.f90 pdb.f90 gettop.f90

OBJS = bestfit_2pdb.o bestfit.o pdb.o gettop.o

MOD_FILES = bestfit.mod quaternion.mod pdb.mod gettop.mod

.SUFFIXES: .f90

All:$(TARGET)

${TARGET}:${OBJS}
	$(FC) $(FLAG) -o $@ $(OBJS) -llapack -lblas -L${HOME}/lib ;

bestfit_2pdb.o: bestfit.o pdb.o gettop.o
bestfit.o :
pdb.o:
gettop.o:

.f90.o:
	$(FC) -c $<;

.PHONY:clean
clean:
	rm $(OBJS) $(MOD_FILES) $(TARGET)
install:
	cp $(TARGET) ~/mybin

