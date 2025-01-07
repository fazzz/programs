#
# Makefile for axis_rotation
#

TARGET = arot

FC   = gfortran

FLAG = -O2

#FLAG = -fbacktrace -ffpe-trap=invalid,zero,overflow -O -Wuninitialized -g -pg -Wall -pedantic -std=f95

SRC = bestfit.f90 pdb.f90 main.f90

OBJS = bestfit.o pdb.o main.o

MOD_FILES = bestfit.mod quaternion.mod pdb.mod

.SUFFIXES: .f90

All:$(TARGET)

${TARGET}:${OBJS}
	$(FC) $(FLAG) -o $@ $(OBJS) -llapack -lblas -L${HOME}/lib ;

main.o: bestfit.o pdb.o
bestfit.o :
pdb.o:

.f90.o:
	$(FC) $(FLAG) -c $<;

.PHONY:clean
clean:
	rm $(OBJS) $(MOD_FILES) $(TARGET)

