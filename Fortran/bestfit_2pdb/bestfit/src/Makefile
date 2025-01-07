#
# Makefile for bestfit
#

TARGET = bestfit

FC   = gfortran

FLAG = -O2

#FLAG = -fbacktrace -ffpe-trap=invalid,zero,overflow -O -Wuninitialized -g -pg -Wall -pedantic -std=f95

SRC = bestfit.f90 ndx.f90 main.f90

OBJS = bestfit.o ndx.o main.o

MOD_FILES = bestfit.mod quaternion.mod ndx.mod

.SUFFIXES: .f90

All:$(TARGET)

${TARGET}:${OBJS}
	$(FC) $(FLAG) -o $@ $(OBJS) -llapack -lblas -L${HOME}/lib ;

bestfit_2pdb_wndx.f90.o: bestfit.o ndx.o
bestfit.o :

.f90.o:
	$(FC) $(FLAG) -c $<;

.PHONY:clean
clean:
	rm $(OBJS) $(MOD_FILES) $(TARGET)
install:
	cp $(TARGET) ../bin

