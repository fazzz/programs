#
# Makefile for Entropy Calculation Program
#

TARGET = calEntropy

FC   = gfortran
prefix = /Users/yamamoriyuu/onedrive/Research/work/programs/Fortran/Entropy
exec_prefix = ${prefix}
bindir = ${exec_prefix}/bin

#FLAG = -fbacktrace -ffpe-trap=invalid,zero,overflow -O -Wuninitialized -g -pg -Wall -pedantic -std=f95
FLAG = -O2

SRC = bestfit.f90 dgeev0823.f90

OBJS = mpiproc.o xdr.o

MOD_FILES = graphical_Lasso.mod standardlize.mod

.SUFFIXES: .f90

All:$(TARGET)

${TARGET}:${OBJS}
	$(FC) $(FLAG) -o $@ $(OBJS) ;

SpStLg_main.o:
graphical_Lasso.o:
standardlize.o:

.f90.o:
	$(FC) -c $<;

.PHONY:clean install
clean:
	rm $(OBJS) $(MOD_FILES) $(TARGET) 
install:
	cp $(TARGET) ${exec_prefix}/bin
