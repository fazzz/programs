#
# Makefile for wham_phi_psi
#

TARGET = wham_phi-psi

FC   = gfortran

#FLAG = -fbacktrace -ffpe-trap=invalid,zero,overflow -O -Wuninitialized -g -pg -Wall -pedantic -std=f95
FLAG = -O2

SRC = wham_phi-psi.f90 reallocate.f90

OBJS =  wham_phi-psi.o reallocate.o

MOD_FILES = reallocate.mod

.SUFFIXES: .f90

all:$(TARGET)

${TARGET}:${OBJS}
	$(FC) $(FLAG) -o $@ $(OBJS) ;

wham_phi-psi.o: reallocate.o
reallocate.o:

.f90.o:
	$(FC) -c $<;

.PHONY:clean
clean:
	rm $(OBJS) $(MOD_FILES) $(TARGET)
