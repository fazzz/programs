#
# Makefile for pdb2gro
#

TARGET = pdb2gro

FC   = gfortran

#FLAG = -fbacktrace -ffpe-trap=invalid,zero,overflow \
-O -Wuninitialized -g -pg -Wall -pedantic -std=f95
FLAG = -O2

SRC = pdb2gro.f90 pdb.f90 gro.f90

OBJS = pdb2gro.o pdb.o gro.o

MOD_FILES = pdb.mod gro.mod

.SUFFIXES: .f90

All:$(TARGET)

${TARGET}:${OBJS}
	$(FC) $(FLAG) -o $@ $(OBJS) ;

pdb2gro.o: pdb.o gro.o
pdb.o:
gro.o:

.f90.o:
	$(FC) -c $<;

.PHONY:clean
clean:
	rm $(OBJS) $(MOD_FILES) $(TARGET)
