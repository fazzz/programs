#
# Makefile for get_num_fgro
#

TARGET = get_num_fgro

FC   = gfortran

#FLAG = -fbacktrace -ffpe-trap=invalid,zero,overflow \
-O -Wuninitialized -g -pg -Wall -pedantic -std=f95
FLAG = -O2

SRC = get_num_fgro.f90 gro.f90

OBJS = get_num_fgro.o gro.o

MOD_FILES = gro.mod

.SUFFIXES: .f90

All:$(TARGET)

${TARGET}:${OBJS}
	$(FC) $(FLAG) -o $@ $(OBJS) ;

get_num_fgro.o: gro.o
gro.o:

.f90.o:
	$(FC) -c $<;

.PHONY:clean
clean:
	rm $(OBJS) $(MOD_FILES) $(TARGET)
