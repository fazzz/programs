#
# Makefile for clust_simple
#

TARGET = clust_simple

FC   = gfortran

#FLAG = -fbacktrace -ffpe-trap=invalid,zero,overflow \
-O -Wuninitialized -g -pg -Wall -pedantic -std=f95
FLAG = -O2

SRC = clustering_threshold.f90

OBJS = clustering_threshold.o

.SUFFIXES: .f90

All:$(TARGET)

${TARGET}:${OBJS}
	$(FC) $(FLAG) -o $@ $(OBJS) ;

clusteing_threshold.o:

.f90.o:
	$(FC) -c $<;

.PHONY:clean
clean:
	rm $(OBJS) $(TARGET)
