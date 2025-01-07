#
# Makefile for Sparse Structure Learning-Principle Component Analysis Program
#

TARGET = g_sslpca

FC   = gfortran
bindir = ../bin

FLAG = -fbacktrace -ffpe-trap=invalid,zero,overflow -O -Wuninitialized -g -pg -Wall -pedantic -std=f95
#FLAG = -O2

SRC = main.f90 bestfit.f90 graphical_Lasso.f90 standardlize.f90 gettop.f90 pca.f90 xdr.f90

OBJS = main.o bestfit.o graphical_Lasso.o standardlize.o gettop.o pca.o xdr.o

MOD_FILES = bestfit.mod graphical_Lasso.mod standardlize.mod gettop.mod pca.mod xdr.mod

.SUFFIXES: .f90

All:$(TARGET)

${TARGET}:${OBJS}
	$(FC) $(FLAG) -o $@ $(OBJS) \
	-llapack -lblas -L${HOME}/lib \
	-L/Users/yamamoriyuu/onedrive/Research/work/programs/Fortran/xdrfile/lib -lxdrfile ;

main.o: bestfit.o graphical_Lasso.o standardlize.o gettop.o pca.o xdr.o
graphical_Lasso.o:
standardlize.o:
bestfit.o:
gettop.o:
pca.o: bestfit.o
xdr.o:

.f90.o:
	$(FC) -c $<;

.PHONY:clean install
clean:
	rm $(OBJS) $(MOD_FILES) $(TARGET) 
install:
	cp $(TARGET) ../bin
