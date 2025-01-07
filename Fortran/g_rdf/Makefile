#
# Makefile for RDF (MolA-MolB)
#

TARGET = rdf-xtc_npt_MolA-MolB

FC   = gfortran

FLAG = -O2

SRC = rdf-xtc_npt_MolA-MolB.f90 xdr.f90

OBJS = rdf-xtc_npt_MolA-MolB.o xdr.o

MOD_FILES = xdr.mod

.SUFFIXES: .f90

All:$(TARGET)

${TARGET}:${OBJS}
	$(FC) $(FLAG) -o $@ $(OBJS) \
	-L/Users/yamamoriyuu/onedrive/Research/work/programs/Fortran/xdrfile/lib -lxdrfile ;

rdf-xtc.o: xdr.o
xdr.o:

.f90.o:
	$(FC) -c $<;

.PHONY:clean
clean:
	rm $(OBJS) $(MOD_FILES) $(TARGET)
