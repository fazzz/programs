#
# Makefile for g_df (Ligand-Solvent)
#

TARGET = g_df

FC   = gfortran -g

FLAG = -O2

SRC = g_df.f90 gettop.f90 xdr.f90

OBJS = g_df.o gettop.o xdr.o

MOD_FILES = gettop.mod xdr.mod

.SUFFIXES: .f90

All:$(TARGET)

${TARGET}:${OBJS}
	$(FC) $(FLAG) -o $@ $(OBJS) \
	-L/Users/yamamoriyuu/onedrive/Research/work/programs/Fortran/xdrfile/lib -lxdrfile ;

g_df.o: gettop.o xdr.o
gettop.o:
xdr.o:

.f90.o:
	$(FC) -c $<;

.PHONY:clean install
clean:
	rm $(OBJS) $(MOD_FILES) $(TARGET)
install:
	cp $(TARGET) ../bin
