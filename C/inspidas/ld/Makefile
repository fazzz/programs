# Makefile : originally written by F.Hirata, revised edition by K.Oda
#
OBJ= \
flxsub.o \
graph.o \
iolibdas.o \
ioroutn.o \
mathlib.o \
minim.o \
mc.o \
nma.o \
nmanal.o \
precep.o \
prein.o \
utilitydas.o \
main.o

TARGET= inspidas
FORT= gfortran
FC  = gfortran
FORT= g77
FC  = g77
FORT= ifort
FC  = ifort
CFLAGS= -f77rtl -intconstant -check all -debug all
CFLAGS= -check all -traceback -g
CFLAGS= -check all -debug minimal
CFLAGS= -traceback -g -fpe0
CFLAGS= -traceback  -fpe0
CFLAGS= -intconstant -fpconstant -f77rtl
CFLAGS=  -traceback  -fpe0 -intconstant -fpconstant -f77rtl
CFLAGS=  -O0 -g
CFLAGS= -check all -debug all -traceback -fpe0
CFLAGS= -debug all -traceback -fpe0
CFLAGS= -traceback -fpe0

$(TARGET): $(OBJ)
	$(FORT) -o $(TARGET) $(OBJ) $(CFLAGS)
.f.o:
	$(FORT) $(CFLAGS) -c $<
clean :
	rm *.o
