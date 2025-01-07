# Makefile : originally written by F.Hirata, revised edition by K.Oda
#
OBJS= \
/labo/qclib/inspidas/ld/flxsub.o \
/labo/qclib/inspidas/ld/graph.o \
/labo/qclib/inspidas/ld/iolibdas.o \
/labo/qclib/inspidas/ld/ioroutn.o \
/labo/qclib/inspidas/ld/mathlib.o \
/labo/qclib/inspidas/ld/minim.o \
/labo/qclib/inspidas/ld/mc.o \
nma.o \
/labo/qclib/inspidas/ld/nmanal.o \
/labo/qclib/inspidas/ld/precep.o \
/labo/qclib/inspidas/ld/prein.o \
/labo/qclib/inspidas/ld/utilitydas.o \
/labo/qclib/inspidas/ld/main.o


TARGET= nma

CFLAGS= -O2

main : $(OBJS)
	fc $(CFLAGS) -o inspidas $(OBJS)

all :
	fc $(CFLAGS) -c *.f

$(TARGET).o : $(TARGET).f
	fc $(CFLAGS) -c $(TARGET).f

$(TARGET).f : $(TARGET)
	cut.cttn $(TARGET) >$(TARGET).f
clear :
	rm *.o
