#
# Makefile for g_sasa
#

TARGET = sasa

FC   = gfortran

FLAG = -O2

SRC = sasa.f90

OBJS = sasa.o

.SUFFIXES: .f90

All:$(TARGET)

${TARGET}:${OBJS}
	$(FC) $(FLAG) -o $@ $(OBJS)

sasa.o:

.f90.o:
	$(FC) -c $<;

.PHONY:clean
clean:
	rm $(OBJS) $(MOD_FILES) $(TARGET)
