#
# Makefile for PEPCA Calculation Program
#

TARGET = g_pepca

FC   = gfortran

FLAG = -O2

SRC = main.f90 peppca.f90

OBJS = main.o pepca.o

MOD_FILES = pepca.mod

.SUFFIXES: .f90

All:$(TARGET)

${TARGET}:${OBJS}
	$(FC) $(FLAG) -o $@ $(OBJS) \
	-llapack -lblas -L${HOME}/lib;

main.o: pepca.o
pepca.o:

.f90.o:
	$(FC) -c $<;

.PHONY:clean install
clean:
	rm $(OBJS) $(MOD_FILES) $(TARGET)
install:
	cp $(TARGET) ../bin
