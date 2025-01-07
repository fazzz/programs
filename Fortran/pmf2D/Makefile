#
# Makefile for pmf2d
#

TARGET = pmf2d

FC   = gfortran

FLAG = -O2

SRC = pmf2d.f90

.SUFFIXES: .f90

All:$(TARGET)

${TARGET}:${OBJS}
	$(FC) $(FLAG) $(SRC) -o $@ $(OBJS);

.f90.o:
	$(FC) $<;

.PHONY:clean
clean:
	rm $(TARGET)
