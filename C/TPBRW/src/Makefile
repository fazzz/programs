#
# Makefile for Trigonometric polynomial base reweigthing program
#

HOME=../
SRC=.
BIN=../bin
LIBDIR=-L/home/yamamori/work/programs/GBRW/lib -L/home/yamamori/work/programs/ABAMD_2014_05_01/lib
INCDIR=-I/home/yamamori/work/programs/GBRW/include -I/home/yamamori/work/programs/ABAMD_2014_05_01/include

# definitions
#CC    = icc -pg -O3 -tpp6 -ipo
#CC    = icc -check=conversion,stack,uninit
#CC    = icc
CC    = gcc -g
#CC     = gcc

CFLAG = #-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -m64 -O2

OBJS =  main_TPBRW.o Optimize_BFGS_TPBRW.o Simpson_integ_TPBRW.o EF.o

LIBS =  Optimize_BFGS_TPBRW.h Simpson_integ_TPBRW.h EF.h

EXLIBS=-llbfgs -lm -lgc

TARGET =TPBRW

.c.o:
	$(CC) $(INCDIR) -c $(CFLAG) $<;

all: TPBRW install

# rules of generations
TPBRW:  $(OBJS) $(LIBS)
	$(CC) $(CFLAG) -o $@ $(OBJS) $(LIBDIR) $(INCDIR) $(EXLIBS) ;

install: 
	cp $(TARGET) $(BIN) ; 

main_TPBRW.o: Optimize_BFGS_TPBRW.h Simpson_integ_TPBRW.h EF.h
Optimize_BFGS_TPBRW.o: Optimize_BFGS_TPBRW.h Simpson_integ_TPBRW.h EF.h
Simpson_integ_TPBRW.o: Simpson_integ_TPBRW.h EF.h
EF.o: EF.h 

clean: 
	rm $(OBJS); \
	rm $(TARGET);

