#
# Makefile for PATH_search_EM_Z_string program
#

HOME=../
SRC=.
BIN=../bin
LIBDIR=-L/home/yamamori/work/programs/SigmoidBaseReweight/lib -L/home/yamamori/work/programs/ABAMD_2014_05_01/lib
INCDIR=-I/home/yamamori/work/programs/SigmoidBaseReweight/include -I/home/yamamori/work/programs/ABAMD_2014_05_01/include

# definitions
#CC    = icc -pg -O3 -tpp6 -ipo
#CC    = icc -check=conversion,stack,uninit
#CC    = icc
CC    = gcc -g
#CC     = gcc

CFLAG = #-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -m64 -O2

OBJS =  main.o Optimiza_BFGS_SBRW.o Simpson_integ_SBRW.o EF.o

LIBS =  Optimiza_BFGS_SBRW.o EF.o

EXLIBS=-llbfgs -lm -lgc

TARGET =SBRW

.c.o:
	$(CC) $(INCDIR) -c $(CFLAG) $<;

all: SBRW install

# rules of generations
SBRW: $(OBJS) $(LIBS)
		$(CC) $(CFLAG) -o $@ $(OBJS) $(LIBDIR) $(INCDIR) $(EXLIBS) ;

install: 
	cp $(TARGET) $(BIN) ; 

main.o: Optimiza_BFGS_SBRW.h EF.h
Optimiza_BFGS_SBRW.o: Optimiza_BFGS_SBRW.h EF.h
Simpson_integ_SGRW.o: Simpson_integ_SGRW.h EF.h
EF.o: EF.h 

clean: 
	rm $(OBJS); \
	rm $(TARGET);
