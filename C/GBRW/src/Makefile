#
# Makefile for PATH_search_EM_Z_string program
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

OBJS =  main.o Optimiza_BFGS.o EF.o

LIBS =  Optimiza_BFGS.o EF.o

EXLIBS=-llbfgs -lm -lgc

TARGET =GBRW

.c.o:
	$(CC) $(INCDIR) -c $(CFLAG) $<;

all: GBRW install

# rules of generations
GBRW:  $(OBJS) $(LIBS)
	$(CC) $(CFLAG) -o $@ $(OBJS) $(LIBDIR) $(INCDIR) $(EXLIBS) ;

install: 
	cp $(TARGET) $(BIN) ; 

main.o: Optimiza_BFGS.h EF.h
Optimiza_BFGS.o: Optimiza_BFGS.h EF.h
EF.o: EF.h 

clean: 
	rm $(OBJS); \
	rm $(TARGET);
