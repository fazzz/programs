#
# Makefile for PATH_search_EM_Z_string program
#

HOME=../
SRC=.
BIN=../bin
LIBDIR=-L/home/yamamori/work/programs/ABAMD_2014_05_01/lib #-L~/mylib -L~/lib
INCDIR=-I/home/yamamori/work/programs/ABAMD_2014_05_01/include #-I~/myinclude -I~/include

# definitions
#CC    = icc -pg -O3 -tpp6 -ipo
#CC    = icc -check=conversion,stack,uninit
#CC    = icc
CC    = gcc -g
#CC     = gcc

CFLAG = #-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -m64 -O2

OBJS =  main_oneD_Gaussian.o EM_reweight.o EM.o Simpson_integ.o EF.o #Gaussian.o

LIBS =  EM_reweight.h EM.h Simpson_integ.h EF.h #Gaussian.h

EXLIBS=-lm -lgc

TARGET =MGaussian_reweight

.c.o:
	$(CC) $(INCDIR) -c $(CFLAG) $<;

all: MGaussian_reweight install

# rules of generations
MGaussian_reweight:  $(OBJS) $(LIBS)
	$(CC) $(CFLAG) -o $@ $(OBJS) $(LIBDIR) $(INCDIR) $(EXLIBS) ;

install: 
	cp $(TARGET) $(BIN) ; 

main_oneD_Gaussian.o:	EM_reweight.h Simpson_integ.h EF.h
EM_reweight.o: EM_reweight.h Simpson_integ.h EF.h
EM.o: EM.h EF.h
#Gaussian.o: Gaussian.h EF.h
Simpson_integ.o: Simpson_integ.h EF.h
EF.o: EF.h 

clean: 
	rm $(OBJS); \
	rm $(TARGET);
