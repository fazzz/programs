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

OBJS =  main_twoD_Gaussian.o EM_twoD.o EF.o Gaussian_twoD.o EM_reweight_twoD.o Simpson_integ_twoD.o

LIBS =  EM_twoD.h Simpson_integ_twoD.h EF.h Gaussian_twoD.h EM_reweight_twoD.h Simpson_integ_twoD.h

EXLIBS=-lm -lgc

TARGET =MGaussian_reweight_twoD

.c.o:
	$(CC) $(INCDIR) -c $(CFLAG) $<;

all: MGaussian_reweight_twoD install

# rules of generations
MGaussian_reweight_twoD:  $(OBJS) $(LIBS)
	$(CC) $(CFLAG) -o $@ $(OBJS) $(LIBDIR) $(INCDIR) $(EXLIBS) ;

install: 
	cp $(TARGET) $(BIN) ; 

main_twoD_Gaussian.o:	EM_reweight_twoD.h Simpson_integ_twoD.h EF.h Gaussian_twoD.h
EM_reweight_twoD.o: EM_reweight_twoD.h Simpson_integ_twoD.h EF.h
EM_twoD.o: EM_twoD.h EF.h
Gaussian_twoD.o: Gaussian_twoD.h EF.h
Simpson_integ_twoD.o: Simpson_integ_twoD.h EF.h
EF.o: EF.h 

clean: 
	rm $(OBJS); \
	rm $(TARGET);
