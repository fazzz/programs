#
# Makefile for MD_NHC_MP1996
#

HOME=/home/yamamori/work/programs/MD_pep
SRC=.
BIN=../bin

LIBDIR=-L../lib/netcdf -L../lib/gc
INCDIR=-I../include/netcdf -I../include/gc

netcdfdir=$(HOME)/lib/netcdf-4.3.0
gcdir=$(HOME)/lib/gc-7.0

MAKE=make --no-print-directory

# definitions
#CC    = icc -pg -O3 -tpp6 -ipo
#CC    = icc -check=conversion,stack,uninit
#CC    = icc
CC    = gcc -g
#CC     = gcc

CFLAG = #-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -m64 -O2

OBJS =  MD_pep_NH_MP1996_AAFF_Amber.o MD_NHC_MP1996.o \
	UMBP.o PTL.o FFL.o TOPO.o LA.o MB.o EF.o \
	netcdf_mineL.o BOXMULL.o mt19937ar.o 

LIBS =  MD_NHC_MP1996.h UMBP.h EF.h PTL.h FFL.h TOPO.h LA.h MB.h netcdf_mineL.h BOXMULL.h RAND.h

EXLIBS=-lm -lgc -lnetcdf

TARGET = MD_NHC_MP1996

.c.o:
	$(CC) $(INCDIR) -c $(CFLAG) $<;

all: MD_NHC_MP1996 install

# rules of generations
MD_NHC_MP1996:  $(OBJS) $(LIBS)
	$(CC) $(CFLAG) -o $@ $(OBJS) $(LIBDIR) $(INCDIR) $(EXLIBS) ;

netcdf:
	cd $(netcdfdir); \
	./configure --prefix=$(HOME) --disable-netcdf-4 ; \
	$(MAKE) ; \
	cp $(netcdfdir)/liblib/.libs/libnetcdf.a ../../lib/netcdf; 
	cp $(netcdfdir)/include/netcdf.h $(HOME)/include/netcdf;

gc:
	cd $(gcdir) ; \
	./configure ; \
	$(MAKE) ; \
	cp $(gcdir)/.libs/libgc.a $(gcdir)/.libs/libgc.la $(gcdir)/.libs/libgc.lai $(gcdir)/.libs/libgc.so \
	$(gcdir)/.libs/libgc.so.1 $(gcdir)/.libs/libgc.so.1.0.2 ../../lib/gc ; \
	cp $(gcdir)/include/cord.h $(gcdir)/include/ec.h $(gcdir)/include/gc.h $(gcdir)/include/gc_alloc.h \
	   $(gcdir)/include/gc_allocator.h $(gcdir)/include/gc_amiga_redirects.h $(gcdir)/include/gc_backptr.h \
	   $(gcdir)/include/gc_config_macros.h $(gcdir)/include/gc_cpp.h $(gcdir)/include/gc_gcj.h \
	   $(gcdir)/include/gc_inl.h $(gcdir)/include/gc_inline.h $(gcdir)/include/gc_local_alloc.h \
	   $(gcdir)/include/gc_mark.h $(gcdir)/include/gc_pthread_redirects.h $(gcdir)/include/gc_typed.h \
	   $(gcdir)/include/javaxfc.h $(gcdir)/include/leak_detector.h $(gcdir)/include/new_gc_alloc.h \
	   $(gcdir)/include/weakpointer.h $(HOME)/include/gc ;

install: 
	cp $(TARGET) $(BIN) ; 

MD_pep_NH_MP1996_AAFF_Amber.o: MD_NHC_MP1996.h UMBP.h PTL.h netcdf_mineL.h BOXMULL.h EF.h
MD_NHC_MP1996.o:               MD_NHC_MP1996.h FFL.h UMBP.h EF.h
UMBP.o:                        UMBP.h PTL.h TOPO.h LA.h
PTL.o:                         PTL.h EF.h
FFL.o:                         FFL.h MB.h PTL.h LA.h TOPO.h EF.h
TOPO.o:                        TOPO.h
MB.o:                          MB.h
LA.o:                          LA.h
netcdf_mineL.o:                netcdf_mineL.h
BOXMULL.o:                     BOXMULL.h RAND.h EF.h
mt19937ar.o:                   RAND.h
EF.o:                          EF.h

clean: 
	rm $(OBJS); \
	rm $(TARGET);
	rm ../bin/$(TARGET);

netcdfclean:
	cd $(netcdfdir); \
	make clean ; \
	rm $(HOME)/lib/libnetcdf.a;
	rm $(HOME)/include/netcdf.h; 

gcclean:
	cd $(gcdir) ; \
	make clean ; \
	rm $(HOME)/lib/gc; \
	rm $(HOME)/include/gc
