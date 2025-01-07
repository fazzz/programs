#
# Makefile for ABAMD program
#

HOME=${prefix}
SRC=.
BIN=${exec_prefix}/bin
LIBDIR=/home/yamamori/work/programs/ABAMD_2014_05_01/lib
INCDIR=/home/yamamori/work/programs/ABAMD_2014_05_01/include
netcdfdir=$(HOME)/lib/netcdf-4.3.0
CLAPACKdir=$(HOME)/lib/CLAPACK-3.2.1
gcdir=$(HOME)/lib/gc6.7

MAKE=make --no-print-directory

# definitions
#CC    = icc -pg -O3 -tpp6 -ipo
#CC    = icc -g -check=conversion,stack,uninit
CC    = icc -O3
#CC     = gcc

#CFLAG = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -m64 -O2
CFLAG = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -m64
#CFLAG =  #-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -m64 -O2

OBJS =  main.o ABAb.o ABA_mainpass.o ABA_backpass.o ABA_prepass.o \
	ABA_Nose-Hoover.o ABA_Nose-Hoover_new.o ABA_Nose-Hoover_new_mvV.o ABA_Nose-Hoover_chain.o  \
	ABA_set_trans.o ABA_set_tmat.o ABA_set_imat.o ABA_set_imatb.o ABA_set_rst.o ABA_set_lref.o  ABA_set_frc.o \
	ABA_update.o ABA_integ.o ABA_pick_data.o ABA_calcKineE.o  \
	EF.o PT.o PTL.o FF.o FFL.o SBFF.o MB.o BOXMULL.o TOPO.o LA.o quaternion.o netcdf_mineL.o mt19937ar.o \

LIBS =  ABAb.h ABA_backpass.h ABA_mainpass.h ABA_prepass.h \
	ABA_Nose-Hoover_new_mvV.h ABA_Nose-Hoover_new.h ABA_Nose-Hoover_chain.h ABA_Nose-Hoover.h \
	ABA_set_trans.h ABA_set_tmat.h ABA_set_rst.h ABA_set_lref.h ABA_set_imatb.h ABA_set_imat.h ABA_set_frc.h \
	ABA_update.h ABA_integ.h ABA_pick_data.h ABA_calcKineE.h \
	EF.h PT.h PTL.h FF.h FFL.h SBFF.h MB.h BOXMULL.h TOPO.h LA.h quaternion.h netcdf_mineL.h RAND.h

EXLIBS=-llapack -lblas -lf2c -lm -lgc -lnetcdf

TARGET =ABAMD_woGC_wcalloc

.c.o:
	$(CC) -I$(INCDIR) -c $(CFLAG) $<;

#all: netcdf clapack gc ABAMD
all: ABAMD

# rules of generations
ABAMD:  $(OBJS) $(LIBS)
#	$(CC) $(CFLAG) -o $@ $(OBJS) -L$(LIBDIR) -I$(INCDIR) $(EXLIBS) ;
	$(CC) $(CFLAG) -o $(TARGET) $(OBJS) -L$(LIBDIR) -I$(INCDIR) $(EXLIBS) ;

netcdf:
	cd $(netcdfdir); \
	./configure --prefix=$(ABAMDHOME) --enable-shared --disable-netcdf-4 --disable-dap ; \
	$(MAKE) ; \
	cp $(netcdfdir)/liblib/.libs/libnetcdf.a $(ABAMDHOME)/lib/; 
	cp $(netcdfdir)/include/netcdf.h $(ABAMDHOME)/include;

clapack:
	cd $(CLAPACKdir) ; \
	$(MAKE) ; \
	cp lapack_LINUX.a $(ABAMDHOME)/lib/liblapack.a ; \
	cp blas_LINUX.a $(ABAMDHOME)/lib/libblas.a ; \
	cp F2CLIBS/libf2c.a $(ABAMDHOME)/lib/ ; \
	cp INCLUDE/blaswrap.h INCLUDE/clapack.h INCLUDE/f2c.h  $(ABAMDHOME)/include ;

gc:
	cd $(gcdir) ; \
	./configure ; \
	$(MAKE) ; \
	cp $(gcdir)/.libs/libgc.a $(gcdir)/.libs/libgc.la $(gcdir)/.libs/libgc.lai $(gcdir)/.libs/libgc.so \
	$(gcdir)/.libs/libgc.so.1 $(gcdir)/.libs/libgc.so.1.0.2 $(ABAMDHOME)/lib ; \
	cp $(gcdir)/include/cord.h $(gcdir)/include/ec.h $(gcdir)/include/gc.h $(gcdir)/include/gc_alloc.h \
	   $(gcdir)/include/gc_allocator.h $(gcdir)/include/gc_amiga_redirects.h $(gcdir)/include/gc_backptr.h \
	   $(gcdir)/include/gc_config_macros.h $(gcdir)/include/gc_cpp.h $(gcdir)/include/gc_gcj.h \
	   $(gcdir)/include/gc_inl.h $(gcdir)/include/gc_inline.h $(gcdir)/include/gc_local_alloc.h \
	   $(gcdir)/include/gc_mark.h $(gcdir)/include/gc_pthread_redirects.h $(gcdir)/include/gc_typed.h \
	   $(gcdir)/include/javaxfc.h $(gcdir)/include/leak_detector.h $(gcdir)/include/new_gc_alloc.h \
	   $(gcdir)/include/weakpointer.h $(ABAMDHOME)/include ;

install: ABAMD \
	 cp @ ${exec_prefix}/bin ; 

main.o:                    ABAb.h PTL.h FFL.h EF.h netcdf_mineL.h EF.h
ABAb.o:                    ABAb.h EF.h
ABA_mainpass.o:            ABAb.h LA.h
ABA_backpass.o:            ABAb.h
ABA_prepass.o:             ABAb.h EF.h
ABA_Nose-Hoover.o:         ABAb.h
ABA_Nose-Hoover_new.o:     ABAb.h
ABA_Nose-Hoover_new_mvV.o: ABAb.h
ABA_Nose-Hoover_chain.o:   ABAb.h
ABA_set_trans.o:           ABAb.h EF.h LA.h
ABA_set_tmat.o:            ABAb.h EF.h
ABA_set_rst.o:             ABAb.h
ABA_set_lref.o:            ABAb.h EF.h
ABA_set_imat.o:            ABAb.h ABAb.h
ABA_set_imatb.o:           ABA_set_imatb.h EF.h
ABA_set_frc.o:             ABAb.h EF.h
ABA_update.o:              ABA_update.h quaternion.h EF.h
ABA_integ.o:               ABAb.h
ABA_pick_data.o:           ABAb.h
ABA_calcKineE.o:           ABAb.h LA.h EF.h RAND.h BOXMULL.h
EF.o:                      EF.h
PT.o:                      PT.h
PTL.o:                     PTL.h EF.h
FF.o:                      FF.h MB.h
FFL.o:                     FFL.h MB.h PTL.h LA.h TOPO.h EF.h
SBFF.o:                    SBFF.h
MB.o:                      MB.h
BOXMULL.o:                 EF.h RAND.h
TOPO.o:                    TOPO.h
LA.o:                      LA.h EF.h
quaternion.o:              quaternion.h
netcdf_mineL.o:            FF.h PTL.h SBFF.h netcdf_mineL.h
mt19937ar.o:               RAND.h

clean: ABAMDclean #netcdfclean clapackclean gcclean

ABAMDclean:
	rm $(OBJS); \
	rm $(TARGET);

netcdfclean:
	cd $(netcdfdir); \
	make clean ; \
	rm $(ABAMDHOME)/lib/libnetcdf.a;
	rm $(ABAMDHOME)/include/netcdf.h; 

clapackclean:
	cd $(CLAPACKdir) ; \
	make clean; \
	rm $(ABAMDHOME)/lib/libf2c.a $(ABAMDHOME)/lib/libblas.a $(ABAMDHOME)/lib/liblapack.a ; \
	rm $(ABAMDHOME)/include/blaswrap.h $(ABAMDHOME)/include/clapack.h $(ABAMDHOME)/include/f2c.h ;

gcclean:
	cd $(gcdir) ; \
	make clean ; \
	rm $(ABAMDHOME)/lib/libgc.a $(ABAMDHOME)/lib/libgc.la $(ABAMDHOME)/lib/libgc.lai $(ABAMDHOME)/lib/libgc.so \
	$(ABAMDHOME)/lib/libgc.so.1 $(ABAMDHOME)/lib/libgc.so.1.0.2 ; \
	rm $(ABAMDHOME)/include/cord.h $(ABAMDHOME)/include/ec.h $(ABAMDHOME)/include/gc.h $(ABAMDHOME)/include/gc_alloc.h \
	   $(ABAMDHOME)/include/gc_allocator.h $(ABAMDHOME)/include/gc_amiga_redirects.h $(ABAMDHOME)/include/gc_backptr.h \
	   $(ABAMDHOME)/include/gc_config_macros.h $(ABAMDHOME)/include/gc_cpp.h $(ABAMDHOME)/include/gc_gcj.h \
	   $(ABAMDHOME)/include/gc_inl.h $(ABAMDHOME)/include/gc_inline.h $(ABAMDHOME)/include/gc_local_alloc.h \
	   $(ABAMDHOME)/include/gc_mark.h $(ABAMDHOME)/include/gc_pthread_redirects.h $(ABAMDHOME)/include/gc_typed.h \
	   $(ABAMDHOME)/include/javaxfc.h $(ABAMDHOME)/include/leak_detector.h $(ABAMDHOME)/include/new_gc_alloc.h \
	   $(ABAMDHOME)/include/weakpointer.h  ;
