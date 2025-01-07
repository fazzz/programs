#  Amber configuration file, created with: ./configure intel

###############################################################################

# (1)  Location of the installation

BASEDIR=/home/yamamori/software/Amber12/amber12
BINDIR=/home/yamamori/software/Amber12/amber12/bin
LIBDIR=/home/yamamori/software/Amber12/amber12/lib
INCDIR=/home/yamamori/software/Amber12/amber12/include
DATDIR=/home/yamamori/software/Amber12/amber12/dat
LOGDIR=/home/yamamori/software/Amber12/amber12/logs

###############################################################################


#  (2) If you want to search additional libraries by default, add them
#      to the FLIBS variable here.  (External libraries can also be linked into
#      NAB programs simply by including them on the command line; libraries
#      included in FLIBS are always searched.)

FLIBS=  -lsff -lpbsa -lrism -lfftw3 -larpack -llapack -lblas  -L$(BASEDIR)/lib -lnetcdf  -L/opt/intel/Compiler/11.0/083/lib/intel64/ -lifport -lifcore -lsvml 
FLIBS_PTRAJ= -larpack -llapack -lblas   -L/opt/intel/Compiler/11.0/083/lib/intel64/ -lifport -lifcore -lsvml
FLIBSF= -larpack -llapack -lblas    -lsvml
FLIBS_FFTW3= -lfftw3
###############################################################################

#  (3)  Modify any of the following if you need to change, e.g. to use gcc
#        rather than cc, etc.

SHELL=/bin/sh
INSTALLTYPE=serial
BUILDAMBER=amber

#  Set the C compiler, etc. 

#  The configure script should be fine, but if you need to hand-edit,
#  here is some info:

#  Example:  CC-->gcc; LEX-->flex; YACC-->yacc (built in byacc)
#     Note: If your lexer is "really" flex, you need to set
#     LEX=flex below.  For example, on some distributions,
#     /usr/bin/lex is really just a pointer to /usr/bin/flex,
#     so LEX=flex is necessary.  In general, gcc seems to need flex.

#   The compiler flags CFLAGS and CXXFLAGS should always be used.
#   By contrast, *OPTFLAGS and *NOOPTFLAGS will only be used with
#   certain files, and usually at compile-time but not link-time.
#   Where *OPTFLAGS and *NOOPTFLAGS are requested (in Makefiles,
#   makedepend and depend), they should come before CFLAGS or
#   CXXFLAGS; this allows the user to override *OPTFLAGS and
#   *NOOPTFLAGS using the BUILDFLAGS variable.
#
CC=icc
CFLAGS= -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -DBINTRAJ  $(CUSTOMBUILDFLAGS) 
CNOOPTFLAGS=
COPTFLAGS=-ip -O3 -xHost -DBINTRAJ -DHASGZ -DHASBZ2 
AMBERCFLAGS= $(AMBERBUILDFLAGS)

CXX=icpc
CPLUSPLUS=icpc
CXXFLAGS=  $(CUSTOMBUILDFLAGS)
CXXNOOPTFLAGS=
CXXOPTFLAGS=-O3
AMBERCXXFLAGS= $(AMBERBUILDFLAGS)

NABFLAGS=
PBSAFLAG=-DFFTW

LDFLAGS=-shared-intel  $(CUSTOMBUILDFLAGS)
AMBERLDFLAGS=$(AMBERBUILDFLAGS)

LEX=   flex
YACC=  $(BINDIR)/yacc
AR=    ar rv
M4=    m4
RANLIB=ranlib

#  Set the C-preprocessor.  Code for a small preprocessor is in
#    ucpp-1.3;  it gets installed as $(BINDIR)/ucpp;
#    this can generally be used (maybe not on 64-bit machines like altix).

CPP=    ucpp -l

#  These variables control whether we will use compiled versions of BLAS
#  and LAPACK (which are generally slower), or whether those libraries are
#  already available (presumably in an optimized form).

LAPACK=install
BLAS=install
F2C=skip

#  These variables determine whether builtin versions of certain components
#  can be used, or whether we need to compile our own versions.

UCPP=install
C9XCOMPLEX=skip

#  For Windows/cygwin, set SFX to ".exe"; for Unix/Linux leave it empty:
#  Set OBJSFX to ".obj" instead of ".o" on Windows:

SFX=
OSFX=.o
MV=mv
RM=rm
CP=cp

#  Information about Fortran compilation:

FC=ifort
FFLAGS=  $(LOCALFLAGS) $(CUSTOMBUILDFLAGS) -I$(INCDIR) $(NETCDFINC) 
FNOOPTFLAGS= -O0
FOPTFLAGS= -ip -O3 -xHost
AMBERFFLAGS=$(AMBERBUILDFLAGS)
FREEFORMAT_FLAG= -FR
LM=-lm
FPP=cpp -traditional -P
FPPFLAGS= -DBINTRAJ  $(CUSTOMBUILDFLAGS)
AMBERFPPFLAGS=$(AMBERBUILDFLAGS)
FCREAL8=

XHOME= /usr
XLIBS= -L/usr/lib64 -L/usr/lib
MAKE_XLEAP=install_xleap

NETCDF=$(BASEDIR)/include/netcdf.mod
NETCDFLIB=-L$(BASEDIR)/lib -lnetcdf
NETCDFINC=-I$(BASEDIR)/include
PNETCDF=
PNETCDFLIB=
FFTWLIB=-lfftw3

ZLIB=-lz
BZLIB=-lbz2

HASFC=yes
MTKPP=install_mtkpp
XBLAS=$(LIBDIR)/libxblas-amb.a
FFTW3=$(LIBDIR)/libfftw3.a
MDGX=yes

COMPILER=intel
MKL=
MKL_PROCESSOR=

#CUDA Specific build flags
NVCC=
PMEMD_CU_INCLUDES=
PMEMD_CU_LIBS=
PMEMD_CU_DEFINES=

#PMEMD Specific build flags
PMEMD_F90=ifort   -DBINTRAJ -DDIRFRC_EFS -DDIRFRC_COMTRANS -DDIRFRC_NOVEC -DFFTLOADBAL_2PROC -DPUBFFT
PMEMD_FOPTFLAGS=-ipo -O3 -no-prec-div -xHost
PMEMD_CC=icc
PMEMD_COPTFLAGS=-ipo -O3 -no-prec-div -xHost -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -DBINTRAJ 
PMEMD_FLIBSF= 
PMEMD_LD= ifort 
LDOUT= -o 

#for NAB:
MPI=

#1D-RISM
RISM=yes

#3D-RISM NAB
RISMSFF=-DRISMSFF
SFF_RISM_INTERFACE=../rism/amber_rism_interface.NAB.o
TESTRISMSFF=testrism

#3D-RISM SANDER
RISMSANDER=
SANDER_RISM_INTERFACE=
FLIBS_RISMSANDER=
TESTRISMSANDER=

#PUPIL
PUPILLIBS=-lrt -lm -lc -L${PUPIL_PATH}/lib -lPUPIL -lPUPILBlind

#Python interpreter we are using
PYTHON=/home/yamamori/bin/python2.7
