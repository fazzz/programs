#  Amber configuration file, created with: ./configure -mpi intel

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

FLIBS=  -lsff_mpi -lpbsa   -larpack -llapack -lblas  -L$(BASEDIR)/lib -lnetcdf  -L/opt/intel/Compiler/11.0/083/lib/intel64/ -lifport -lifcore -lsvml 
FLIBS_PTRAJ= -larpack -llapack -lblas   -L/opt/intel/Compiler/11.0/083/lib/intel64/ -lifport -lifcore -lsvml
FLIBSF= -larpack -llapack -lblas    -lsvml
FLIBS_FFTW3= 
###############################################################################

#  (3)  Modify any of the following if you need to change, e.g. to use gcc
#        rather than cc, etc.

SHELL=/bin/sh
INSTALLTYPE=parallel
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
CC=mpicc
CFLAGS= -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -DBINTRAJ -DMPI  $(CUSTOMBUILDFLAGS) 
CNOOPTFLAGS=
COPTFLAGS=-ip -O3 -xHost -DBINTRAJ -DHASGZ -DHASBZ2 
AMBERCFLAGS= $(AMBERBUILDFLAGS)

CXX=icpc
CPLUSPLUS=icpc
CXXFLAGS= -DMPI  $(CUSTOMBUILDFLAGS)
CXXNOOPTFLAGS=
CXXOPTFLAGS=-O3
AMBERCXXFLAGS= $(AMBERBUILDFLAGS)

NABFLAGS=
PBSAFLAG=

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

FC=mpif90
FFLAGS=  $(LOCALFLAGS) $(CUSTOMBUILDFLAGS) -I$(INCDIR) $(NETCDFINC) 
FNOOPTFLAGS= -O0
FOPTFLAGS= -ip -O3 -xHost
AMBERFFLAGS=$(AMBERBUILDFLAGS)
FREEFORMAT_FLAG= -FR
LM=-lm
FPP=cpp -traditional -P
FPPFLAGS= -DBINTRAJ -DMPI  $(CUSTOMBUILDFLAGS)
AMBERFPPFLAGS=$(AMBERBUILDFLAGS)
FCREAL8=

XHOME= /usr
XLIBS= -L/usr/lib64 -L/usr/lib
MAKE_XLEAP=install_xleap

NETCDF=$(BASEDIR)/include/netcdf.mod
NETCDF=$(BASEDIR)/include/netcdf.mod
NETCDFLIB=-L$(BASEDIR)/lib -lnetcdf
NETCDFINC=-I$(BASEDIR)/include
PNETCDF=
PNETCDFLIB=
FFTWLIB=

ZLIB=-lz
BZLIB=-lbz2

HASFC=yes
MTKPP=
XBLAS=
FFTW3=
MDGX=no

COMPILER=intel
MKL=
MKL_PROCESSOR=

#CUDA Specific build flags
NVCC=
PMEMD_CU_INCLUDES=
PMEMD_CU_LIBS=
PMEMD_CU_DEFINES=

#PMEMD Specific build flags
PMEMD_F90=mpif90 -DMPI   -DBINTRAJ -DDIRFRC_EFS -DDIRFRC_COMTRANS -DDIRFRC_NOVEC -DFFTLOADBAL_2PROC -DPUBFFT
PMEMD_FOPTFLAGS=#-ipo -O3 # -no-prec-div -xHost 2014-11-07 # 2014-11-07
PMEMD_CC=mpicc
PMEMD_COPTFLAGS=#-ipo -O3 -no-prec-div -xHost -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -DBINTRAJ -DMPI # 2014-11-07
PMEMD_FLIBSF= 
PMEMD_LD= mpif90 
LDOUT= -o 

#for NAB:
MPI=mpi

#1D-RISM
RISM=no

#3D-RISM NAB
RISMSFF=
SFF_RISM_INTERFACE=
TESTRISMSFF=

#3D-RISM SANDER
RISMSANDER=
SANDER_RISM_INTERFACE=
FLIBS_RISMSANDER=
TESTRISMSANDER=

#PUPIL
PUPILLIBS=-lrt -lm -lc -L${PUPIL_PATH}/lib -lPUPIL -lPUPILBlind

#Python interpreter we are using
PYTHON=/home/yamamori/bin/python2.7
