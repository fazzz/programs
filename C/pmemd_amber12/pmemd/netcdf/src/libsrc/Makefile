# Makefile.in generated by automake 1.9.6 from Makefile.am.
# libsrc/Makefile.  Generated from Makefile.in by configure.

# Copyright (C) 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002,
# 2003, 2004, 2005  Free Software Foundation, Inc.
# This Makefile.in is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.



# Copyright 2005, see the COPYRIGHT file for more information.

# THis automake file is in charge of building the libsrc directory.

# $Id: Makefile.in,v 10.0 2008/04/15 23:23:12 case Exp $


srcdir = .
top_srcdir = ..

pkgdatadir = $(datadir)/netcdf
pkglibdir = $(libdir)/netcdf
pkgincludedir = $(includedir)/netcdf
top_builddir = ..
am__cd = CDPATH="$${ZSH_VERSION+.}$(PATH_SEPARATOR)" && cd
INSTALL = /usr/bin/install -c
install_sh_DATA = $(install_sh) -c -m 644
install_sh_PROGRAM = $(install_sh) -c
install_sh_SCRIPT = $(install_sh) -c
INSTALL_HEADER = $(INSTALL_DATA)
transform = $(program_transform_name)
NORMAL_INSTALL = :
PRE_INSTALL = :
POST_INSTALL = :
NORMAL_UNINSTALL = :
PRE_UNINSTALL = :
POST_UNINSTALL = :
build_triplet = x86_64-unknown-linux-gnu
host_triplet = x86_64-unknown-linux-gnu

# Does the user want to use ffio or posixio?
#am__append_1 = ffio.c
am__append_2 = posixio.c

# Does  the user want the V2 API?

# This will create a convenience library for all the netcdf-2
# functions.
am__append_3 = libnetcdf2.la 

# If the user enabled netcdf-4 in the configure, then the netcdf-3
# library will only be built as a non-installing library.

# If the user specified a root for HDF5, use it.
##am__append_4 = -I@HDF5DIR@/include

# This will create a convenience library for all the netcdf-3
# functions.
#am__append_5 = libnetcdf3.la

# Include netcdf-2 convenience library.
am__append_6 = libnetcdf2.la

# If we are not building separate fortran libraries, then include
# their functions in the C library.

# If F77 is built, include its convenience library. If F90 is built,
# it will already be part of the F77 convenience library.
am__append_7 = ../fortran/libnetcdff.la

# Verson for the netcdf-3 library. Don't mess with this number
# lightly!! This is not the same as the netCDF version, it is a
# completely different versioning system for shared library binary
# files, which is used in a specific way by OSs which are using shared
# libraries. Do don't dork around with it or you will break netCDF for
# shared library users, and everyone will be very confused. We are
# starting with 4 because some package distributors (and other users)
# have already build netCDF with shared libraries, and we need to have
# a larger number than any they used.
am__append_8 = $(AM_LDFLAGS) -version-info 4:0:0 -no-undefined

# If building the DLL on mingw, do some extra stuff. For now, this
# works with netCDF-3 only. (Ed 12/3/6)

# We need extra arguments to produce the def file, which is needed by
# MS tools to create the VC++ import library. Add -avoid-version to
# avoid the version number in the .dll file name. 
#am__append_9 = -Wl,--output-def,.libs/libnetcdf.def 
#am__append_10 = -DNC_DLL_EXPORT
check_PROGRAMS = $(am__EXEEXT_1)
subdir = libsrc
DIST_COMMON = $(am__include_HEADERS_DIST) $(srcdir)/Makefile.am \
	$(srcdir)/Makefile.in
ACLOCAL_M4 = $(top_srcdir)/aclocal.m4
am__aclocal_m4_deps = $(top_srcdir)/acinclude.m4 \
	$(top_srcdir)/configure.ac
am__configure_deps = $(am__aclocal_m4_deps) $(CONFIGURE_DEPENDENCIES) \
	$(ACLOCAL_M4)
mkinstalldirs = $(install_sh) -d
CONFIG_HEADER = $(top_builddir)/config.h \
	$(top_builddir)/fortran/nfconfig1.inc
CONFIG_CLEAN_FILES =
am__vpath_adj_setup = srcdirstrip=`echo "$(srcdir)" | sed 's|.|.|g'`;
am__vpath_adj = case $$p in \
    $(srcdir)/*) f=`echo "$$p" | sed "s|^$$srcdirstrip/||"`;; \
    *) f=$$p;; \
  esac;
am__strip_dir = `echo $$p | sed -e 's|^.*/||'`;
am__installdirs = "$(DESTDIR)$(libdir)" "$(DESTDIR)$(man3dir)" \
	"$(DESTDIR)$(includedir)"
libLTLIBRARIES_INSTALL = $(INSTALL)
LTLIBRARIES = $(lib_LTLIBRARIES) $(noinst_LTLIBRARIES)
am__DEPENDENCIES_1 = libnetcdf2.la
am__DEPENDENCIES_2 = ../fortran/libnetcdff.la
libnetcdf_la_DEPENDENCIES = $(am__DEPENDENCIES_1) \
	$(am__DEPENDENCIES_2)
am__libnetcdf_la_SOURCES_DIST = attr.m4 ncx.m4 putget.m4 dim.c error.c \
	libvers.c nc.c string.c v1hpg.c var.c fbits.h netcdf3.h nc.h \
	ncio.h ncx.h netcdf.h onstack.h rnd.h ffio.c posixio.c
#am__objects_1 = ffio.lo
am__objects_2 = posixio.lo
am__objects_3 = attr.lo ncx.lo putget.lo dim.lo error.lo libvers.lo \
	nc.lo string.lo v1hpg.lo var.lo $(am__objects_1) \
	$(am__objects_2)
am_libnetcdf_la_OBJECTS = $(am__objects_3)
libnetcdf_la_OBJECTS = $(am_libnetcdf_la_OBJECTS)
am_libnetcdf_la_rpath = -rpath $(libdir)
libnetcdf2_la_LIBADD =
am__libnetcdf2_la_SOURCES_DIST = v2i.c
am_libnetcdf2_la_OBJECTS = v2i.lo
libnetcdf2_la_OBJECTS = $(am_libnetcdf2_la_OBJECTS)
am_libnetcdf2_la_rpath =
libnetcdf3_la_LIBADD =
am__libnetcdf3_la_SOURCES_DIST = attr.m4 ncx.m4 putget.m4 dim.c \
	error.c libvers.c nc.c string.c v1hpg.c var.c fbits.h \
	netcdf3.h nc.h ncio.h ncx.h netcdf.h onstack.h rnd.h ffio.c \
	posixio.c nc3convert.h
#am_libnetcdf3_la_OBJECTS = $(am__objects_3)
libnetcdf3_la_OBJECTS = $(am_libnetcdf3_la_OBJECTS)
#am_libnetcdf3_la_rpath =
am__EXEEXT_1 = t_nc$(EXEEXT)
t_nc_SOURCES = t_nc.c
t_nc_OBJECTS = t_nc-t_nc.$(OBJEXT)
t_nc_LDADD = $(LDADD)
am__DEPENDENCIES_3 = libnetcdf.la
t_nc_DEPENDENCIES = $(am__DEPENDENCIES_3)
#t_nc_DEPENDENCIES = libnetcdf3.la
DEFAULT_INCLUDES = -I. -I$(srcdir) -I$(top_builddir) -I$(top_builddir)/fortran
depcomp = $(SHELL) $(top_srcdir)/depcomp
am__depfiles_maybe = depfiles
COMPILE = $(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) \
	$(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS)
LTCOMPILE = $(LIBTOOL) --tag=CC --mode=compile $(CC) $(DEFS) \
	$(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
	$(AM_CFLAGS) $(CFLAGS)
CCLD = $(CC)
LINK = $(LIBTOOL) --tag=CC --mode=link $(CCLD) $(AM_CFLAGS) $(CFLAGS) \
	$(AM_LDFLAGS) $(LDFLAGS) -o $@
SOURCES = $(libnetcdf_la_SOURCES) $(libnetcdf2_la_SOURCES) \
	$(libnetcdf3_la_SOURCES) t_nc.c
DIST_SOURCES = $(am__libnetcdf_la_SOURCES_DIST) \
	$(am__libnetcdf2_la_SOURCES_DIST) \
	$(am__libnetcdf3_la_SOURCES_DIST) t_nc.c
man3dir = $(mandir)/man3
NROFF = nroff
MANS = $(man_MANS)
am__include_HEADERS_DIST = netcdf.h
includeHEADERS_INSTALL = $(INSTALL_HEADER)
HEADERS = $(include_HEADERS)
ETAGS = etags
CTAGS = ctags
DISTFILES = $(DIST_COMMON) $(DIST_SOURCES) $(TEXINFOS) $(EXTRA_DIST)
ACLOCAL = ${SHELL} /home/yamamori/software/Amber12/amber12/AmberTools/src/netcdf/src/missing --run aclocal-1.9
ALLOCA = 
AMDEP_FALSE = #
AMDEP_TRUE = 
AMTAR = ${SHELL} /home/yamamori/software/Amber12/amber12/AmberTools/src/netcdf/src/missing --run tar
AR = ar
AS = as
AUTOCONF = ${SHELL} /home/yamamori/software/Amber12/amber12/AmberTools/src/netcdf/src/missing --run autoconf
AUTOHEADER = ${SHELL} /home/yamamori/software/Amber12/amber12/AmberTools/src/netcdf/src/missing --run autoheader
AUTOMAKE = ${SHELL} /home/yamamori/software/Amber12/amber12/AmberTools/src/netcdf/src/missing --run automake-1.9
AWK = gawk
BUILD_CXX_FALSE = 
BUILD_CXX_TRUE = #
BUILD_C_FALSE = #
BUILD_C_TRUE = 
BUILD_DLL_FALSE = 
BUILD_DLL_TRUE = #
BUILD_EXAMPLES_FALSE = #
BUILD_EXAMPLES_TRUE = 
BUILD_F77_FALSE = #
BUILD_F77_TRUE = 
BUILD_F90_FALSE = #
BUILD_F90_TRUE = 
BUILD_SEPARATE_FORTRAN_FALSE = 
BUILD_SEPARATE_FORTRAN_TRUE = #
BUILD_SHARED_FALSE = 
BUILD_SHARED_TRUE = #
BUILD_UTILITIES_FALSE = #
BUILD_UTILITIES_TRUE = 
BUILD_V2_FALSE = #
BUILD_V2_TRUE = 
CC = icc
CCDEPMODE = depmode=gcc3
CFLAGS =  -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE 
CPP = icc -E
CPPFLAGS = -DNDEBUG
CXX = icpc
CXXCPP = icpc -E
CXXDEPMODE = depmode=gcc3
CXXFLAGS =  
CYGPATH_W = echo
DEFS = -DHAVE_CONFIG_H
DEPDIR = .deps
DLLTOOL = dlltool
DSYMUTIL = 
ECHO = echo
ECHO_C = 
ECHO_N = -n
ECHO_T = 
EGREP = /bin/grep -E
EXEEXT = 
F77 = ifort
FC = ifort
FCFLAGS = -shared-intel  -ip -O3 -xHost
FCFLAGS_f90 = 
FFLAGS = -shared-intel  -ip -O3 -xHost
FTPBINDIR = /home/ftp/pub/binary/dummy_system
GREP = /bin/grep
INSTALL_DATA = ${INSTALL} -m 644
INSTALL_DOCS_FALSE = 
INSTALL_DOCS_TRUE = #
INSTALL_PROGRAM = ${INSTALL}
INSTALL_SCRIPT = ${INSTALL}
INSTALL_STRIP_PROGRAM = ${SHELL} $(install_sh) -c -s
LARGE_FILE_TESTS_FALSE = 
LARGE_FILE_TESTS_TRUE = #
LDFLAGS = -shared-intel 
LIBOBJS = 
LIBS = 
LIBTOOL = $(SHELL) $(top_builddir)/libtool
LN_S = ln -s
LTLIBOBJS = 
MAINT = #
MAINTAINER_MODE_FALSE = 
MAINTAINER_MODE_TRUE = #
MAKEINFO = ${SHELL} /home/yamamori/software/Amber12/amber12/AmberTools/src/netcdf/src/missing --run makeinfo
MAKEWHATIS_CMD = 
MOD_FLAG = -I
NETCDF_C_LIB = 
NM = /usr/bin/nm -B
NMEDIT = 
NMFLAGS = 
OBJDUMP = objdump
OBJEXT = o
PACKAGE = netcdf
PACKAGE_BUGREPORT = support@unidata.ucar.edu
PACKAGE_NAME = netCDF
PACKAGE_STRING = netCDF 3.6.2
PACKAGE_TARNAME = netcdf
PACKAGE_VERSION = 3.6.2
PATH_SEPARATOR = :
RANLIB = ranlib
SED = /bin/sed
SET_MAKE = 
SHAREEXT_DYLIB_FALSE = 
SHAREEXT_DYLIB_TRUE = #
SHAREEXT_SL_FALSE = 
SHAREEXT_SL_TRUE = #
SHELL = /bin/sh
STRIP = strip
TEST_F77_FALSE = #
TEST_F77_TRUE = 
TEST_PARALLEL_FALSE = 
TEST_PARALLEL_TRUE = #
UPPER_CASE_MOD_FALSE = 
UPPER_CASE_MOD_TRUE = #
USE_FFIO_FALSE = 
USE_FFIO_TRUE = #
USE_HDF5_DIR_FALSE = 
USE_HDF5_DIR_TRUE = #
USE_NETCDF4_FALSE = 
USE_NETCDF4_TRUE = #
USE_NETCDF_C_LIB_DIR_FALSE = 
USE_NETCDF_C_LIB_DIR_TRUE = #
VERSION = 3.6.2
WHATIS = whatis
ac_ct_CC = 
ac_ct_CXX = 
ac_ct_F77 = 
ac_ct_FC = 
am__fastdepCC_FALSE = #
am__fastdepCC_TRUE = 
am__fastdepCXX_FALSE = #
am__fastdepCXX_TRUE = 
am__include = include
am__leading_dot = .
am__quote = 
am__tar = ${AMTAR} chof - "$$tardir"
am__untar = ${AMTAR} xf -
bindir = ${exec_prefix}/bin
build = x86_64-unknown-linux-gnu
build_alias = 
build_cpu = x86_64
build_os = linux-gnu
build_vendor = unknown
datadir = ${datarootdir}
datarootdir = ${prefix}/share
docdir = ${datarootdir}/doc/${PACKAGE_TARNAME}
dvidir = ${docdir}
exec_prefix = ${prefix}
host = x86_64-unknown-linux-gnu
host_alias = 
host_cpu = x86_64
host_os = linux-gnu
host_vendor = unknown
htmldir = ${docdir}
includedir = ${prefix}/include
infodir = ${datarootdir}/info
install_sh = /home/yamamori/software/Amber12/amber12/AmberTools/src/netcdf/src/install-sh
libdir = ${exec_prefix}/lib
libexecdir = ${exec_prefix}/libexec
localedir = ${datarootdir}/locale
localstatedir = ${prefix}/var
mandir = ${datarootdir}/man
mkdir_p = mkdir -p --
oldincludedir = /usr/include
pdfdir = ${docdir}
prefix = /home/yamamori/software/Amber12/amber12/AmberTools
prog = 
program_transform_name = s,x,x,
psdir = ${docdir}
sbindir = ${exec_prefix}/sbin
sharedstatedir = ${prefix}/com
sysconfdir = ${prefix}/etc
target_alias = 

# The C API man page.
man_MANS = netcdf.3 

# These files are part of the netCDF-3 library.
NC3_SOURCES = attr.m4 ncx.m4 putget.m4 dim.c error.c libvers.c nc.c \
	string.c v1hpg.c var.c fbits.h netcdf3.h nc.h ncio.h ncx.h \
	netcdf.h onstack.h rnd.h $(am__append_1) $(am__append_2)

# We may have to add to these later.
libnetcdf_la_LIBADD = $(am__append_6) $(am__append_7)
noinst_LTLIBRARIES = $(am__append_3) $(am__append_5)
AM_CPPFLAGS = $(am__append_4) $(am__append_10)
libnetcdf_la_LDFLAGS = $(am__append_8) $(am__append_9)

# This is the v2 API source code.
libnetcdf2_la_SOURCES = v2i.c
##AM_LDFLAGS = -L@HDF5DIR@/lib

# This is a list of the code netcdf-3 library.
#libnetcdf3_la_SOURCES = $(NC3_SOURCES) nc3convert.h

# Compile in the library when building the test.
LDADD = ${lib_LTLIBRARIES}

# Tell the test program to link to netcdf 3 library.
#LDADD = libnetcdf3.la

# This will cause the netcdf-3 header file to be installed.
include_HEADERS = netcdf.h

# This will cause the netcdf-3 library to be built.
lib_LTLIBRARIES = libnetcdf.la

# This is a list of the code netcdf-3 library.
libnetcdf_la_SOURCES = $(NC3_SOURCES)

# Test the netCDF-3 library.
TESTPROGRAMS = t_nc
TESTS = $(TESTPROGRAMS)

# If building a DLL, the test program needs to define DLL_NETCDF.
#t_nc_CPPFLAGS = -DDLL_NETCDF

# These files are cleaned on developer workstations (and then rebuilt
# with m4), but they are included in the distribution so that the user
# does not have to have m4.
MAINTAINERCLEANFILES = attr.c ncx.c putget.c Makefile.in $(man_MANS)
EXTRA_DIST = attr.c ncx.c putget.c $(man_MANS)

# This file is created by the test.
CLEANFILES = test.nc
all: all-am

.SUFFIXES:
.SUFFIXES: .c .lo .m4 .o .obj
$(srcdir)/Makefile.in: # $(srcdir)/Makefile.am  $(am__configure_deps)
	@for dep in $?; do \
	  case '$(am__configure_deps)' in \
	    *$$dep*) \
	      cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh \
		&& exit 0; \
	      exit 1;; \
	  esac; \
	done; \
	echo ' cd $(top_srcdir) && $(AUTOMAKE) --foreign  libsrc/Makefile'; \
	cd $(top_srcdir) && \
	  $(AUTOMAKE) --foreign  libsrc/Makefile
.PRECIOUS: Makefile
Makefile: $(srcdir)/Makefile.in $(top_builddir)/config.status
	@case '$?' in \
	  *config.status*) \
	    cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh;; \
	  *) \
	    echo ' cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@ $(am__depfiles_maybe)'; \
	    cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@ $(am__depfiles_maybe);; \
	esac;

$(top_builddir)/config.status: $(top_srcdir)/configure $(CONFIG_STATUS_DEPENDENCIES)
	cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh

$(top_srcdir)/configure: # $(am__configure_deps)
	cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh
$(ACLOCAL_M4): # $(am__aclocal_m4_deps)
	cd $(top_builddir) && $(MAKE) $(AM_MAKEFLAGS) am--refresh
install-libLTLIBRARIES: $(lib_LTLIBRARIES)
	@$(NORMAL_INSTALL)
	test -z "$(libdir)" || $(mkdir_p) "$(DESTDIR)$(libdir)"
	@list='$(lib_LTLIBRARIES)'; for p in $$list; do \
	  if test -f $$p; then \
	    f=$(am__strip_dir) \
	    echo " $(LIBTOOL) --mode=install $(libLTLIBRARIES_INSTALL) $(INSTALL_STRIP_FLAG) '$$p' '$(DESTDIR)$(libdir)/$$f'"; \
	    $(LIBTOOL) --mode=install $(libLTLIBRARIES_INSTALL) $(INSTALL_STRIP_FLAG) "$$p" "$(DESTDIR)$(libdir)/$$f"; \
	  else :; fi; \
	done

uninstall-libLTLIBRARIES:
	@$(NORMAL_UNINSTALL)
	@set -x; list='$(lib_LTLIBRARIES)'; for p in $$list; do \
	  p=$(am__strip_dir) \
	  echo " $(LIBTOOL) --mode=uninstall rm -f '$(DESTDIR)$(libdir)/$$p'"; \
	  $(LIBTOOL) --mode=uninstall rm -f "$(DESTDIR)$(libdir)/$$p"; \
	done

clean-libLTLIBRARIES:
	-test -z "$(lib_LTLIBRARIES)" || rm -f $(lib_LTLIBRARIES)
	@list='$(lib_LTLIBRARIES)'; for p in $$list; do \
	  dir="`echo $$p | sed -e 's|/[^/]*$$||'`"; \
	  test "$$dir" != "$$p" || dir=.; \
	  echo "rm -f \"$${dir}/so_locations\""; \
	  rm -f "$${dir}/so_locations"; \
	done

clean-noinstLTLIBRARIES:
	-test -z "$(noinst_LTLIBRARIES)" || rm -f $(noinst_LTLIBRARIES)
	@list='$(noinst_LTLIBRARIES)'; for p in $$list; do \
	  dir="`echo $$p | sed -e 's|/[^/]*$$||'`"; \
	  test "$$dir" != "$$p" || dir=.; \
	  echo "rm -f \"$${dir}/so_locations\""; \
	  rm -f "$${dir}/so_locations"; \
	done
libnetcdf.la: $(libnetcdf_la_OBJECTS) $(libnetcdf_la_DEPENDENCIES) 
	$(LINK) $(am_libnetcdf_la_rpath) $(libnetcdf_la_LDFLAGS) $(libnetcdf_la_OBJECTS) $(libnetcdf_la_LIBADD) $(LIBS)
libnetcdf2.la: $(libnetcdf2_la_OBJECTS) $(libnetcdf2_la_DEPENDENCIES) 
	$(LINK) $(am_libnetcdf2_la_rpath) $(libnetcdf2_la_LDFLAGS) $(libnetcdf2_la_OBJECTS) $(libnetcdf2_la_LIBADD) $(LIBS)
libnetcdf3.la: $(libnetcdf3_la_OBJECTS) $(libnetcdf3_la_DEPENDENCIES) 
	$(LINK) $(am_libnetcdf3_la_rpath) $(libnetcdf3_la_LDFLAGS) $(libnetcdf3_la_OBJECTS) $(libnetcdf3_la_LIBADD) $(LIBS)

clean-checkPROGRAMS:
	@list='$(check_PROGRAMS)'; for p in $$list; do \
	  f=`echo $$p|sed 's/$(EXEEXT)$$//'`; \
	  echo " rm -f $$p $$f"; \
	  rm -f $$p $$f ; \
	done
t_nc$(EXEEXT): $(t_nc_OBJECTS) $(t_nc_DEPENDENCIES) 
	@rm -f t_nc$(EXEEXT)
	$(LINK) $(t_nc_LDFLAGS) $(t_nc_OBJECTS) $(t_nc_LDADD) $(LIBS)

mostlyclean-compile:
	-rm -f *.$(OBJEXT)

distclean-compile:
	-rm -f *.tab.c

include ./$(DEPDIR)/attr.Plo
include ./$(DEPDIR)/dim.Plo
include ./$(DEPDIR)/error.Plo
include ./$(DEPDIR)/ffio.Plo
include ./$(DEPDIR)/libvers.Plo
include ./$(DEPDIR)/nc.Plo
include ./$(DEPDIR)/ncx.Plo
include ./$(DEPDIR)/posixio.Plo
include ./$(DEPDIR)/putget.Plo
include ./$(DEPDIR)/string.Plo
include ./$(DEPDIR)/t_nc-t_nc.Po
include ./$(DEPDIR)/v1hpg.Plo
include ./$(DEPDIR)/v2i.Plo
include ./$(DEPDIR)/var.Plo

.c.o:
	depbase=`echo $@ | sed 's|[^/]*$$|$(DEPDIR)/&|;s|\.o$$||'`; \
	if $(COMPILE) -MT $@ -MD -MP -MF "$$depbase.Tpo" -c -o $@ $<; \
	then mv -f "$$depbase.Tpo" "$$depbase.Po"; else rm -f "$$depbase.Tpo"; exit 1; fi
#	source='$<' object='$@' libtool=no \
#	DEPDIR=$(DEPDIR) $(CCDEPMODE) $(depcomp) \
#	$(COMPILE) -c -o $@ $<

.c.obj:
	depbase=`echo $@ | sed 's|[^/]*$$|$(DEPDIR)/&|;s|\.obj$$||'`; \
	if $(COMPILE) -MT $@ -MD -MP -MF "$$depbase.Tpo" -c -o $@ `$(CYGPATH_W) '$<'`; \
	then mv -f "$$depbase.Tpo" "$$depbase.Po"; else rm -f "$$depbase.Tpo"; exit 1; fi
#	source='$<' object='$@' libtool=no \
#	DEPDIR=$(DEPDIR) $(CCDEPMODE) $(depcomp) \
#	$(COMPILE) -c -o $@ `$(CYGPATH_W) '$<'`

.c.lo:
	depbase=`echo $@ | sed 's|[^/]*$$|$(DEPDIR)/&|;s|\.lo$$||'`; \
	if $(LTCOMPILE) -MT $@ -MD -MP -MF "$$depbase.Tpo" -c -o $@ $<; \
	then mv -f "$$depbase.Tpo" "$$depbase.Plo"; else rm -f "$$depbase.Tpo"; exit 1; fi
#	source='$<' object='$@' libtool=yes \
#	DEPDIR=$(DEPDIR) $(CCDEPMODE) $(depcomp) \
#	$(LTCOMPILE) -c -o $@ $<

t_nc-t_nc.o: t_nc.c
	if $(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(t_nc_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS) -MT t_nc-t_nc.o -MD -MP -MF "$(DEPDIR)/t_nc-t_nc.Tpo" -c -o t_nc-t_nc.o `test -f 't_nc.c' || echo '$(srcdir)/'`t_nc.c; \
	then mv -f "$(DEPDIR)/t_nc-t_nc.Tpo" "$(DEPDIR)/t_nc-t_nc.Po"; else rm -f "$(DEPDIR)/t_nc-t_nc.Tpo"; exit 1; fi
#	source='t_nc.c' object='t_nc-t_nc.o' libtool=no \
#	DEPDIR=$(DEPDIR) $(CCDEPMODE) $(depcomp) \
#	$(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(t_nc_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS) -c -o t_nc-t_nc.o `test -f 't_nc.c' || echo '$(srcdir)/'`t_nc.c

t_nc-t_nc.obj: t_nc.c
	if $(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(t_nc_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS) -MT t_nc-t_nc.obj -MD -MP -MF "$(DEPDIR)/t_nc-t_nc.Tpo" -c -o t_nc-t_nc.obj `if test -f 't_nc.c'; then $(CYGPATH_W) 't_nc.c'; else $(CYGPATH_W) '$(srcdir)/t_nc.c'; fi`; \
	then mv -f "$(DEPDIR)/t_nc-t_nc.Tpo" "$(DEPDIR)/t_nc-t_nc.Po"; else rm -f "$(DEPDIR)/t_nc-t_nc.Tpo"; exit 1; fi
#	source='t_nc.c' object='t_nc-t_nc.obj' libtool=no \
#	DEPDIR=$(DEPDIR) $(CCDEPMODE) $(depcomp) \
#	$(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(t_nc_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS) -c -o t_nc-t_nc.obj `if test -f 't_nc.c'; then $(CYGPATH_W) 't_nc.c'; else $(CYGPATH_W) '$(srcdir)/t_nc.c'; fi`

mostlyclean-libtool:
	-rm -f *.lo

clean-libtool:
	-rm -rf .libs _libs

distclean-libtool:
	-rm -f libtool
uninstall-info-am:
install-man3: $(man3_MANS) $(man_MANS)
	@$(NORMAL_INSTALL)
	test -z "$(man3dir)" || $(mkdir_p) "$(DESTDIR)$(man3dir)"
	@list='$(man3_MANS) $(dist_man3_MANS) $(nodist_man3_MANS)'; \
	l2='$(man_MANS) $(dist_man_MANS) $(nodist_man_MANS)'; \
	for i in $$l2; do \
	  case "$$i" in \
	    *.3*) list="$$list $$i" ;; \
	  esac; \
	done; \
	for i in $$list; do \
	  if test -f $(srcdir)/$$i; then file=$(srcdir)/$$i; \
	  else file=$$i; fi; \
	  ext=`echo $$i | sed -e 's/^.*\\.//'`; \
	  case "$$ext" in \
	    3*) ;; \
	    *) ext='3' ;; \
	  esac; \
	  inst=`echo $$i | sed -e 's/\\.[0-9a-z]*$$//'`; \
	  inst=`echo $$inst | sed -e 's/^.*\///'`; \
	  inst=`echo $$inst | sed '$(transform)'`.$$ext; \
	  echo " $(INSTALL_DATA) '$$file' '$(DESTDIR)$(man3dir)/$$inst'"; \
	  $(INSTALL_DATA) "$$file" "$(DESTDIR)$(man3dir)/$$inst"; \
	done
uninstall-man3:
	@$(NORMAL_UNINSTALL)
	@list='$(man3_MANS) $(dist_man3_MANS) $(nodist_man3_MANS)'; \
	l2='$(man_MANS) $(dist_man_MANS) $(nodist_man_MANS)'; \
	for i in $$l2; do \
	  case "$$i" in \
	    *.3*) list="$$list $$i" ;; \
	  esac; \
	done; \
	for i in $$list; do \
	  ext=`echo $$i | sed -e 's/^.*\\.//'`; \
	  case "$$ext" in \
	    3*) ;; \
	    *) ext='3' ;; \
	  esac; \
	  inst=`echo $$i | sed -e 's/\\.[0-9a-z]*$$//'`; \
	  inst=`echo $$inst | sed -e 's/^.*\///'`; \
	  inst=`echo $$inst | sed '$(transform)'`.$$ext; \
	  echo " rm -f '$(DESTDIR)$(man3dir)/$$inst'"; \
	  rm -f "$(DESTDIR)$(man3dir)/$$inst"; \
	done
install-includeHEADERS: $(include_HEADERS)
	@$(NORMAL_INSTALL)
	test -z "$(includedir)" || $(mkdir_p) "$(DESTDIR)$(includedir)"
	@list='$(include_HEADERS)'; for p in $$list; do \
	  if test -f "$$p"; then d=; else d="$(srcdir)/"; fi; \
	  f=$(am__strip_dir) \
	  echo " $(includeHEADERS_INSTALL) '$$d$$p' '$(DESTDIR)$(includedir)/$$f'"; \
	  $(includeHEADERS_INSTALL) "$$d$$p" "$(DESTDIR)$(includedir)/$$f"; \
	done

uninstall-includeHEADERS:
	@$(NORMAL_UNINSTALL)
	@list='$(include_HEADERS)'; for p in $$list; do \
	  f=$(am__strip_dir) \
	  echo " rm -f '$(DESTDIR)$(includedir)/$$f'"; \
	  rm -f "$(DESTDIR)$(includedir)/$$f"; \
	done

ID: $(HEADERS) $(SOURCES) $(LISP) $(TAGS_FILES)
	list='$(SOURCES) $(HEADERS) $(LISP) $(TAGS_FILES)'; \
	unique=`for i in $$list; do \
	    if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | \
	  $(AWK) '    { files[$$0] = 1; } \
	       END { for (i in files) print i; }'`; \
	mkid -fID $$unique
tags: TAGS

TAGS:  $(HEADERS) $(SOURCES)  $(TAGS_DEPENDENCIES) \
		$(TAGS_FILES) $(LISP)
	tags=; \
	here=`pwd`; \
	list='$(SOURCES) $(HEADERS)  $(LISP) $(TAGS_FILES)'; \
	unique=`for i in $$list; do \
	    if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | \
	  $(AWK) '    { files[$$0] = 1; } \
	       END { for (i in files) print i; }'`; \
	if test -z "$(ETAGS_ARGS)$$tags$$unique"; then :; else \
	  test -n "$$unique" || unique=$$empty_fix; \
	  $(ETAGS) $(ETAGSFLAGS) $(AM_ETAGSFLAGS) $(ETAGS_ARGS) \
	    $$tags $$unique; \
	fi
ctags: CTAGS
CTAGS:  $(HEADERS) $(SOURCES)  $(TAGS_DEPENDENCIES) \
		$(TAGS_FILES) $(LISP)
	tags=; \
	here=`pwd`; \
	list='$(SOURCES) $(HEADERS)  $(LISP) $(TAGS_FILES)'; \
	unique=`for i in $$list; do \
	    if test -f "$$i"; then echo $$i; else echo $(srcdir)/$$i; fi; \
	  done | \
	  $(AWK) '    { files[$$0] = 1; } \
	       END { for (i in files) print i; }'`; \
	test -z "$(CTAGS_ARGS)$$tags$$unique" \
	  || $(CTAGS) $(CTAGSFLAGS) $(AM_CTAGSFLAGS) $(CTAGS_ARGS) \
	     $$tags $$unique

GTAGS:
	here=`$(am__cd) $(top_builddir) && pwd` \
	  && cd $(top_srcdir) \
	  && gtags -i $(GTAGS_ARGS) $$here

distclean-tags:
	-rm -f TAGS ID GTAGS GRTAGS GSYMS GPATH tags

check-TESTS: $(TESTS)
	@failed=0; all=0; xfail=0; xpass=0; skip=0; \
	srcdir=$(srcdir); export srcdir; \
	list='$(TESTS)'; \
	if test -n "$$list"; then \
	  for tst in $$list; do \
	    if test -f ./$$tst; then dir=./; \
	    elif test -f $$tst; then dir=; \
	    else dir="$(srcdir)/"; fi; \
	    if $(TESTS_ENVIRONMENT) $${dir}$$tst; then \
	      all=`expr $$all + 1`; \
	      case " $(XFAIL_TESTS) " in \
	      *" $$tst "*) \
		xpass=`expr $$xpass + 1`; \
		failed=`expr $$failed + 1`; \
		echo "XPASS: $$tst"; \
	      ;; \
	      *) \
		echo "PASS: $$tst"; \
	      ;; \
	      esac; \
	    elif test $$? -ne 77; then \
	      all=`expr $$all + 1`; \
	      case " $(XFAIL_TESTS) " in \
	      *" $$tst "*) \
		xfail=`expr $$xfail + 1`; \
		echo "XFAIL: $$tst"; \
	      ;; \
	      *) \
		failed=`expr $$failed + 1`; \
		echo "FAIL: $$tst"; \
	      ;; \
	      esac; \
	    else \
	      skip=`expr $$skip + 1`; \
	      echo "SKIP: $$tst"; \
	    fi; \
	  done; \
	  if test "$$failed" -eq 0; then \
	    if test "$$xfail" -eq 0; then \
	      banner="All $$all tests passed"; \
	    else \
	      banner="All $$all tests behaved as expected ($$xfail expected failures)"; \
	    fi; \
	  else \
	    if test "$$xpass" -eq 0; then \
	      banner="$$failed of $$all tests failed"; \
	    else \
	      banner="$$failed of $$all tests did not behave as expected ($$xpass unexpected passes)"; \
	    fi; \
	  fi; \
	  dashes="$$banner"; \
	  skipped=""; \
	  if test "$$skip" -ne 0; then \
	    skipped="($$skip tests were not run)"; \
	    test `echo "$$skipped" | wc -c` -le `echo "$$banner" | wc -c` || \
	      dashes="$$skipped"; \
	  fi; \
	  report=""; \
	  if test "$$failed" -ne 0 && test -n "$(PACKAGE_BUGREPORT)"; then \
	    report="Please report to $(PACKAGE_BUGREPORT)"; \
	    test `echo "$$report" | wc -c` -le `echo "$$banner" | wc -c` || \
	      dashes="$$report"; \
	  fi; \
	  dashes=`echo "$$dashes" | sed s/./=/g`; \
	  echo "$$dashes"; \
	  echo "$$banner"; \
	  test -z "$$skipped" || echo "$$skipped"; \
	  test -z "$$report" || echo "$$report"; \
	  echo "$$dashes"; \
	  test "$$failed" -eq 0; \
	else :; fi

distdir: $(DISTFILES)
	@srcdirstrip=`echo "$(srcdir)" | sed 's|.|.|g'`; \
	topsrcdirstrip=`echo "$(top_srcdir)" | sed 's|.|.|g'`; \
	list='$(DISTFILES)'; for file in $$list; do \
	  case $$file in \
	    $(srcdir)/*) file=`echo "$$file" | sed "s|^$$srcdirstrip/||"`;; \
	    $(top_srcdir)/*) file=`echo "$$file" | sed "s|^$$topsrcdirstrip/|$(top_builddir)/|"`;; \
	  esac; \
	  if test -f $$file || test -d $$file; then d=.; else d=$(srcdir); fi; \
	  dir=`echo "$$file" | sed -e 's,/[^/]*$$,,'`; \
	  if test "$$dir" != "$$file" && test "$$dir" != "."; then \
	    dir="/$$dir"; \
	    $(mkdir_p) "$(distdir)$$dir"; \
	  else \
	    dir=''; \
	  fi; \
	  if test -d $$d/$$file; then \
	    if test -d $(srcdir)/$$file && test $$d != $(srcdir); then \
	      cp -pR $(srcdir)/$$file $(distdir)$$dir || exit 1; \
	    fi; \
	    cp -pR $$d/$$file $(distdir)$$dir || exit 1; \
	  else \
	    test -f $(distdir)/$$file \
	    || cp -p $$d/$$file $(distdir)/$$file \
	    || exit 1; \
	  fi; \
	done
check-am: all-am
	$(MAKE) $(AM_MAKEFLAGS) $(check_PROGRAMS)
	$(MAKE) $(AM_MAKEFLAGS) check-TESTS
check: check-am
all-am: Makefile $(LTLIBRARIES) $(MANS) $(HEADERS)
installdirs:
	for dir in "$(DESTDIR)$(libdir)" "$(DESTDIR)$(man3dir)" "$(DESTDIR)$(includedir)"; do \
	  test -z "$$dir" || $(mkdir_p) "$$dir"; \
	done
install: install-am
install-exec: install-exec-am
install-data: install-data-am
uninstall: uninstall-am

install-am: all-am
	@$(MAKE) $(AM_MAKEFLAGS) install-exec-am install-data-am

installcheck: installcheck-am
install-strip:
	$(MAKE) $(AM_MAKEFLAGS) INSTALL_PROGRAM="$(INSTALL_STRIP_PROGRAM)" \
	  install_sh_PROGRAM="$(INSTALL_STRIP_PROGRAM)" INSTALL_STRIP_FLAG=-s \
	  `test -z '$(STRIP)' || \
	    echo "INSTALL_PROGRAM_ENV=STRIPPROG='$(STRIP)'"` install
mostlyclean-generic:

clean-generic:
	-test -z "$(CLEANFILES)" || rm -f $(CLEANFILES)

distclean-generic:
	-test -z "$(CONFIG_CLEAN_FILES)" || rm -f $(CONFIG_CLEAN_FILES)

maintainer-clean-generic:
	@echo "This command is intended for maintainers to use"
	@echo "it deletes files that may require special tools to rebuild."
	-test -z "$(MAINTAINERCLEANFILES)" || rm -f $(MAINTAINERCLEANFILES)
clean: clean-am

clean-am: clean-checkPROGRAMS clean-generic clean-libLTLIBRARIES \
	clean-libtool clean-noinstLTLIBRARIES mostlyclean-am

distclean: distclean-am
	-rm -rf ./$(DEPDIR)
	-rm -f Makefile
distclean-am: clean-am distclean-compile distclean-generic \
	distclean-libtool distclean-tags

dvi: dvi-am

dvi-am:

html: html-am

info: info-am

info-am:

install-data-am: install-includeHEADERS install-man

install-exec-am: install-libLTLIBRARIES

install-info: install-info-am

install-man: install-man3

installcheck-am:

maintainer-clean: maintainer-clean-am
	-rm -rf ./$(DEPDIR)
	-rm -f Makefile
maintainer-clean-am: distclean-am maintainer-clean-generic

mostlyclean: mostlyclean-am

mostlyclean-am: mostlyclean-compile mostlyclean-generic \
	mostlyclean-libtool

pdf: pdf-am

pdf-am:

ps: ps-am

ps-am:

uninstall-am: uninstall-includeHEADERS uninstall-info-am \
	uninstall-libLTLIBRARIES uninstall-man

uninstall-man: uninstall-man3

.PHONY: CTAGS GTAGS all all-am check check-TESTS check-am clean \
	clean-checkPROGRAMS clean-generic clean-libLTLIBRARIES \
	clean-libtool clean-noinstLTLIBRARIES ctags distclean \
	distclean-compile distclean-generic distclean-libtool \
	distclean-tags distdir dvi dvi-am html html-am info info-am \
	install install-am install-data install-data-am install-exec \
	install-exec-am install-includeHEADERS install-info \
	install-info-am install-libLTLIBRARIES install-man \
	install-man3 install-strip installcheck installcheck-am \
	installdirs maintainer-clean maintainer-clean-generic \
	mostlyclean mostlyclean-compile mostlyclean-generic \
	mostlyclean-libtool pdf pdf-am ps ps-am tags uninstall \
	uninstall-am uninstall-includeHEADERS uninstall-info-am \
	uninstall-libLTLIBRARIES uninstall-man uninstall-man3


# This rule generates the C manpage.
netcdf.3: $(top_srcdir)/man/netcdf.m4 
	m4 $(M4FLAGS) -DAPI=C $? >$@  || rm $@

# This tells make how to turn .m4 files into .c files.
.m4.c:
	m4 $(AM_M4FLAGS) $(M4FLAGS) $< >$@

test: check
# Tell versions [3.59,3.63) of GNU make to not export all variables.
# Otherwise a system limit (for SysV at least) may be exceeded.
.NOEXPORT:
