## Process this file with automake to produce Makefile.in

AM_CFLAGS          =
AM_CPPFLAGS        = -I$(top_srcdir)/lib -I$(top_builddir)/lib
AM_FCFLAGS         = -I$(top_srcdir)/lib -I$(top_builddir)/lib
AM_LDFLAGS         = -L$(top_builddir)/lib
bin_PROGRAMS       = # Append below
check_PROGRAMS     = # Append below
lib_LTLIBRARIES    = # Append below
noinst_HEADERS     = # Append below
noinst_LTLIBRARIES = # Append below
pkginclude_HEADERS = # Append below
BUILT_SOURCES      = # Append below
CLEANFILES         = # Append below

# Use MPI compiler wrappers since HDF5 wrappers break autodependencies.
CC=@MPICC@
FC=@MPIFC@

# Headers shipped with the library
pkginclude_HEADERS += error.h
pkginclude_HEADERS += esio.h
pkginclude_HEADERS += version.h
# Headers internal to the library
noinst_HEADERS += layout.h
noinst_HEADERS += metadata.h

# Public ESIO library containing C and Fortran interfaces
# Need distinct C and Fortran sublibraries for libtool --tag purposes
lib_LTLIBRARIES     += libesio.la
libesio_la_LDFLAGS   = -release $(PACKAGE_VERSION)
libesio_la_SOURCES   = # Required even if empty
libesio_la_LIBADD    = libesio_c.la libesio_f.la

# C sources to place into the public library
noinst_LTLIBRARIES    += libesio_c.la
libesio_c_la_CFLAGS    = $(AM_CFLAGS) $(HDF5_CFLAGS)
libesio_c_la_CPPFLAGS  = $(AM_CPPFLAGS) $(HDF5_CPPFLAGS)
libesio_c_la_SOURCES   = error.c
libesio_c_la_SOURCES  += esio.c
libesio_c_la_SOURCES  += layout.c
libesio_c_la_SOURCES  += metadata.c

# Fortran sources to place into the public library
noinst_LTLIBRARIES        += libesio_f.la
libesio_f_la_FCFLAGS       = $(AM_FCFLAGS) $(HDF5_FFLAGS)
libesio_f_la_SOURCES       = esio.f90
libesio_f_la_LIBTOOLFLAGS  = --tag=FC

# Extra things to clean up
clean-local:
	rm -f *.mod *.gcno

# License updates
include $(top_srcdir)/license.am

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
