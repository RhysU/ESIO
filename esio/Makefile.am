## Process this file with automake to produce Makefile.in

AM_CFLAGS     = -I$(top_srcdir) -I$(top_builddir)
AM_CFLAGS    += -I$(top_srcdir)/lib -I$(top_builddir)/lib
AM_FCFLAGS    = -I$(top_srcdir) -I$(top_builddir)
AM_FCFLAGS   += -I$(top_srcdir)/lib -I$(top_builddir)/lib
AM_LDFLAGS         =
bin_PROGRAMS       = # Append below
check_PROGRAMS     = # Append below
lib_LTLIBRARIES    = # Append below
noinst_HEADERS     = # Append below
noinst_LTLIBRARIES = # Append below
pkginclude_HEADERS = # Append below
CLEANFILES         = # Append below

# Use parallel HDF5 compiler wrappers per configure.ac
CC=@H5CC@
FC=@H5FC@

# Headers installed with the library
pkginclude_HEADERS += esio.h
pkginclude_HEADERS += error.h

# Public ESIO library containing C and Fortran interfaces
# Need distinct C and Fortran sublibraries for libtool --tag purposes
lib_LTLIBRARIES     += libesio.la
libesio_la_LDFLAGS   = -release $(PACKAGE_VERSION)
libesio_la_SOURCES   = # Required even if empty
libesio_la_LIBADD    = libesio_c.la libesio_f.la

# C sources to place into the public library
noinst_LTLIBRARIES   += libesio_c.la
libesio_c_la_SOURCES  = esiow.c esiofb.c error.c

# Fortran sources to place into the public library
noinst_LTLIBRARIES        += libesio_f.la
libesio_f_la_SOURCES       = esiof.f90
libesio_f_la_LIBTOOLFLAGS  = --tag=FC

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
