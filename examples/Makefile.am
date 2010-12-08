## Process this file with automake to produce Makefile.in

AM_CFLAGS       =
AM_CPPFLAGS     = -I$(top_srcdir)     -I$(top_builddir)
AM_CPPFLAGS    += -I$(top_srcdir)/lib -I$(top_builddir)/lib
AM_FCFLAGS      = -I$(top_srcdir)     -I$(top_builddir)
AM_FCFLAGS     += -I$(top_srcdir)/lib -I$(top_builddir)/lib
AM_FCFLAGS     += -I$(top_builddir)/esio
AM_LDFLAGS      = -L$(top_builddir)/lib
noinst_PROGRAMS = # Append below
BUILT_DOCS      = # Append below
BUILT_SOURCES   = # Append below
CLEANFILES      = # Append below
EXTRA_DIST      = # Append below

# Use MPI compiler wrappers since HDF5 wrappers break autodependencies.
CC=@MPICC@
FC=@MPIFC@

## Pick up HDF5 library includes
AM_CFLAGS   += $(HDF5_CFLAGS)
AM_CPPFLAGS += $(HDF5_CPPFLAGS)
AM_FCFLAGS  += $(HDF5_FFLAGS)

# Sample progams
noinst_PROGRAMS        += concepts1
concepts1_SOURCES       = concepts1.c
concepts1_LDADD         = ../esio/libesio.la $(HDF5_LDFLAGS)

noinst_PROGRAMS        += concepts2
concepts2_SOURCES       = concepts2.f90
concepts2_LDADD         = ../esio/libesiof.la $(HDF5_LDFLAGS)
concepts2_LIBTOOLFLAGS  = --tag=FC

# Additional files to clean
CLEANFILES += *.mod *.gcno *.gcda

# License updates
include $(top_srcdir)/license.am

# Automatic dependency generation for Fortran
include $(top_srcdir)/sfmakedepend.am

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
