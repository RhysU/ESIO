## Process this file with automake to produce Makefile.in

# Note this Makefile.am lacks Gnulib-related -I options like
#    AM_CPPFLAGS     = -I$(top_builddir)/gnulib -I$(top_srcdir)/gnulib
#    AM_FCFLAGS      = -I$(top_builddir)/gnulib -I$(top_srcdir)/gnulib
# This example code deliberately lacks #include <config.h> dependence
# and therefore cannot rely on Gnulib per
# http://www.gnu.org/software/gnulib/manual/gnulib.html#Source-changes

AM_CFLAGS       =
AM_CPPFLAGS     = -I$(top_builddir) -I$(top_srcdir)
AM_FCFLAGS      = -I$(top_builddir) -I$(top_srcdir)
AM_FCFLAGS     += -I$(top_builddir)/esio
AM_LDFLAGS      = # Append below
noinst_PROGRAMS = # Append below
BUILT_DOCS      = # Append below
BUILT_SOURCES   = # Append below
CLEANFILES      = # Append below

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
concepts2_SOURCES       = concepts2.F90
concepts2_LDADD         = ../esio/libesiof.la $(HDF5_LDFLAGS)
concepts2_LIBTOOLFLAGS  = --tag=FC

noinst_PROGRAMS       += restart1
restart1_SOURCES       = restart1.c
restart1_LDADD         = ../esio/libesio.la $(HDF5_LDFLAGS)

# Additional files to clean
CLEANFILES += *.mod *.gcno *.gcda

# License updates
include $(top_srcdir)/license.am

# Automatic dependency generation for Fortran
SFMAKEDEPEND_FLAGS = -i esio
include $(top_srcdir)/sfmakedepend.am

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
