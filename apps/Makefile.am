## Process this file with automake to produce Makefile.in

AM_CFLAGS       =
AM_CPPFLAGS     = -I$(top_srcdir)     -I$(top_builddir)
AM_CPPFLAGS    += -I$(top_srcdir)/lib -I$(top_builddir)/lib
AM_FCFLAGS      = -I$(top_srcdir)     -I$(top_builddir)
AM_FCFLAGS     += -I$(top_srcdir)/lib -I$(top_builddir)/lib
AM_FCFLAGS     += -I$(top_builddir)/esio
AM_LDFLAGS      = -L$(top_builddir)/lib
bin_PROGRAMS    = # Append below
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

## Pick up GRVY library includes
AM_CFLAGS  += $(GRVY_CFLAGS)
AM_FCFLAGS += $(GRVY_FCFLAGS)

# Sample progams
bin_PROGRAMS       += hellompi_c
hellompi_c_SOURCES  = hellompi_c.c

bin_PROGRAMS              += hellompi_f90
hellompi_f90_SOURCES       = hellompi_f90.f90
hellompi_f90_LIBTOOLFLAGS  = --tag=FC

# PSDNS iokernel test program
noinst_PROGRAMS   += main
main_SOURCES       = header.f90 main.f90
main_LDADD         = ../esio/libesiof.la $(GRVY_LIBS) $(HDF5_LDFLAGS)
main_LIBTOOLFLAGS  = --tag=FC
main.f90: header.o                       # Explicit module dependence

# Additional files to clean
CLEANFILES += *.mod *.gcno

# License updates
include $(top_srcdir)/license.am

# Automatic dependency generation for Fortran
include $(top_srcdir)/sfmakedepend.am

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
