## Process this file with automake to produce Makefile.in

AM_CFLAGS      = -I$(top_srcdir) -I$(top_builddir)
AM_CFLAGS     += -I$(top_srcdir)/lib -I$(top_builddir)/lib
AM_FCFLAGS     = -I$(top_srcdir) -I$(top_builddir)
AM_FCFLAGS    += -I$(top_srcdir)/lib -I$(top_builddir)/lib
AM_FCFLAGS    += -I$(top_srcdir)/esio -I$(top_builddir)/esio
AM_LDFLAGS     =
bin_PROGRAMS   =
BUILT_SOURCES  =
BUILT_DOCS     =
EXTRA_DIST     =

# Use parallel HDF5 compiler wrappers per configure.ac
CC=@H5CC@
FC=@H5FC@

## Pick up GRVY library installation
AM_CFLAGS  += $(GRVY_CFLAGS)
AM_FCFLAGS += $(GRVY_FCFLAGS)
AM_LDFLAGS += $(GRVY_LIBS)

# Sample progams
bin_PROGRAMS       += hellompi_c
hellompi_c_SOURCES  = hellompi_c.c

bin_PROGRAMS              += hellompi_f90
hellompi_f90_SOURCES       = hellompi_f90.f90
hellompi_f90_LIBTOOLFLAGS  = --tag=FC

# PSDNS iokernel test program
bin_PROGRAMS      += main
main_SOURCES       = header.f90 main.f90
main_LDADD         = ../esio/libesio.la
main_LIBTOOLFLAGS  = --tag=FC
# Fortran module dependencies
main.o: header.f90 main.f90

# Clean up Fortran modules on 'make clean'
CLEANFILES = header.mod

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
