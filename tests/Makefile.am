## Process this file with automake to produce Makefile.in

AM_CFLAGS       = -I$(top_srcdir) -I$(top_builddir)
AM_CFLAGS      += -I$(top_srcdir)/lib -I$(top_builddir)/lib
AM_FCFLAGS      = -I$(top_srcdir) -I$(top_builddir)
AM_FCFLAGS     += -I$(top_srcdir)/lib -I$(top_builddir)/lib
AM_FCFLAGS     += -I$(top_srcdir)/esio -I$(top_builddir)/esio
AM_LDFLAGS      = # Append below
bin_PROGRAMS    = # Append below
check_PROGRAMS  = # Append below
noinst_HEADERS  = fct.h
noinst_PROGRAMS = # Append below

# Use parallel HDF5 compiler wrappers per configure.ac
CC=@H5CC@
FC=@H5FC@

## Pick up GRVY library installation
AM_CFLAGS  += $(GRVY_CFLAGS)
AM_FCFLAGS += $(GRVY_FCFLAGS)
AM_LDFLAGS += $(GRVY_LIBS)

# FIXME This should not be necessary for C-only tests
## Pick up Fortran support libraries
AM_LDFLAGS += $(FCLIBS)

## Sample test program
check_PROGRAMS += sample
sample_SOURCES  = sample.c
sample_LDADD    = ../esio/libesio.la

## Layout test program
check_PROGRAMS += basic
basic_SOURCES   = basic.c
basic_LDADD     = ../esio/libesio.la

# FIXME Disabled Fortran test
## Fortran test program to read/write a field
#noinst_PROGRAMS        += simplef90
#simplef90_SOURCES       = header.f90 simplef90.f90
#simplef90_LDADD         = ../esio/libesio.la
#simplef90_LIBTOOLFLAGS  = --tag=FC

TESTS_ENVIRONMENT  = ESIO_TEST_INPUT_DIR=$(srcdir)
TESTS_ENVIRONMENT += ESIO_TEST_OUPUT_DIR=.
TESTS = $(check_PROGRAMS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
