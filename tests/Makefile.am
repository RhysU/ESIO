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
EXTRA_DIST      = layout_template.c

TESTS              = # Append below
TESTS_ENVIRONMENT  = ESIO_TEST_INPUT_DIR=$(srcdir)
TESTS_ENVIRONMENT += ESIO_TEST_OUPUT_DIR=.

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

## Basic C API sanity checks
TESTS          += basic
check_PROGRAMS += basic
basic_SOURCES   = basic.c
basic_LDADD     = ../esio/libesio.la

## Double precision fields from C
TESTS                  += layout0_double.sh    # NB Wrapper script
EXTRA_DIST             += layout0_double.sh    # NB Wrapper script
check_PROGRAMS         += layout0_double
layout0_double_SOURCES  = layout0_double.c
layout0_double_LDADD    = ../esio/libesio.la
layout0_double.o: layout0_double.c layout_template.c

## Single precision fields from C
TESTS                  += layout0_float.sh          # NB Wrapper script
EXTRA_DIST             += layout0_float.sh          # NB Wrapper script
check_PROGRAMS         += layout0_float
layout0_float_SOURCES   = layout0_float.c
layout0_float_LDADD     = ../esio/libesio.la
layout0_float.o: layout0_float.c layout_template.c

# FIXME Disabled Fortran test
## Fortran test program to read/write a field
#noinst_PROGRAMS        += simplef90
#simplef90_SOURCES       = header.f90 simplef90.f90
#simplef90_LDADD         = ../esio/libesio.la
#simplef90_LIBTOOLFLAGS  = --tag=FC

# Build test programs without running 'make check'
.PHONY: check-programs
check-programs: $(check_PROGRAMS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
