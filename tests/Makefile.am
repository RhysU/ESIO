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
EXTRA_DIST      = layout_template.c plane_template.c

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

#############################################################################
## PLANE PLANE PLANE PLANE PLANE PLANE PLANE PLANE PLANE PLANE PLANE PLANE ##
#############################################################################

## Plane double precision tests in C
TESTS                += plane_double.sh    # NB Wrapper script
EXTRA_DIST           += plane_double.sh    # NB Wrapper script
check_PROGRAMS       += plane_double
plane_double_SOURCES  = plane_double.c
plane_double_LDADD    = ../esio/libesio.la
plane_double.o: plane_double.c plane_template.c

## Plane single precision tests in C
TESTS                += plane_float.sh          # NB Wrapper script
EXTRA_DIST           += plane_float.sh          # NB Wrapper script
check_PROGRAMS       += plane_float
plane_float_SOURCES   = plane_float.c
plane_float_LDADD     = ../esio/libesio.la
plane_float.o: plane_float.c plane_template.c

## Plane integer tests in C
TESTS              += plane_int.sh          # NB Wrapper script
EXTRA_DIST         += plane_int.sh          # NB Wrapper script
check_PROGRAMS     += plane_int
plane_int_SOURCES   = plane_int.c
plane_int_LDADD     = ../esio/libesio.la
plane_int.o: plane_int.c plane_template.c

#############################################################################
## LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0 ##
#############################################################################

## Layout 0 double precision tests in C
TESTS                  += layout0_double.sh    # NB Wrapper script
EXTRA_DIST             += layout0_double.sh    # NB Wrapper script
check_PROGRAMS         += layout0_double
layout0_double_SOURCES  = layout0_double.c
layout0_double_LDADD    = ../esio/libesio.la
layout0_double.o: layout0_double.c layout_template.c

## Layout 0 single precision tests in C
TESTS                  += layout0_float.sh          # NB Wrapper script
EXTRA_DIST             += layout0_float.sh          # NB Wrapper script
check_PROGRAMS         += layout0_float
layout0_float_SOURCES   = layout0_float.c
layout0_float_LDADD     = ../esio/libesio.la
layout0_float.o: layout0_float.c layout_template.c

## Layout 0 integer tests in C
TESTS                += layout0_int.sh          # NB Wrapper script
EXTRA_DIST           += layout0_int.sh          # NB Wrapper script
check_PROGRAMS       += layout0_int
layout0_int_SOURCES   = layout0_int.c
layout0_int_LDADD     = ../esio/libesio.la
layout0_int.o: layout0_int.c layout_template.c

#############################################################################
## LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 ##
#############################################################################

## Layout 1 double precision tests in C
TESTS                  += layout1_double.sh    # NB Wrapper script
EXTRA_DIST             += layout1_double.sh    # NB Wrapper script
check_PROGRAMS         += layout1_double
layout1_double_SOURCES  = layout1_double.c
layout1_double_LDADD    = ../esio/libesio.la
layout1_double.o: layout1_double.c layout_template.c

## Layout 1 single precision tests in C
TESTS                  += layout1_float.sh          # NB Wrapper script
EXTRA_DIST             += layout1_float.sh          # NB Wrapper script
check_PROGRAMS         += layout1_float
layout1_float_SOURCES   = layout1_float.c
layout1_float_LDADD     = ../esio/libesio.la
layout1_float.o: layout1_float.c layout_template.c

## Layout 1 integer tests in C
TESTS                += layout1_int.sh          # NB Wrapper script
EXTRA_DIST           += layout1_int.sh          # NB Wrapper script
check_PROGRAMS       += layout1_int
layout1_int_SOURCES   = layout1_int.c
layout1_int_LDADD     = ../esio/libesio.la
layout1_int.o: layout1_int.c layout_template.c

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
