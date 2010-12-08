## Process this file with automake to produce Makefile.in

AM_CFLAGS       =
AM_CPPFLAGS     = -I$(top_srcdir)     -I$(top_builddir)
AM_CPPFLAGS    += -I$(top_srcdir)/lib -I$(top_builddir)/lib
AM_FCFLAGS      = -I$(top_srcdir)     -I$(top_builddir)
AM_FCFLAGS     += -I$(top_srcdir)/lib -I$(top_builddir)/lib
AM_FCFLAGS     += -I$(top_builddir)/esio
AM_LDFLAGS      = -L$(top_builddir)/lib
check_PROGRAMS  = # Append below
noinst_HEADERS  = fct.h
noinst_PROGRAMS = # Append below
BUILT_SOURCES   = # Append below
CLEANFILES      = # Append below
EXTRA_DIST      = layout_template.c \
                  plane_template.c  \
                  line_template.c   \
                  attribute_template.c

TESTS              = # Append below
TESTS_ENVIRONMENT  = ESIO_TEST_INPUT_DIR=$(srcdir)
TESTS_ENVIRONMENT += ESIO_TEST_OUPUT_DIR=.

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

#############################################################################
## BASIC BASIC BASIC BASIC BASIC BASIC BASIC BASIC BASIC BASIC BASIC BASIC ##
#############################################################################

## Basic C API sanity checks
TESTS          += basic
check_PROGRAMS += basic
basic_SOURCES   = basic.c
basic_LDADD     = ../esio/libesio.la $(GRVY_LIBS) $(HDF5_LDFLAGS)

## Basic Fortran API sanity checks
TESTS           += basic_f
check_PROGRAMS  += basic_f
basic_f_SOURCES  = basic_f.f90 testframework.f90
basic_f_LDADD    = ../esio/libesiof.la $(GRVY_LIBS) $(HDF5_LDFLAGS)

###########################################################################
## STRING STRING STRING STRING STRING STRING STRING STRING STRING STRING ##
###########################################################################

## String tests...
TESTS          += string.sh    # NB Wrapper script
EXTRA_DIST     += string.sh    # NB Wrapper script

## ...in C
check_PROGRAMS += string
string_SOURCES  = string.c
string_LDADD    = ../esio/libesio.la $(GRVY_LIBS) $(HDF5_LDFLAGS)

## ...in Fortran
check_PROGRAMS   += string_f
string_f_SOURCES  = string_f.f90 testframework.f90
string_f_LDADD    = ../esio/libesiof.la $(GRVY_LIBS) $(HDF5_LDFLAGS)

###########################################################################
## ATTRIBUTE ATTRIBUTE ATTRIBUTE ATTRIBUTE ATTRIBUTE ATTRIBUTE ATTRIBUTE ##
###########################################################################

## Attribute double precision tests in C
TESTS                    += attribute_double.sh    # NB Wrapper script
EXTRA_DIST               += attribute_double.sh    # NB Wrapper script
check_PROGRAMS           += attribute_double
attribute_double_SOURCES  = attribute_double.c
attribute_double_LDADD    = ../esio/libesio.la $(GRVY_LIBS) $(HDF5_LDFLAGS)

## Attribute single precision tests in C
TESTS                   += attribute_float.sh     # NB Wrapper script
EXTRA_DIST              += attribute_float.sh     # NB Wrapper script
check_PROGRAMS          += attribute_float
attribute_float_SOURCES  = attribute_float.c
attribute_float_LDADD    = ../esio/libesio.la $(GRVY_LIBS) $(HDF5_LDFLAGS)

## Attribute integer tests...
TESTS                 += attribute_int.sh         # NB Wrapper script
EXTRA_DIST            += attribute_int.sh         # NB Wrapper script

## ...in C
check_PROGRAMS        += attribute_int
attribute_int_SOURCES  = attribute_int.c
attribute_int_LDADD    = ../esio/libesio.la $(GRVY_LIBS) $(HDF5_LDFLAGS)

## ...in Fortran
check_PROGRAMS          += attribute_int_f
attribute_int_f_SOURCES  = attribute_int_f.f90 testframework.f90
attribute_int_f_LDADD    = ../esio/libesiof.la $(GRVY_LIBS) $(HDF5_LDFLAGS)

###########################################################################
## LINE LINE LINE LINE LINE LINE LINE LINE LINE LINE LINE LINE LINE LINE ##
###########################################################################

## Line double precision tests in C
TESTS               += line_double.sh    # NB Wrapper script
EXTRA_DIST          += line_double.sh    # NB Wrapper script
check_PROGRAMS      += line_double
line_double_SOURCES  = line_double.c
line_double_LDADD    = ../esio/libesio.la $(GRVY_LIBS) $(HDF5_LDFLAGS)

## Line single precision tests in C
TESTS               += line_float.sh     # NB Wrapper script
EXTRA_DIST          += line_float.sh     # NB Wrapper script
check_PROGRAMS      += line_float
line_float_SOURCES   = line_float.c
line_float_LDADD     = ../esio/libesio.la $(GRVY_LIBS) $(HDF5_LDFLAGS)

## Line integer tests...
TESTS             += line_int.sh         # NB Wrapper script
EXTRA_DIST        += line_int.sh         # NB Wrapper script

## ...in C
check_PROGRAMS    += line_int
line_int_SOURCES   = line_int.c
line_int_LDADD     = ../esio/libesio.la $(GRVY_LIBS) $(HDF5_LDFLAGS)

## ...in Fortran
check_PROGRAMS     += line_int_f
line_int_f_SOURCES  = line_int_f.f90 testframework.f90
line_int_f_LDADD    = ../esio/libesiof.la $(GRVY_LIBS) $(HDF5_LDFLAGS)

#############################################################################
## PLANE PLANE PLANE PLANE PLANE PLANE PLANE PLANE PLANE PLANE PLANE PLANE ##
#############################################################################

## Plane double precision tests in C
TESTS                += plane_double.sh    # NB Wrapper script
EXTRA_DIST           += plane_double.sh    # NB Wrapper script
check_PROGRAMS       += plane_double
plane_double_SOURCES  = plane_double.c
plane_double_LDADD    = ../esio/libesio.la $(GRVY_LIBS) $(HDF5_LDFLAGS)

## Plane single precision tests in C
TESTS                += plane_float.sh     # NB Wrapper script
EXTRA_DIST           += plane_float.sh     # NB Wrapper script
check_PROGRAMS       += plane_float
plane_float_SOURCES   = plane_float.c
plane_float_LDADD     = ../esio/libesio.la $(GRVY_LIBS) $(HDF5_LDFLAGS)

## Plane integer tests...
TESTS              += plane_int.sh         # NB Wrapper script
EXTRA_DIST         += plane_int.sh         # NB Wrapper script

## ...in C
check_PROGRAMS     += plane_int
plane_int_SOURCES   = plane_int.c
plane_int_LDADD     = ../esio/libesio.la $(GRVY_LIBS) $(HDF5_LDFLAGS)

## ...in Fortran
check_PROGRAMS      += plane_int_f
plane_int_f_SOURCES  = plane_int_f.f90 testframework.f90
plane_int_f_LDADD    = ../esio/libesiof.la $(GRVY_LIBS) $(HDF5_LDFLAGS)

#############################################################################
## LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0 ##
#############################################################################

## Layout 0 double precision tests in C
TESTS                  += layout0_double.sh    # NB Wrapper script
EXTRA_DIST             += layout0_double.sh    # NB Wrapper script
check_PROGRAMS         += layout0_double
layout0_double_SOURCES  = layout0_double.c
layout0_double_LDADD    = ../esio/libesio.la $(GRVY_LIBS) $(HDF5_LDFLAGS)

## Layout 0 single precision tests in C
TESTS                  += layout0_float.sh     # NB Wrapper script
EXTRA_DIST             += layout0_float.sh     # NB Wrapper script
check_PROGRAMS         += layout0_float
layout0_float_SOURCES   = layout0_float.c
layout0_float_LDADD     = ../esio/libesio.la $(GRVY_LIBS) $(HDF5_LDFLAGS)

## Layout 0 integer tests...
TESTS                += layout0_int.sh         # NB Wrapper script
EXTRA_DIST           += layout0_int.sh         # NB Wrapper script

## ...in C
check_PROGRAMS       += layout0_int
layout0_int_SOURCES   = layout0_int.c
layout0_int_LDADD     = ../esio/libesio.la $(GRVY_LIBS) $(HDF5_LDFLAGS)

## ...in Fortran
check_PROGRAMS        += layout0_int_f
layout0_int_f_SOURCES  = layout0_int_f.f90 testframework.f90
layout0_int_f_LDADD    = ../esio/libesiof.la $(GRVY_LIBS) $(HDF5_LDFLAGS)

#############################################################################
## LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 ##
#############################################################################

## Layout 1 double precision tests in C
TESTS                  += layout1_double.sh    # NB Wrapper script
EXTRA_DIST             += layout1_double.sh    # NB Wrapper script
check_PROGRAMS         += layout1_double
layout1_double_SOURCES  = layout1_double.c
layout1_double_LDADD    = ../esio/libesio.la $(GRVY_LIBS) $(HDF5_LDFLAGS)

## Layout 1 single precision tests in C
TESTS                  += layout1_float.sh     # NB Wrapper script
EXTRA_DIST             += layout1_float.sh     # NB Wrapper script
check_PROGRAMS         += layout1_float
layout1_float_SOURCES   = layout1_float.c
layout1_float_LDADD     = ../esio/libesio.la $(GRVY_LIBS) $(HDF5_LDFLAGS)

## Layout 1 integer tests in C
TESTS                += layout1_int.sh         # NB Wrapper script
EXTRA_DIST           += layout1_int.sh         # NB Wrapper script
check_PROGRAMS       += layout1_int
layout1_int_SOURCES   = layout1_int.c
layout1_int_LDADD     = ../esio/libesio.la $(GRVY_LIBS) $(HDF5_LDFLAGS)

#############################################################################
### MISCELLANEOUS MISCELLANEOUS MISCELLANEOUS MISCELLANEOUS MISCELLANEOUS ###
#############################################################################

# Build test programs without running 'make check'
.PHONY: check-programs
check-programs: $(check_PROGRAMS)

# Extra things to clean up
clean-local:
	rm -f *.mod *.gcno *.gcda

# License updates
include $(top_srcdir)/license.am

# Automatic dependency generation for Fortran
include $(top_srcdir)/sfmakedepend.am

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
