## Process this file with automake to produce Makefile.in

AM_CFLAGS            = # Append below
AM_CPPFLAGS         = -I$(top_builddir)/gnulib -I$(top_srcdir)/gnulib
AM_CPPFLAGS        += -I$(top_builddir)        -I$(top_srcdir)
AM_FCFLAGS          = -I$(top_builddir)        -I$(top_srcdir)
AM_FCFLAGS         += -I$(top_builddir)/esio
AM_LDFLAGS          = # Append below
check_PROGRAMS      = # Append below
dist_check_SCRIPTS  = # Append below
noinst_PROGRAMS     = # Append below
BUILT_SOURCES       = # Append below
CLEANFILES          = # Append below
noinst_HEADERS      = fct.h                  \
                      testframework_assert.h \
                      testutils.h
EXTRA_DIST          = layout_template.c      \
                      plane_template.c       \
                      line_template.c        \
                      attribute_template.c   \
                      empty.h5

TESTS              = # Append below
XFAIL_TESTS        = # Append below
TESTS_ENVIRONMENT  = ESIO_TEST_INPUT_DIR=$(srcdir)
TESTS_ENVIRONMENT += ESIO_TEST_OUTPUT_DIR=$(builddir)

# Use MPI compiler wrappers since HDF5 wrappers break autodependencies.
CC=@MPICC@
FC=@MPIFC@

## Use HDF5 includes as many tests use HDF5 API calls to verify functionality
AM_CFLAGS   += $(HDF5_CFLAGS)
AM_CPPFLAGS += $(HDF5_CPPFLAGS)
AM_FCFLAGS  += $(HDF5_FFLAGS)
AM_LDFLAGS  += $(HDF5_LDFLAGS)

#############################################################################
## INTERNAL INTERNAL INTERNAL INTERNAL INTERNAL INTERNAL INTERNAL INTERNAL ##
#############################################################################

## Internal restart renaming logic tests
TESTS                   += restart_helpers
check_PROGRAMS          += restart_helpers
restart_helpers_SOURCES  = restart_helpers.c ../esio/restart-rename.c
restart_helpers_LDADD    = ../esio/libesio.la

## Internal restart renaming logic tests
TESTS             += uri_tests
check_PROGRAMS    += uri_tests
uri_tests_SOURCES  = uri_tests.c ../esio/uri.c
uri_tests_LDADD    = ../esio/libesio.la

#############################################################################
## BASIC BASIC BASIC BASIC BASIC BASIC BASIC BASIC BASIC BASIC BASIC BASIC ##
#############################################################################

## Basic C API sanity checks
TESTS              += basic.sh
dist_check_SCRIPTS += basic.sh
check_PROGRAMS     += basic
basic_SOURCES       = basic.c testutils.c
basic_LDADD         = ../esio/libesio.la

## Ensure Fortran testframework can return non-zero exit status
TESTS            += sanity_f
XFAIL_TESTS      += sanity_f
check_PROGRAMS   += sanity_f
sanity_f_SOURCES  = sanity_f.F90 testutils.c testframework.F90
sanity_f_LDADD    = ../esio/libesiof.la

## Basic Fortran API sanity checks
TESTS           += basic_f
check_PROGRAMS  += basic_f
basic_f_SOURCES  = basic_f.F90 testutils.c testframework.F90
basic_f_LDADD    = ../esio/libesiof.la

###########################################################################
## STRING STRING STRING STRING STRING STRING STRING STRING STRING STRING ##
###########################################################################

## String tests...
TESTS              += string.sh
dist_check_SCRIPTS += string.sh

## ...in C
check_PROGRAMS += string
string_SOURCES  = string.c testutils.c
string_LDADD    = ../esio/libesio.la

## ...in Fortran
check_PROGRAMS   += string_f
string_f_SOURCES  = string_f.F90 testutils.c testframework.F90
string_f_LDADD    = ../esio/libesiof.la

###########################################################################
## ATTRIBUTE ATTRIBUTE ATTRIBUTE ATTRIBUTE ATTRIBUTE ATTRIBUTE ATTRIBUTE ##
###########################################################################

## Attribute double precision tests in C
TESTS                    += attribute_double.sh
dist_check_SCRIPTS       += attribute_double.sh
check_PROGRAMS           += attribute_double
attribute_double_SOURCES  = attribute_double.c testutils.c
attribute_double_LDADD    = ../esio/libesio.la

## Attribute single precision tests in C
TESTS                   += attribute_float.sh
dist_check_SCRIPTS      += attribute_float.sh
check_PROGRAMS          += attribute_float
attribute_float_SOURCES  = attribute_float.c testutils.c
attribute_float_LDADD    = ../esio/libesio.la

## Attribute integer tests...
TESTS                 += attribute_int.sh
dist_check_SCRIPTS    += attribute_int.sh

## ...in C
check_PROGRAMS        += attribute_int
attribute_int_SOURCES  = attribute_int.c testutils.c
attribute_int_LDADD    = ../esio/libesio.la

## ...in Fortran
check_PROGRAMS          += attribute_int_f
attribute_int_f_SOURCES  = attribute_int_f.F90 testutils.c testframework.F90
attribute_int_f_LDADD    = ../esio/libesiof.la

###########################################################################
## LINE LINE LINE LINE LINE LINE LINE LINE LINE LINE LINE LINE LINE LINE ##
###########################################################################

## Line double precision tests in C
TESTS               += line_double.sh
dist_check_SCRIPTS  += line_double.sh
check_PROGRAMS      += line_double
line_double_SOURCES  = line_double.c testutils.c
line_double_LDADD    = ../esio/libesio.la

## Line single precision tests in C
TESTS               += line_float.sh
dist_check_SCRIPTS  += line_float.sh
check_PROGRAMS      += line_float
line_float_SOURCES   = line_float.c testutils.c
line_float_LDADD     = ../esio/libesio.la

## Line integer tests...
TESTS              += line_int.sh
dist_check_SCRIPTS += line_int.sh

## ...in C
check_PROGRAMS    += line_int
line_int_SOURCES   = line_int.c testutils.c
line_int_LDADD     = ../esio/libesio.la

## ...in Fortran
check_PROGRAMS     += line_int_f
line_int_f_SOURCES  = line_int_f.F90 testutils.c testframework.F90
line_int_f_LDADD    = ../esio/libesiof.la

#############################################################################
## PLANE PLANE PLANE PLANE PLANE PLANE PLANE PLANE PLANE PLANE PLANE PLANE ##
#############################################################################

## Plane double precision tests in C
TESTS                += plane_double.sh
dist_check_SCRIPTS   += plane_double.sh
check_PROGRAMS       += plane_double
plane_double_SOURCES  = plane_double.c testutils.c
plane_double_LDADD    = ../esio/libesio.la

## Plane single precision tests in C
TESTS                += plane_float.sh
dist_check_SCRIPTS   += plane_float.sh
check_PROGRAMS       += plane_float
plane_float_SOURCES   = plane_float.c testutils.c
plane_float_LDADD     = ../esio/libesio.la

## Plane integer tests...
TESTS              += plane_int.sh
dist_check_SCRIPTS += plane_int.sh

## ...in C
check_PROGRAMS     += plane_int
plane_int_SOURCES   = plane_int.c testutils.c
plane_int_LDADD     = ../esio/libesio.la

## ...in Fortran
check_PROGRAMS      += plane_int_f
plane_int_f_SOURCES  = plane_int_f.F90 testutils.c testframework.F90
plane_int_f_LDADD    = ../esio/libesiof.la

#############################################################################
## LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0 LAYOUT 0 ##
#############################################################################

## Layout 0 double precision tests in C
TESTS                  += layout0_double.sh
dist_check_SCRIPTS     += layout0_double.sh
check_PROGRAMS         += layout0_double
layout0_double_SOURCES  = layout0_double.c testutils.c
layout0_double_LDADD    = ../esio/libesio.la

## Layout 0 single precision tests in C
TESTS                  += layout0_float.sh
dist_check_SCRIPTS     += layout0_float.sh
check_PROGRAMS         += layout0_float
layout0_float_SOURCES   = layout0_float.c testutils.c
layout0_float_LDADD     = ../esio/libesio.la

## Layout 0 integer tests...
TESTS                += layout0_int.sh
dist_check_SCRIPTS   += layout0_int.sh

## ...in C
check_PROGRAMS       += layout0_int
layout0_int_SOURCES   = layout0_int.c testutils.c
layout0_int_LDADD     = ../esio/libesio.la

## ...in Fortran
check_PROGRAMS        += layout0_int_f
layout0_int_f_SOURCES  = layout0_int_f.F90 testutils.c testframework.F90
layout0_int_f_LDADD    = ../esio/libesiof.la

#############################################################################
## LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 LAYOUT 1 ##
#############################################################################

## Layout 1 double precision tests in C
TESTS                  += layout1_double.sh
dist_check_SCRIPTS     += layout1_double.sh
check_PROGRAMS         += layout1_double
layout1_double_SOURCES  = layout1_double.c testutils.c
layout1_double_LDADD    = ../esio/libesio.la

## Layout 1 single precision tests in C
TESTS                  += layout1_float.sh
dist_check_SCRIPTS     += layout1_float.sh
check_PROGRAMS         += layout1_float
layout1_float_SOURCES   = layout1_float.c testutils.c
layout1_float_LDADD     = ../esio/libesio.la

## Layout 1 integer tests in C
TESTS                += layout1_int.sh
dist_check_SCRIPTS   += layout1_int.sh
check_PROGRAMS       += layout1_int
layout1_int_SOURCES   = layout1_int.c testutils.c
layout1_int_LDADD     = ../esio/libesio.la

#############################################################################
## LAYOUT 2 LAYOUT 2 LAYOUT 2 LAYOUT 2 LAYOUT 2 LAYOUT 2 LAYOUT 2 LAYOUT 2 ##
#############################################################################

## Layout 2 double precision tests in C
TESTS                  += layout2_double.sh
dist_check_SCRIPTS     += layout2_double.sh
check_PROGRAMS         += layout2_double
layout2_double_SOURCES  = layout2_double.c testutils.c
layout2_double_LDADD    = ../esio/libesio.la

## Layout 2 single precision tests in C
TESTS                  += layout2_float.sh
dist_check_SCRIPTS     += layout2_float.sh
check_PROGRAMS         += layout2_float
layout2_float_SOURCES   = layout2_float.c testutils.c
layout2_float_LDADD     = ../esio/libesio.la

## Layout 2 integer tests in C
TESTS                += layout2_int.sh
dist_check_SCRIPTS   += layout2_int.sh
check_PROGRAMS       += layout2_int
layout2_int_SOURCES   = layout2_int.c testutils.c
layout2_int_LDADD     = ../esio/libesio.la

#############################################################################
### INSTALLED APPLICATIONS INSTALLED APPLICATIONS INSTALLED APPLICATIONS  ###
#############################################################################

## Externally observable restart renaming tets
TESTS              += esio_rename.sh
dist_check_SCRIPTS += esio_rename.sh

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
SFMAKEDEPEND_FLAGS = -i esio -i esio_c_binding
include $(top_srcdir)/sfmakedepend.am

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
