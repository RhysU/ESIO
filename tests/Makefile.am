## Process this file with automake to produce Makefile.in

AM_CFLAGS      =
AM_FCFLAGS     =
AM_LDFLAGS     =
bin_PROGRAMS   = # Append below
check_PROGRAMS = # Append below
noinst_HEADERS = # Append below

# Use MPI compilers
CC=@MPICC@
FC=@MPIFC@

## Sample test program
#check_PROGRAMS         += test_something
#test_something_SOURCES  = test_something.c
#test_something_LDADD    = ../esio/libesio.la

TESTS_ENVIRONMENT =
TESTS = $(check_PROGRAMS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
