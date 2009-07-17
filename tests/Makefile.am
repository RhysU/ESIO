## Process this file with automake to produce Makefile.in

AM_CFLAGS     = -I$(top_srcdir) -I$(top_builddir)
AM_CFLAGS    += -I$(top_srcdir)/lib -I$(top_builddir)/lib
AM_FCFLAGS    = -I$(top_srcdir) -I$(top_builddir)
AM_FCFLAGS   += -I$(top_srcdir)/lib -I$(top_builddir)/lib
bin_PROGRAMS   = # Append below
check_PROGRAMS = # Append below
noinst_HEADERS = # Append below

# Use MPI compilers
CC=@MPICC@
FC=@MPIFC@

## Sample test program
check_PROGRAMS += sample
sample_SOURCES  = sample.c
sample_LDADD    = ../esio/libesio.la

TESTS_ENVIRONMENT =
TESTS = $(check_PROGRAMS)

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
