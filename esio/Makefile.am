## Process this file with automake to produce Makefile.in

AM_CFLAGS          =
AM_FCFLAGS         =
AM_LDFLAGS         =
bin_PROGRAMS       = # Append below
check_PROGRAMS     = # Append below
noinst_HEADERS     = # Append below
noinst_LTLIBRARIES = # Append below
pkginclude_HEADERS = # Append below

# Use MPI compilers
CC=@MPICC@
FC=@MPIFC@

# Headers installed with the library
pkginclude_HEADERS += esio.h

lib_LTLIBRARIES     += libesio.la
libesio_la_SOURCES   = esio.c

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
