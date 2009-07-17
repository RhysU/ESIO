## Process this file with automake to produce Makefile.in

AM_CFLAGS     = -I$(top_srcdir) -I$(top_builddir)
AM_CFLAGS    += -I$(top_srcdir)/lib -I$(top_builddir)/lib
AM_FCFLAGS    = -I$(top_srcdir) -I$(top_builddir)
AM_FCFLAGS   += -I$(top_srcdir)/lib -I$(top_builddir)/lib
AM_LDFLAGS         =
bin_PROGRAMS       = # Append below
check_PROGRAMS     = # Append below
noinst_HEADERS     = # Append below
lib_LTLIBRARIES    = # Append below
pkginclude_HEADERS = # Append below

# Use MPI compilers
CC=@MPICC@
FC=@MPIFC@

# Headers installed with the library
pkginclude_HEADERS += esio.h

# Source files built into the library
lib_LTLIBRARIES     += libesio.la
libesio_la_SOURCES   = esio.c

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
