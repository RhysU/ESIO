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

# Use MPI compilers
CC=@MPICC@
FC=@MPIFC@

# Sample progams
bin_PROGRAMS       += hellompi_c
hellompi_c_SOURCES  = hellompi_c.c

bin_PROGRAMS              += hellompi_f90
hellompi_f90_SOURCES       = hellompi_f90.f90
hellompi_f90_LIBTOOLFLAGS  = --tag=FC

# test program
bin_PROGRAMS       += npmtest
npmtest_SOURCES   = npmtest.c

# Nick's test program
bin_PROGRAMS      += main
main_SOURCES       = header.f90 main.f90
main_LDADD         = ../esio/libesio.la
main_LIBTOOLFLAGS  = --tag=FC

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
