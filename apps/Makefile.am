## Process this file with automake to produce Makefile.in

AM_CFLAGS       =
AM_CPPFLAGS     = -I$(top_builddir)/gnulib -I$(top_srcdir)/gnulib
AM_CPPFLAGS    += -I$(top_builddir)        -I$(top_srcdir)
AM_FCFLAGS      = -I$(top_builddir)        -I$(top_srcdir)
AM_FCFLAGS     += -I$(top_builddir)/esio
AM_LDFLAGS      = # Append below
noinst_PROGRAMS = # Append below
bin_PROGRAMS    = # Append below
BUILT_DOCS      = # Append below
BUILT_SOURCES   = # Append below
CLEANFILES      = # Append below
EXTRA_DIST      = # Append below

# Use MPI compiler wrappers since HDF5 wrappers break autodependencies.
CC=@MPICC@
FC=@MPIFC@

## Pick up GRVY library includes
AM_CPPFLAGS += $(GRVY_CFLAGS)
AM_FCFLAGS  += $(GRVY_FCFLAGS)

# C-based benchmarking tool
bin_PROGRAMS        += esio_bench
esio_bench_SOURCES   = esio_bench.c mpi_argp.c mpi_argp.h
esio_bench_CPPFLAGS  = $(AM_CPPFLAGS) $(HDF5_CPPFLAGS)
esio_bench_LDADD     = ../esio/libesio.la
esio_bench_LDADD    += ../gnulib/libgnu.la $(LTLIBINTL)
esio_bench_LDADD    += $(GRVY_LIBS)

# Restart renaming tool which has few formal ESIO dependencies
bin_PROGRAMS        += esio_rename
esio_rename_SOURCES  = esio_rename.c ../esio/restart-rename.c ../esio/error.c
esio_rename_LDADD    = ../gnulib/libgnu.la $(LTLIBINTL) -lm

# PSDNS iokernel test program based on P3DFFT v2.3 loading patterns
# See http://code.google.com/p/p3dfft/ for more details on P3DFFT
noinst_PROGRAMS          += p3dfft_like
p3dfft_like_SOURCES       = p3dfft_like_header.f90 p3dfft_like.f90
p3dfft_like_LDADD         = ../esio/libesiof.la $(GRVY_LIBS)
p3dfft_like_LIBTOOLFLAGS  = --tag=FC
p3dfft_like.f90: p3dfft_like_header.o    # Explicit module dependence

# Additional files to clean
CLEANFILES += *.mod *.gcno *.gcda

# License updates
include $(top_srcdir)/license.am

# Automatic dependency generation for Fortran
SFMAKEDEPEND_FLAGS = -i esio
include $(top_srcdir)/sfmakedepend.am

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
