## Process this file with automake to produce Makefile.in

AM_CFLAGS          = # Append below
AM_CPPFLAGS        = -I$(top_builddir)/gnulib -I$(top_srcdir)/gnulib
AM_CPPFLAGS       += -I$(top_builddir)        -I$(top_srcdir)
AM_FCFLAGS         = -I$(top_builddir)        -I$(top_srcdir)
AM_LDFLAGS         = # Append below
BUILT_SOURCES      = # Append below
check_PROGRAMS     = # Append below
CLEANFILES         = # Append below
EXTRA_DIST         = # Append below
include_HEADERS    = # Append below
lib_LTLIBRARIES    = # Append below
noinst_HEADERS     = # Append below
noinst_LTLIBRARIES = # Append below
moduledir          = @libdir@

# Use MPI compiler wrappers since HDF5 wrappers break autodependencies.
CC=@MPICC@
FC=@MPIFC@

# Make HDF5 C APIs visible during compilation
# Notice HDF5 Fortran APIs are not used and so AM_FCFLAGS is not modified
AM_CFLAGS   += $(HDF5_CFLAGS)
AM_CPPFLAGS += $(HDF5_CPPFLAGS)

# Make public C APIs visible in resulting libraries-- see visibility.h
AM_CPPFLAGS += -DESIO_SHARED_EXPORTS

# Internal ESIO C library consumed by both the public C and Fortran libraries
noinst_LTLIBRARIES                += libesio_internal.la
libesio_internal_la_SOURCES        = chunksize.c      chunksize.h
libesio_internal_la_SOURCES       += error.c          error.h
libesio_internal_la_SOURCES       += esio.c           esio.h
libesio_internal_la_SOURCES       += file-copy.c      file-copy.h
libesio_internal_la_SOURCES       += h5utils.c        h5utils.h
libesio_internal_la_SOURCES       += layout.c         layout.h
libesio_internal_la_SOURCES       += metadata.c       metadata.h
libesio_internal_la_SOURCES       += restart-rename.c restart-rename.h
libesio_internal_la_SOURCES       += uri.c            uri.h
libesio_internal_la_LDFLAGS        = $(LTLIBINTL)         # gnulib requirement
libesio_internal_la_LIBADD         = ../gnulib/libgnu.la  # gnulib requirement
EXTRA_libesio_internal_la_SOURCES  = layout0.c  # x-macro template
EXTRA_libesio_internal_la_SOURCES += layout1.c  # x-macro template
EXTRA_libesio_internal_la_SOURCES += layout2.c  # x-macro template
EXTRA_libesio_internal_la_SOURCES += line.c     # x-macro template
EXTRA_libesio_internal_la_SOURCES += plane.c    # x-macro template

# Public ESIO C library
lib_LTLIBRARIES    += libesio.la
libesio_la_SOURCES  =
libesio_la_LDFLAGS  = -release $(PACKAGE_VERSION)
libesio_la_LIBADD   = libesio_internal.la

# Public ESIO C headers
pkginclude_HEADERS  = error.h
pkginclude_HEADERS += esio.h
pkginclude_HEADERS += version.h
pkginclude_HEADERS += visibility.h
pkginclude_HEADERS += esio-config.h
BUILT_SOURCES      += esio-config.h # See AX_PREFIX_CONFIG_H

# Public ESIO Fortran library and module files
lib_LTLIBRARIES           += libesiof.la
libesiof_la_LDFLAGS        = -release $(PACKAGE_VERSION)
libesiof_la_LIBTOOLFLAGS   = --tag=FC
libesiof_la_SOURCES        = esio_c_binding.f90 esio.f90
module_DATA                = esio_c_binding.mod esio.mod
libesiof_la_LIBADD         = libesio_internal.la
EXTRA_libesiof_la_SOURCES  = attribute.f90   # x-macro template
EXTRA_libesiof_la_SOURCES += field.f90       # x-macro template
EXTRA_libesiof_la_SOURCES += line.f90        # x-macro template
EXTRA_libesiof_la_SOURCES += plane.f90       # x-macro template

# Extra things to clean up
CLEANFILES += *.mod *.gcno *.gcna

# License updates
include $(top_srcdir)/license.am

# Automatic dependency generation for Fortran
SFMAKEDEPEND_FLAGS =
include $(top_srcdir)/sfmakedepend.am

# Required for AX_AM_MACROS
@INC_AMINCLUDE@
