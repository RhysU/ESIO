dnl ----------------------------------------------------
dnl Package name and release number
dnl ---------------------------------------------------
AC_INIT(esio, [0.0.1], [rhys@ices.utexas.edu])
AC_REVISION([$Id$])

dnl ------------------------------------------------
dnl Initialization macros
dnl ------------------------------------------------
AC_CONFIG_SRCDIR([esio/esio.h])
AM_CONFIG_HEADER(esio/config.h)
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE(1.9 no-define tar-ustar -Wall -Werror foreign)
AC_PREREQ([2.61])
AX_AM_MACROS
AX_ADD_AM_CTAGSFLAGS
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl Check for compiler toolchain
dnl ------------------------------------------------
AX_CHECK_ENABLE_DEBUG
AC_PROG_CC
AC_PROG_FC
AC_FC_WRAPPERS
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl Enable GNU libtool
dnl Current version checks also in Makefile.am
dnl ------------------------------------------------
LT_INIT
LT_PREREQ([2.2.6])
AC_SUBST(LIBTOOL_DEPS)
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl C strictness checks
dnl ------------------------------------------------
AC_LANG([C])
AX_CFLAGS_WARN_ALL_ANSI
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl HPC Toolkit checks
dnl ------------------------------------------------
AC_LANG([C])
PKG_CHECK_MODULES(GRVY, grvy >= 0.27.0)
AC_SUBST(GRVY_CFLAGS)
AC_SUBST(GRVY_LIBS)
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl MPI toolchain checks
dnl ------------------------------------------------
AC_LANG([C])
ACX_MPI(,AC_MSG_ERROR(Could not find MPI compilation tools for C, e.g. mpicc.))
AC_LANG([Fortran])
ACX_MPI(,AC_MSG_ERROR(Could not find MPI compilation tools for Fortran, e.g. mpif90.))
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl HDF5 checks
dnl ------------------------------------------------
AC_LANG([C])
AX_LIB_HDF5([parallel])
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl Documentation generation
dnl ------------------------------------------------
DX_DOXYGEN_FEATURE([ON])
DX_HTML_FEATURE([ON])
DX_CHM_FEATURE([OFF])
DX_CHI_FEATURE([OFF])
DX_MAN_FEATURE([OFF])
DX_RTF_FEATURE([OFF])
DX_XML_FEATURE([OFF])
DX_PDF_FEATURE([OFF])
DX_PS_FEATURE([OFF])
DX_INIT_DOXYGEN([esio],[doxygen.cfg],[docs])

dnl ------------------------------------------------
dnl Generate Makefiles
dnl ------------------------------------------------
AC_CONFIG_FILES([
    Makefile
    apps/Makefile
    lib/Makefile
    esio/Makefile
    tests/Makefile
    esio.pc
])

AC_OUTPUT()
