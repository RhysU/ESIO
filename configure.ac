dnl ----------------------------------------------------
dnl Package name, release number, description, and URL
dnl ----------------------------------------------------
AC_INIT([esio], [0.0.1], [rhys@ices.utexas.edu])
AC_REVISION([$Id$])
PACKAGE_DESCRIPTION="ExaScale IO library for turbulence simulation restart files"
AC_SUBST([PACKAGE_DESCRIPTION])
PACKAGE_URL="http://pecos.ices.utexas.edu/"
AC_SUBST([PACKAGE_URL])

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
AX_AM_JOBSERVER([no])
AX_SPLIT_VERSION
AC_SUBST([AX_MAJOR_VERSION])
AC_SUBST([AX_MINOR_VERSION])
AC_SUBST([AX_POINT_VERSION])
AC_PATH_PROG([PERL],[perl],[no])
AM_CONDITIONAL([HAVE_PERL],[test x$PERL != xno])
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl Manipulate C settings and find MPI/HDF5
dnl ------------------------------------------------
AX_CHECK_ENABLE_DEBUG
AC_PROG_CC
AC_LANG_PUSH([C])
if test "x$ax_enable_debug" != "xno"; then
    AX_TRACEBACK()
fi
AX_PATH_GRVY_NEW([0.27],[no])
AX_CFLAGS_WARN_ALL
AX_PROG_CC_C99_CFLAGS(,AC_MSG_ERROR([Compiler cannot support ISO C99]))
ACX_MPI(,AC_MSG_ERROR(Could not find MPI compilation tools for C, e.g. mpicc.))
MPICC_LOCATION=$(eval command -v $MPICC)
AC_SUBST([MPICC_LOCATION])
AX_MPIP()
AX_LIB_HDF5([parallel])
AC_LANG_POP([C])
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl Manipulate Fortran settings and find MPI
dnl ------------------------------------------------
AC_PROG_FC
AC_LANG_PUSH([Fortran])
AC_FC_SRCEXT([f90])
AC_FC_FREEFORM()
dnl BEGIN: Enable the preprocessor for all Fortran sources
dnl Man this is stupidly complicated to get correct
ax_save_FCFLAGS=$FCFLAGS
AC_LANG_WERROR
AC_MSG_CHECKING([for flags required to enable Fortran preprocessing])
for fppflag in '' -fpp -cpp none
do
    if test "x$fppflag" == xnone; then
        AC_MSG_RESULT([unknown])
        AC_MSG_ERROR([Could not enable Fortran preprocessing])
        break
    fi
    ax_fpp_save_FCFLAGS=$FCFLAGS
    FCFLAGS="$FCFLAGS $fppflag"
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[#define TESTDEFINE])],[
                         FCFLAGS=$ax_fpp_save_FCFLAGS
                         break
                      ],[
                         FCFLAGS=$ax_fpp_save_FCFLAGS
                      ])
done
AC_MSG_RESULT([$fppflag])
FCFLAGS="$ax_save_FCFLAGS $fppflag"
dnl END: Enable the preprocessor for all Fortran sources
AC_FC_WRAPPERS
AX_FCFLAGS_WARN_ALL
if test "x$ax_enable_debug" != "xno"; then
    AX_TRACEBACK()
fi
ACX_MPI(,AC_MSG_ERROR(Could not find MPI compilation tools for Fortran, e.g. mpf90.))
MPIFC_LOCATION=$(eval command -v $MPIFC)
AC_SUBST([MPIFC_LOCATION])
if test "x$with_hdf5_fortran" != "xyes"; then
    dnl TODO Build C-only portions when Fortran not detected
    AC_MSG_WARN([Did not detect HDF5 Fortran support!])
fi
AC_LANG_POP([Fortran])
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
DX_INIT_DOXYGEN([esio],[doxygen/doxygen.cfg],[docs])

#--------------------------
# Checks for code coverage
#-------------------------
AX_CODE_COVERAGE

dnl ------------------------------------------------
dnl Generate Makefiles
dnl ------------------------------------------------
AC_CONFIG_FILES([
    LICENSE
    Makefile
    esio.pc
    esio-hdf5.pc
    apps/Makefile
    doxygen/Makefile
    esio/Makefile
    esio/version.h
    examples/Makefile
    lib/Makefile
    tests/Makefile
])

AC_OUTPUT()
