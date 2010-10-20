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
AX_SPLIT_VERSION
AC_SUBST([AX_MAJOR_VERSION])
AC_SUBST([AX_MINOR_VERSION])
AC_SUBST([AX_POINT_VERSION])
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl Manipulate C compiler settings
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
AC_LANG_POP([C])
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl Look for HDF5 installation and compiler wrappers
dnl ------------------------------------------------
AC_LANG([C])
AX_LIB_HDF5([parallel])
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl Manipulate Fortran compiler settings
dnl ------------------------------------------------
AC_PROG_FC
AC_LANG_PUSH([Fortran])
if test "x$ax_enable_debug" != "xno"; then
    AX_TRACEBACK()
fi
AC_FC_WRAPPERS
AX_FCFLAGS_WARN_ALL
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
DX_INIT_DOXYGEN([esio],[doxygen.cfg],[docs])

#--------------------------
# Checks for code coverage
#-------------------------
AX_CODE_COVERAGE

dnl ------------------------------------------------
dnl Generate Makefiles
dnl ------------------------------------------------
AC_CONFIG_FILES([
    Makefile
    apps/Makefile
    lib/Makefile
    esio/Makefile
    esio/version.h
    tests/Makefile
    esio.pc
])

AC_OUTPUT()
