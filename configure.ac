dnl ----------------------------------------------------
dnl Package name, release number, description, and URL
dnl ----------------------------------------------------
AC_INIT([ESIO], [0.1.9], [rhys@ices.utexas.edu], [esio])
AC_REVISION([$Id$])
PACKAGE_DESCRIPTION="ExaScale IO library for turbulence simulation restart files"
AC_SUBST([PACKAGE_DESCRIPTION])
PACKAGE_URL="http://red.ices.utexas.edu/projects/esio/"
AC_SUBST([PACKAGE_URL])

dnl ------------------------------------------------
dnl Initialization macros
dnl ------------------------------------------------
AC_CONFIG_COMMANDS_PRE([m4_ifdef([AH_HEADER],
    [AC_SUBST([CONFIG_INCLUDE], m4_defn([AH_HEADER]))])])
AC_CONFIG_SRCDIR([esio/esio.h])
AM_CONFIG_HEADER(esio/config.h)
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AX_PREFIX_CONFIG_H([esio/esio-config.h])
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE(1.11 no-define tar-ustar color-tests -Wall foreign)
AC_PREREQ([2.64])
AX_AM_MACROS
AX_ADD_AM_CTAGSFLAGS
AX_AM_JOBSERVER([no])
AX_SPLIT_VERSION
AC_SUBST([AX_MAJOR_VERSION])
AC_SUBST([AX_MINOR_VERSION])
AC_SUBST([AX_POINT_VERSION])
AC_PROG_SED
AC_PATH_PROG([PERL],[perl],[no])
AM_CONDITIONAL([HAVE_PERL],[test x$PERL != xno])
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl Manipulate C settings and find MPI/HDF5
dnl ------------------------------------------------
AX_CHECK_ENABLE_DEBUG([info])
AC_PROG_CC
gl_EARLY
gl_INIT
AX_APPEND_LINK_FLAGS([-Wl,--as-needed],[LDFLAGS])
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
if test "x$with_hdf5" != "xyes"; then
    AC_MSG_ERROR([Parallel HDF5 installation not detected.])
fi
AX_VISIBILITY([hidden],[:],[:])
AX_TLS()
AX_COMPILER_VENDOR()
AX_WARNINGS_SANITIZE()
AC_LANG_POP([C])
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl Manipulate Fortran settings and find MPI
dnl ------------------------------------------------
AC_PROG_FC
AC_LANG_PUSH([Fortran])
dnl BEGIN: Ensure that preprocessor is run on free-form Fortran .F90 sources
dnl Stupidly complicated to get correct.  Modify FCFLAGS_F90 with result.
dnl When Autoconf 2.69 is more common, AC_FC_PP_SRCEXT([F90]) should be used.
AC_FC_SRCEXT([F90])
AC_FC_FREEFORM()
ax_save_FCFLAGS=$FCFLAGS
dnl AC_LANG_WERROR # IBM XL Fortran emits success strings on stderr and breaks us!
AC_MSG_CHECKING([for flags required to enable Fortran preprocessing])
for fppflag in '' -fpp -cpp '-x f95-cpp-input' none
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
FCFLAGS_F90="$FCFLAGS_F90 $fppflag"
dnl END: Ensure that preprocessor is run on free-form Fortran .F90 sources
AC_FC_WRAPPERS
AX_FCFLAGS_WARN_ALL
if test "x$ax_enable_debug" != "xno"; then
    AX_TRACEBACK()
fi
ACX_MPI(,AC_MSG_ERROR(Could not find MPI compilation tools for Fortran, e.g. mpf90.))
MPIFC_LOCATION=$(eval command -v $MPIFC)
AC_SUBST([MPIFC_LOCATION])
if test "x$with_hdf5_fortran" != "xyes"; then
    AC_MSG_WARN([HDF5 Fortran support not detected])
    AC_MSG_NOTICE([ESIO Fortran support does not require HDF5 Fortran support])
fi
AC_LANG_POP([Fortran])
AC_CACHE_SAVE

dnl ------------------------------------------------
dnl Enable GNU libtool
dnl Current version checks also in Makefile.am
dnl ------------------------------------------------
LT_INIT()
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
DX_MAN_FEATURE([ON])
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
AC_SUBST([CONFIG_STATUS_DEPENDENCIES], ['$(top_srcdir)/LICENSE.in'])
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
    gnulib/Makefile
    tests/Makefile
])

AC_OUTPUT()

dnl See also libtool target in toplevel Makefile.am
AC_MSG_NOTICE([Working around linking futziness by adding -Wl to libtool])
AC_MSG_NOTICE([See http://mail.hdfgroup.org/pipermail/hdf-forum_hdfgroup.org/2011-March/008868.html])
mv -v libtool libtool.bak
$SED libtool.bak -e 's/^wl=""$/wl="-Wl,"/' > libtool && rm libtool.bak
