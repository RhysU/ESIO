# SYNOPSIS
#
#   AX_WARNINGS_SANITIZE
#
# DESCRIPTION
#
#   Some compiler vendors, in particular Intel, offer no -pedantic option and
#   have -Wall complain about more than necessary.  Some compiler vendors also
#   offer diagnostic information above and beyond -Wall.  This macro enables
#   the nicer diagnostics and also disables some useless warnings via CFLAGS
#   and CXXFLAGS.
#
#   Requires macros: AX_COMPILER_VENDOR, AX_APPEND_COMPILE_FLAGS
#
# LAST MODIFICATION
#
#   2010-09-21
#
# COPYLEFT
#
#   Copyright (c) 2010 Rhys Ulerich <rhys.ulerich@gmail.com>
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2008 Matteo Frigo
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Macro Archive. When you make and
#   distribute a modified version of the Autoconf Macro, you may extend this
#   special exception to the GPL to apply to your modified version as well.

AC_DEFUN([AX_WARNINGS_SANITIZE],
[
AC_LANG_PUSH([C])
case $ax_cv_c_compiler_vendor in #(
  dec)   ;;#(
  sun)   ;;#(
  hp)    ;;#(
  ibm)   ;;#(
  intel) # Enable -Wcheck for more diagnostic information
         AX_APPEND_COMPILE_FLAGS([-Wcheck])
         # Enable diagnostics about questionable pointer arithmetic
         AX_APPEND_COMPILE_FLAGS([-Wpointer-arith])
         # Disable remark #424: extra ";" ignored
         AX_APPEND_COMPILE_FLAGS([-wd424])
         ;;#(
  gnu)   AX_APPEND_COMPILE_FLAGS([-Wextra])
         ;;#(
  *)     # Problem may occur if AX_COMPILER_VENDOR not called prior to AX_WARNINGS_SANITIZE
         AC_MSG_WARN([AX_WARNINGS[]_SANITIZE: ax_cv_c_compiler_vendor = $ax_cv_c_compiler_vendor unknown])
         ;;
esac
# Disable warnings on unrecognized pragmas for any compiler
AX_APPEND_COMPILE_FLAGS([-Wno-unknown-pragmas])
AC_LANG_POP()
AC_LANG_PUSH([C++])
case $ax_cv_cxx_compiler_vendor in #(
  dec)   ;;#(
  sun)   ;;#(
  hp)    ;;#(
  ibm)   ;;#(
  intel) # Enable -Wcheck for more diagnostic information
         AX_APPEND_COMPILE_FLAGS([-Wcheck])
         # Enable -Weffc++ for warnings related to C++ programming guidelines
dnl      AX_APPEND_COMPILE_FLAGS([-Weffc++])
dnl      # Enable diagnostics about questionable pointer arithmetic
         AX_APPEND_COMPILE_FLAGS([-Wpointer-arith])
         # Enable warnings if generated code is not C++ ABI compliant
         AX_APPEND_COMPILE_FLAGS([-Wabi])
         # Disable remark #383: value copied to temporary, reference to temporary used
         AX_APPEND_COMPILE_FLAGS([-wd383])
         # Disable remark #424: extra ";" ignored
         AX_APPEND_COMPILE_FLAGS([-wd424])
         # Disable remark #981: operands are evaluated in unspecified order
         AX_APPEND_COMPILE_FLAGS([-wd981])
dnl      # Disable warning #2012: Effective C++ Item 1 prefer const and inline to #define
dnl      AX_APPEND_COMPILE_FLAGS([-wd2012])
dnl      # Disable warning #2015: Effective C++ Item 4 prefer C++ style comments
dnl      AX_APPEND_COMPILE_FLAGS([-wd2015])
         ;;#(
  gnu)   AX_APPEND_COMPILE_FLAGS([-Wextra])
         ;;#(
  *)     # Problem may occur if AX_COMPILER_VENDOR not called prior to AX_WARNINGS_SANITIZE
         AC_MSG_WARN([AX_WARNINGS[]_SANITIZE: ax_cv_cxx_compiler_vendor = $ax_cv_cxx_compiler_vendor unknown])
         ;;
esac
# Disable warnings on unrecognized pragmas for any compiler
AX_APPEND_COMPILE_FLAGS([-Wno-unknown-pragmas])
AC_LANG_POP()
])
