# SYNOPSIS
#
#   AX_WARNINGS_SANITIZE
#
# DESCRIPTION
#
#   Some compiler vendors, in particular Intel, offer no -pedantic option and
#   have -Wall complain about more than necessary.  Some compiler vendors also
#   offer diagnostic information above and beyond -Wall.  This macro enables
#   the nicer diagnostics and also disables some useless warnings via CFLAGS.
#
#   Requires macros: AX_COMPILER_VENDOR, AX_CHECK_COMPILER_FLAGS
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
         AX_CHECK_COMPILER_FLAGS([-Wcheck], [CFLAGS="$CFLAGS -Wcheck"])
         # Enable diagnostics about questionable pointer arithmetic
         AX_CHECK_COMPILER_FLAGS([-Wpointer-arith], [CFLAGS="$CFLAGS -Wpointer-arith"])
         # Disable remark #981: operands are evaluated in unspecified order
         AX_CHECK_COMPILER_FLAGS([-wd981], [CFLAGS="$CFLAGS -wd981"])
         # Enable information about vectorization
         AX_CHECK_COMPILER_FLAGS([-vec-report1], [CFLAGS="$CFLAGS -vec-report1"])
         ;;#(
  gnu)   AX_CHECK_COMPILER_FLAGS([-Wextra], [CFLAGS="$CFLAGS -Wextra"])
         ;;#(
  *)     # Problem may occur if AX_COMPILER_VENDOR not called prior to AX_WARNINGS_SANITIZE
         AC_MSG_WARN([AX_WARNINGS[]_SANITIZE: ax_cv_c_compiler_vendor = $ax_cv_c_compiler_vendor unknown])
         ;;
esac
# Disable warnings on unrecognized pragmas for any compiler
AX_CHECK_COMPILER_FLAGS([-Wno-unknown-pragmas],
                        [CFLAGS="$CFLAGS -Wno-unknown-pragmas"])
AC_LANG_POP()
])
