# SYNOPSIS
#
#   AX_TRACEBACK([ACTION-SUCCESS], [ACTION-FAILURE])
#
# DESCRIPTION
#
#   Determine if the current language's compiler supports a "traceback" flag and run
#   the appropriate action.  On success, set TRACEBACK_[]_AC_LANG_PREFIX[]FLAGS
#   to the flag and run ACTION-SUCCESS.  The default ACTION-SUCCESS appends
#   the flag to _AC_LANG_PREFIX[]FLAGS.  The default ACTION-FAILURE does nothing.
#
# REQUIRES
#
#   AX_CHECK_COMPILER_FLAGS
#
# LAST MODIFICATION
#
#   2010-03-04
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

AC_DEFUN([AX_TRACEBACK], [
AC_PREREQ(2.59) dnl for _AC_LANG_PREFIX
AC_REQUIRE([AX_CHECK_COMPILER_FLAGS])

ax_traceback_[]_AC_LANG_ABBREV[]_flag=
for possibility in "-traceback"
do
  AX_CHECK_COMPILER_FLAGS([$possibility],[
    ax_traceback_[]_AC_LANG_ABBREV[]_flag=$possibility
    break
  ],[])
done

if test x"$ax_traceback_[]_AC_LANG_ABBREV[]_flag" != x; then
  TRACEBACK_[]_AC_LANG_PREFIX[]FLAGS=$ax_traceback_[]_AC_LANG_ABBREV[]_flag
  AC_SUBST(TRACEBACK_[]_AC_LANG_PREFIX[]FLAGS)
  m4_default([$1],[
    _AC_LANG_PREFIX[]FLAGS="$_AC_LANG_PREFIX[]FLAGS $TRACEBACK_[]_AC_LANG_PREFIX[]FLAGS"
  ])
else
  m4_default([$2],[:])
fi
])dnl AX_TRACEBACK
