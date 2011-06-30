# SYNOPSIS
#
#   AX_VISIBILITY([setting=hidden/default[,ACTION-SUCCESS[,ACTION-FAILURE]]])
#
# DESCRIPTION
#
#   Determine if the current language's compiler supports a symbol visibility
#   flag at the given setting (e.g. -fvisibility=hidden for GCC) and run the
#   appropriate action.  The default setting is 'hidden'.
#
#   On success, the macro sets VISIBILITY_[]_AC_LANG_PREFIX[]FLAGS to contain
#   the desired setting and runs ACTION-SUCCESS.  On failure, the macro runs
#   ACTION-FAILURE.  In either case, an AM_CONDITIONAL named HAVE_VISIBILITY
#   is made available to Automake to facilitate conditional build behavior.
#
#   The default ACTION-SUCCESS appends the flag to _AC_LANG_PREFIX[]FLAGS and
#   AC_DEFINEs HAVE_VISIBILITY.  The default ACTION-FAILURE does nothing.
#
# REQUIRES
#
#   AX_CHECK_COMPILE_FLAG
#
# LAST MODIFICATION
#
#   2010-12-13
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

AC_DEFUN([AX_VISIBILITY], [
AC_PREREQ(2.59) dnl for _AC_LANG_PREFIX

m4_define(ax_visibility_setting,[m4_tolower(m4_normalize(ifelse([$1],,[hidden],[$1])))])

ax_visibility_[]_AC_LANG_ABBREV[]_flag=
dnl Simulate a PUSH/POP on enabling AC_LANG_WERROR
ax_visibility_[]_AC_LANG_ABBREV[]_werror_save=$ac_[]_AC_LANG_ABBREV[]_werror_flag
ac_[]_AC_LANG_ABBREV[]_werror_flag=yes
for possibility in "-fvisibility=ax_visibility_setting"
do
  AX_CHECK_COMPILE_FLAG([$possibility],[
    ax_visibility_[]_AC_LANG_ABBREV[]_flag=$possibility
    break
  ],[])
done
ac_[]_AC_LANG_ABBREV[]_werror_flag=$ax_visibility_[]_AC_LANG_ABBREV[]_werror_save

AM_CONDITIONAL([HAVE_VISIBILITY],dnl
               [test x"$ax_visibility_[]_AC_LANG_ABBREV[]_flag" != x])

if test x"$ax_visibility_[]_AC_LANG_ABBREV[]_flag" != x; then
  VISIBILITY_[]_AC_LANG_PREFIX[]FLAGS=$ax_visibility_[]_AC_LANG_ABBREV[]_flag
  AC_SUBST(VISIBILITY_[]_AC_LANG_PREFIX[]FLAGS)
  m4_default([$2],[
    _AC_LANG_PREFIX[]FLAGS="$_AC_LANG_PREFIX[]FLAGS $VISIBILITY_[]_AC_LANG_PREFIX[]FLAGS"
  ])
else
  m4_default([$3],[:])
fi
])dnl AX_VISIBILITY
