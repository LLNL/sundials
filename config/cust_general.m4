# This file is part of Autoconf.                       -*- Autoconf -*-
# Parameterized macros.
# Copyright (C) 1992, 1993, 1994, 1995, 1996, 1998, 1999, 2000, 2001,
# 2002, 2003, Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# As a special exception, the Free Software Foundation gives unlimited
# permission to copy, distribute and modify the configure scripts that
# are the output of Autoconf.  You need not follow the terms of the GNU
# General Public License when using or distributing such scripts, even
# though portions of the text of Autoconf appear in them.  The GNU
# General Public License (GPL) does govern all other use of the material
# that constitutes the Autoconf program.
#
# Certain portions of the Autoconf source text are designed to be copied
# (in certain cases, depending on the input) into the output of
# Autoconf.  We call these the "data" portions.  The rest of the Autoconf
# source text consists of comments plus executable code that decides which
# of the data portions to output in any given case.  We call these
# comments and executable code the "non-data" portions.  Autoconf never
# copies any of the non-data portions into its output.
#
# This special exception to the GPL applies to versions of Autoconf
# released by the Free Software Foundation.  When you make and
# distribute a modified version of Autoconf, you may extend this special
# exception to the GPL to apply to your modified version as well, *unless*
# your modified version has the potential to copy into its output some
# of the text that was the non-data portion of the version that you started
# with.  (In other words, unless your change moves or copies text from
# the non-data portions to the data portions.)  If your modification has
# such potential, you must delete any notice of this special exception
# to the GPL from your modified version.
#
# Written by David MacKenzie, with help from
# Franc,ois Pinard, Karl Berry, Richard Pixley, Ian Lance Taylor,
# Roland McGrath, Noah Friedman, david d zuhn, and many others.


# _AC_MSG_LOG_CONFTEST_GENERAL
# ----------------------------
m4_define([_AC_MSG_LOG_CONFTEST_GENERAL],
[echo "$as_me: failed program was:" >&AS_MESSAGE_LOG_FD
if test -f conftest.c ; then
  sed 's/^/| /' conftest.c >&AS_MESSAGE_LOG_FD
elif test -f conftest.cc ; then
  sed 's/^/| /' conftest.cc >&AS_MESSAGE_LOG_FD
elif test -f conftest.f ; then
  sed 's/^/| /' conftest.f >&AS_MESSAGE_LOG_FD
elif test -f conftest.${FC_SRCEXT-f} ; then
  sed 's/^/| /' conftest.${FC_SRCEXT-f} >&AS_MESSAGE_LOG_FD
fi
])


# _AC_LINKONLY_IFELSE(PROGRAM, [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# ------------------------------------------------------------------
# Try to link PROGRAM (empty).
# This macro can be used during the selection of a compiler.
m4_define([_AC_LINKONLY_IFELSE],
[m4_ifvaln([$1], [AC_LANG_CONFTEST([$1])])dnl
AS_IF([_AC_EVAL_STDERR($ac_linkonly) &&
	 AC_TRY_COMMAND([test -z "$ac_[]_AC_LANG_ABBREV[]_werror_flag"
			 || test ! -s conftest.err]) &&
	 AC_TRY_COMMAND([test -s conftest$ac_exeext])],
      [$2],
      [_AC_MSG_LOG_CONFTEST_GENERAL
m4_ifvaln([$3], [$3])dnl])[]dnl
rm -f conftest.err conftest.$ac_objext \
      conftest$ac_exeext m4_ifval([$1], [conftest.$ac_ext])[]dnl
])# _AC_LINKONLY_IFELSE


# AC_LINKONLY_IFELSE(PROGRAM, [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# -----------------------------------------------------------------
# Try to link PROGRAM.  Requires that the compiler for the current
# language was checked for, hence do not use this macro in macros looking
# for a compiler.
AC_DEFUN([AC_LINKONLY_IFELSE],
[AC_LANG_COMPILER_REQUIRE()dnl
_AC_LINKONLY_IFELSE($@)])
