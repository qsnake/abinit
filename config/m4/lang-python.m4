# -*- Autoconf -*-
#
# Copyright (C) 2009-2012 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Python support
#



# _ABI_CHECK_NUMPY_HEADERS()
# --------------------------
#
# Checks for the existence of NumPy headers.
#
AC_DEFUN([_ABI_CHECK_NUMPY_HEADERS],[
  dnl Init
  abi_numpy_new="no"
  abi_numpy_old="no"
  abi_save_CPPFLAGS="${CPPFLAGS}"
  CPPFLAGS="${PYTHON_CPPFLAGS} ${CPPFLAGS}"

  dnl Look for a recent implementation
  AC_MSG_CHECKING([for Python NumPy headers])
  if test "${abi_numpy_new}" = "no"; then
    AC_COMPILE_IFELSE(
      [AC_LANG_PROGRAM([[
#include <Python.h>
#include <numpy/arrayobject.h>]],[[
int main(int argc, char** argv)
{
    return 0;
}]])],[abi_numpy_new="yes"])
    if test "${abi_numpy_new}" = "yes"; then
      AC_DEFINE([HAVE_NUMPY],1,[Define to 1 if you have a modern implementation of NumPy.])
      AC_MSG_RESULT([found])
    else
      AC_MSG_RESULT([not found])
    fi
  fi
  
  if test "${abi_numpy_new}" = "no"; then
    AC_CHECK_HEADER([numarray/arrayobject.h],[abi_numpy_old="yes"])
    if test "${abi_numpy_old}" = "yes"; then
      AC_DEFINE([HAVE_NUMPY_OLD],1,[Define to 1 if you have an old implementation of NumPy.])
    fi
  fi

  dnl Restore environment
  CPPFLAGS="${abi_save_CPPFLAGS}"
]) # _ABI_CHECK_NUMPY_HEADERS



# ABI_CHECK_PYTHON()
# ------------------
#
# Checks whether the Python environment satisfies the requirements of Abinit.
#
AC_DEFUN([ABI_CHECK_PYTHON],
[
  _ABI_CHECK_NUMPY_HEADERS
]) # ABI_CHECK_PYTHON
