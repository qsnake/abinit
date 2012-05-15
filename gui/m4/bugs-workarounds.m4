# -*- Autoconf -*-
#
# Copyright (C) 2006-2012 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# File I/O
#



# ABI_PROG_MKDIR_P()
# ------------------
#
# Wrapper for the bugged AC_PROG_MKDIR_P macro.
#
AC_DEFUN([ABI_PROG_MKDIR_P],[
  _AC_SRCDIRS(["."])
  AC_PROG_MKDIR_P
  abi_tmp_mkdir_p=`echo "${MKDIR_P}" | awk '{print [$]1}'`
  if test "${abi_tmp_mkdir_p}" = "gnu/install-sh"; then
    AC_MSG_NOTICE([fixing wrong path to mkdir replacement])
    MKDIR_P="${ac_abs_top_srcdir}/${MKDIR_P}"
  fi
  unset abi_tmp_mkdir_p
]) # ABI_PROG_MKDIR_P
