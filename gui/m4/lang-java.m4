# -*- Autoconf -*-
#
# Copyright (C) 2011-2012 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Java support
#



# ABI_JAVA_CHECK_OVERRIDE()
# -------------------------
#
# Check whether the current Java implementation supports @Override.
#
AC_DEFUN([ABI_JAVA_CHECK_OVERRIDE],[
  dnl Init
  abi_java_override="no"

  AC_MSG_CHECKING([whether Java supports '@Override'])
  AC_TRY_COMPILE_JAVA(,
    [
        private String field;
        private String attribute;

        @Override
        public int hashCode() {
            return field.hashCode() + attribute.hashCode();
        }

        @Override
        public String toString() {
            return field + " " + attribute;
        }
    ],[abi_java_override="yes"])
  AC_MSG_RESULT([${abi_java_override}])
]) # ABI_JAVA_CHECK_OVERRIDE
