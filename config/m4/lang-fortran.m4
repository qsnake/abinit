# -*- Autoconf -*-
#
# Copyright (C) 2005-2012 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Fortran compilers support
#



# _ABI_CHECK_FC_ABSOFT(COMPILER)
# ------------------------------
#
# Checks whether the specified Fortran compiler is the ABSoft Fortran compiler.
# If yes, tries to determine its version number and sets the abi_fc_vendor
# and abi_fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_ABSOFT],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the ABSoft Fortran compiler])
  fc_info_string=`$1 -V 2>/dev/null | head -n 1`
  abi_result=`echo "${fc_info_string}" | grep '^Pro Fortran'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    fc_info_string=""
    abi_fc_vendor="unknown"
    abi_fc_version="unknown"
  else
    AC_DEFINE([FC_ABSOFT],1,
      [Define to 1 if you are using the ABSOFT Fortran compiler.])
    abi_fc_vendor="absoft"
    abi_fc_version=`echo "${abi_result}" | sed -e 's/Pro Fortran //'`
    if test "${abi_fc_version}" = "${abi_result}"; then
      abi_fc_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_ABSOFT



# _ABI_CHECK_FC_COMPAQ(COMPILER)
# ------------------------------
#
# Checks whether the specified Fortran compiler is the COMPAQ Fortran compiler.
# If yes, tries to determine its version number and sets the abi_fc_vendor
# and abi_fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_COMPAQ],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Compaq Fortran compiler])
  fc_info_string=`$1 -version 2>&1 | sed -e 's/^	//' | grep '^Compaq Fortran Compiler'`
  abi_result="${fc_info_string}"
  if test "${abi_result}" = ""; then
    fc_info_string=`$1 -version 2>&1 | sed -e 's/^	//' | grep '^HP Fortran Compiler'`
    abi_result="${fc_info_string}"
  fi
  if test "${abi_result}" = ""; then
    abi_result="no"
    fc_info_string=""
    abi_fc_vendor="unknown"
    abi_fc_version="unknown"
  else
    AC_DEFINE([FC_COMPAQ],1,
      [Define to 1 if you are using the COMPAQ Fortran compiler.])
    abi_fc_vendor="compaq"
    abi_fc_version=`echo "${abi_result}" | sed -e 's/.* V//;s/-.*//'`
    if test "${abi_fc_version}" = "${abi_result}"; then
      abi_fc_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_COMPAQ



# _ABI_CHECK_FC_FUJITSU(COMPILER)
# -------------------------------
#
# Checks whether the specified Fortran compiler is the Fujitsu Fortran compiler.
# If yes, tries to determine its version number and sets the abi_fc_vendor
# and abi_fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_FUJITSU],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Fujitsu Fortran compiler])
  fc_info_string=`$1 -V 2>/dev/null | head -n 1`
  abi_result=`echo "${fc_info_string}" | grep '^Fujitsu Fortran'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    fc_info_string=""
    abi_fc_vendor="unknown"
    abi_fc_version="unknown"
  else
    AC_DEFINE([FC_FUJITSU],1,
      [Define to 1 if you are using the Fujitsu Fortran compiler.])
    abi_fc_vendor="fujitsu"
    abi_fc_version=`echo "${abi_result}" | sed -e 's/.*Driver //;s/ .*//'`
    if test "${abi_fc_version}" = "${abi_result}"; then
      abi_fc_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_FUJITSU



# _ABI_CHECK_FC_G95(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the G95 Fortran compiler.
# If yes, tries to determine its version number and sets the abi_fc_vendor
# and abi_fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_G95],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the G95 Fortran compiler])
  fc_info_string=`$1 --version 2>/dev/null | head -n 1`
  abi_result=`echo "${fc_info_string}" | grep '^G95'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    fc_info_string=""
    abi_fc_vendor="unknown"
    abi_fc_version="unknown"
  else
    AC_DEFINE([FC_G95],1,
      [Define to 1 if you are using the G95 Fortran compiler.])
    abi_fc_vendor="g95"
    abi_fc_version=`echo ${abi_result} | sed -e 's/.*GCC //; s/ .*//'`
    if test "${abi_fc_version}" = "${abi_result}"; then
      abi_fc_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_G95



# _ABI_CHECK_FC_GNU(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the GNU Fortran compiler.
# If yes, tries to determine its version number and sets the abi_fc_vendor
# and abi_fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_GNU],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the GNU Fortran compiler])
  fc_info_string=`$1 --version 2>/dev/null | head -n 1`
  abi_result=`echo "${fc_info_string}" | grep '^GNU Fortran'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    fc_info_string=""
    abi_fc_vendor="unknown"
    abi_fc_version="unknown"
  else
    AC_DEFINE([FC_GNU],1,
      [Define to 1 if you are using the GNU Fortran compiler.])
    AC_DEFINE([HAVE_FORTRAN2003],1,
      [Define to 1 if your Fortran compiler supports Fortran 2003.])
    abi_fc_vendor="gnu"
    abi_fc_version=`echo ${abi_result} | sed -e 's/^[[^(]]*([[^)]]*) //; s/ .*//'`
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_GNU



# _ABI_CHECK_FC_HITACHI(COMPILER)
# -------------------------------
#
# Checks whether the specified Fortran compiler is the Hitachi Fortran compiler.
# If yes, tries to determine its version number and sets the abi_fc_vendor
# and abi_fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_HITACHI],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Hitachi Fortran compiler])
  fc_info_string=`$1 -V 2>/dev/null | head -n 1`
  abi_result=`echo "${fc_info_string}" | grep '^Hitachi Fortran'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    fc_info_string=""
    abi_fc_vendor="unknown"
    abi_fc_version="unknown"
  else
    AC_DEFINE([FC_HITACHI],1,
      [Define to 1 if you are using the Hitachi Fortran compiler.])
    abi_fc_vendor="hitachi"
    abi_fc_version=`echo "${abi_result}" | sed -e 's/.*Driver //;s/ .*//'`
    if test "${abi_fc_version}" = "${abi_result}"; then
      abi_fc_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_HITACHI



# _ABI_CHECK_FC_IBM(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the IBM XL Fortran compiler.
# If yes, tries to determine its version number and sets the abi_fc_vendor
# and abi_fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_IBM],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the IBM XL Fortran compiler])
  fc_info_string=`$1 -qversion 2>&1 | head -n 1`
  fc_garbage=`$1 -qversion 2>&1 | wc -l | sed -e 's/ //g'`
  abi_result=`echo "${fc_info_string}" | grep 'IBM XL Fortran'`
  if test "${abi_result}" = ""; then
    abi_result=`echo "${fc_info_string}" | grep 'IBM(R) XL Fortran'`
  fi
  if test "${abi_result}" = ""; then
    abi_result="no"
    fc_info_string=""
    abi_fc_vendor="unknown"
    abi_fc_version="unknown"
    if test "${fc_garbage}" -gt 50; then
      AC_DEFINE([FC_IBM],1,
        [Define to 1 if you are using the IBM XL Fortran compiler.])
      abi_fc_vendor="ibm"
      abi_fc_version="unknown"
      abi_result="yes"
    fi
  else
    AC_DEFINE([FC_IBM],1,
      [Define to 1 if you are using the IBM XL Fortran compiler.])
    abi_fc_vendor="ibm"
    abi_fc_version=`echo "${abi_result}" | sed -e 's/.* V//; s/ .*//'`
    if test "${abi_fc_version}" = "${abi_result}"; then
      abi_fc_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_IBM



# _ABI_CHECK_FC_INTEL(COMPILER)
# -----------------------------
#
# Checks whether the specified Fortran compiler is the Intel Fortran compiler.
# If yes, tries to determine its version number and sets the abi_fc_vendor
# and abi_fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_INTEL],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Intel Fortran compiler])
  fc_info_string=`$1 -V 2>&1 | head -n 1`
  abi_result=`echo "${fc_info_string}" | grep '^Intel(R) Fortran'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    fc_info_string=""
    abi_fc_vendor="unknown"
    abi_fc_version="unknown"
  else
    AC_DEFINE([FC_INTEL],1,
      [Define to 1 if you are using the Intel Fortran compiler.])
    abi_fc_vendor="intel"
    abi_fc_version=`echo "${fc_info_string}" | sed -e 's/.*Version //;s/ .*//'`
    if test "${abi_fc_version}" = ""; then
      abi_fc_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_INTEL



# _ABI_CHECK_FC_MIPSPRO(COMPILER)
# -------------------------------
#
# Checks whether the specified Fortran compiler is the MIPSpro Fortran
# compiler.
# If yes, tries to determine its version number and sets the abi_fc_vendor
# and abi_fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_MIPSPRO],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the MIPSpro Fortran compiler])
  fc_info_string=`$1 -version 2>&1 | head -n 1`
  abi_result=`echo "${fc_info_string}" | grep '^MIPSpro'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    fc_info_string=""
    abi_fc_vendor="unknown"
    abi_fc_version="unknown"
  else
    AC_DEFINE([FC_MIPSPRO],1,
      [Define to 1 if you are using the MIPSpro Fortran compiler.])
    abi_fc_vendor="mipspro"
    abi_fc_version=`echo "${abi_result}" | sed -e 's/.*Version //'`
    if test "${abi_fc_version}" = "${abi_result}"; then
      abi_fc_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_MIPSPRO



# _ABI_CHECK_FC_NAG(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the NAGWare Fortran 95
# compiler. If yes, tries to determine its version number and sets the
# abi_fc_vendor and abi_fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_NAG],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the NAGWare Fortran 95 compiler])
  fc_info_string=`$1 -V 2>&1 | head -n 1`
  abi_result=`echo "${fc_info_string}" | grep '^NAG'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    fc_info_string=""
    abi_fc_vendor="unknown"
    abi_fc_version="unknown"
  else
    AC_DEFINE([FC_NAG],1,
      [Define to 1 if you are using the NAGWare Fortran 95 compiler.])
    abi_fc_vendor="nag"
    abi_fc_version=`echo "${fc_info_string}" | sed -e 's/.*Release //;s/[[( ]].*//'`
    if test "${abi_fc_version}" = ""; then
      abi_fc_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_NAG



# _ABI_CHECK_FC_OPEN64(COMPILER)
# ------------------------------
#
# Checks whether the specified Fortran compiler is the Open64
# Fortran compiler.
# If yes, tries to determine its version number and sets the abi_fc_vendor
# and abi_fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_OPEN64],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Open64 Fortran compiler])
  fc_info_string=`$1 --version 2>&1 | head -n 1`
  abi_result=`echo "${fc_info_string}" | grep '^Open64'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    fc_info_string=""
    abi_fc_vendor="unknown"
    abi_fc_version="unknown"
  else
    AC_DEFINE([FC_OPEN64],1,
      [Define to 1 if you are using the Open64 Fortran compiler.])
    abi_fc_vendor="open64"
    abi_fc_version=`echo "${abi_result}" | sed -e 's/.* Version //; s/ .*//'`
    if test "${abi_fc_version}" = "${abi_result}"; then
      abi_fc_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_OPEN64



# _ABI_CHECK_FC_PATHSCALE(COMPILER)
# ---------------------------------
#
# Checks whether the specified Fortran compiler is the PathScale
# Fortran compiler.
# If yes, tries to determine its version number and sets the abi_fc_vendor
# and abi_fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_PATHSCALE],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the PathScale Fortran compiler])
  fc_info_string=`$1 -version 2>&1 | head -n 1`
  abi_result=`echo "${fc_info_string}" | grep '^PathScale'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    fc_info_string=""
    abi_fc_vendor="unknown"
    abi_fc_version="unknown"
  else
    AC_DEFINE([FC_PATHSCALE],1,
      [Define to 1 if you are using the PathScale Fortran compiler.])
    abi_fc_vendor="pathscale"
    abi_fc_version=`echo "${abi_result}" | sed -e 's/.* Version //; s/ .*//'`
    if test "${abi_fc_version}" = "${abi_result}"; then
      abi_fc_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_PATHSCALE



# _ABI_CHECK_FC_PGI(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the Portland Group
# Fortran compiler.
# If yes, tries to determine its version number and sets the abi_fc_vendor
# and abi_fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_PGI],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Portland Group Fortran compiler])
  fc_info_string=`$1 -V 2>&1 | head -n 1`
  abi_result=`echo "${fc_info_string}" | grep '^pgf9[[05]]'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    fc_info_string=""
    abi_fc_vendor="unknown"
    abi_fc_version="unknown"
  else
    AC_DEFINE([FC_PGI],1,
      [Define to 1 if you are using the Portland Group Fortran compiler.])
    abi_fc_vendor="pgi"
    abi_fc_version=`echo "${abi_result}" | sed -e 's/^pgf9[[05]] //' | sed -e 's/-.*//'`
    if test "${abi_fc_version}" = "${abi_result}"; then
      abi_fc_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_PGI



# _ABI_CHECK_FC_SUN(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the Sun WorkShop Fortran compiler.
# If yes, tries to determine its version number and sets the abi_fc_vendor
# and abi_fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_SUN],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Sun Fortran compiler])
  fc_info_string=`$1 -V 2>&1 | head -n 1`
  abi_result=`echo "${fc_info_string}" | grep 'Sun' | grep 'Fortran 95'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    fc_info_string=""
    abi_fc_vendor="unknown"
    abi_fc_version="unknown"
  else
    AC_DEFINE([FC_SUN],1,
      [Define to 1 if you are using the Sun Fortran compiler.])
    abi_fc_vendor="sun"
    abi_fc_version=`echo "${abi_result}" | sed -e 's/.* Fortran 95 //;s/ .*//'`
    if test "${abi_fc_version}" = "${abi_result}" -o "${abi_fc_version}" = ""; then
      abi_fc_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_SUN



 ##############################################################################



# _ABI_CHECK_FC_EXIT()
# --------------------
#
# Checks whether the Fortran compiler supports the exit() subroutine.
#
AC_DEFUN([_ABI_CHECK_FC_EXIT],[
  dnl Init
  fc_has_exit="no"

  AC_MSG_CHECKING([whether the Fortran compiler accepts exit()])

  dnl Try to compile a program calling exit()
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
            call exit(1)
    ]])], [fc_has_exit="yes"])
  AC_LANG_POP()

  if test "${fc_has_exit}" = "yes"; then
    AC_DEFINE([HAVE_FC_EXIT],1,
      [Define to 1 if your Fortran compiler supports exit().])
  fi

  AC_MSG_RESULT(${fc_has_exit})
]) # _ABI_CHECK_FC_EXIT



# _ABI_CHECK_FC_FLUSH()
# ---------------------
#
# Checks whether the Fortran compiler supports the flush() subroutine.
#
AC_DEFUN([_ABI_CHECK_FC_FLUSH],[
  dnl Init
  fc_has_flush="no"

  AC_MSG_CHECKING([whether the Fortran compiler accepts flush()])

  dnl Try to compile a program calling flush()
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
            call flush()
    ]])], [fc_has_flush="yes"])
  AC_LANG_POP()

  if test "${fc_has_flush}" = "yes"; then
    AC_DEFINE([HAVE_FC_FLUSH],1,
      [Define to 1 if your Fortran compiler supports flush().])
  fi

  AC_MSG_RESULT(${fc_has_flush})
]) # _ABI_CHECK_FC_FLUSH


# _ABI_CHECK_FC_FLUSH_()
# ----------------------
#   
# Checks whether the Fortran compiler supports the flush_() subroutine.
# 
AC_DEFUN([_ABI_CHECK_FC_FLUSH_],[
  dnl Init
  fc_has_flush_="no"

  AC_MSG_CHECKING([whether the Fortran compiler accepts flush_()])

  dnl Try to compile a program calling flush_()
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
            call flush_()
    ]])], [fc_has_flush_="yes"])
  AC_LANG_POP()

  if test "${fc_has_flush_}" = "yes"; then
    AC_DEFINE([HAVE_FC_FLUSH_],1,
      [Define to 1 if your Fortran compiler supports flush_().])
  fi

  AC_MSG_RESULT(${fc_has_flush_})
]) # _ABI_CHECK_FC_FLUSH_


# _ABI_CHECK_FC_GAMMA()
# ---------------------
#
# Checks whether the Fortran compiler supports the gamma() intrinsic
# (Fortran 2003 and later).
#
AC_DEFUN([_ABI_CHECK_FC_GAMMA],[
  dnl Init
  fc_has_gamma="no"

  AC_MSG_CHECKING([whether the Fortran compiler accepts gamma()])

  dnl Try to compile a program using gamma()
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
            real :: x
            x = gamma(1.5)
    ]])], [fc_has_gamma="yes"])
  AC_LANG_POP()

  if test "${fc_has_gamma}" = "yes"; then
    AC_DEFINE([HAVE_FC_GAMMA],1,
      [Define to 1 if your Fortran compiler supports gamma().])
  fi

  AC_MSG_RESULT(${fc_has_gamma})
]) # _ABI_CHECK_FC_GAMMA



# _ABI_CHECK_FC_GETENV()
# ----------------------
#
# Checks whether the Fortran compiler supports GET_ENVIRONMENT_VARIABLE
# (Fortran 2003 and later).
#
AC_DEFUN([_ABI_CHECK_FC_GETENV],[
  dnl Init
  fc_has_getenv="no"

  AC_MSG_CHECKING([whether the Fortran compiler accepts getenv()])

  dnl Try to compile a call to getenv
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], 
    [[
      character(len=255) :: homedir
      call getenv("HOME", homedir)

    ]])], [fc_has_getenv="yes"])
  AC_LANG_POP()

  if test "${fc_has_getenv}" = "yes"; then
    AC_DEFINE([HAVE_FC_GETENV],1, 
      [Define to 1 if your Fortran compiler supports getenv().])
  fi

  AC_MSG_RESULT(${fc_has_getenv})
]) # _ABI_CHECK_FC_GETENV



# _ABI_CHECK_FC_INT_QUAD()
# ------------------------
#
# Checks whether the Fortran compiler supports quadruple integers.
#
AC_DEFUN([_ABI_CHECK_FC_INT_QUAD],[
  dnl Init
  fc_has_int_quad="no"

  AC_MSG_CHECKING([whether the Fortran compiler accepts quadruple integers])

  dnl Try to compile a program defining a quadruple integer
  dnl Note: xlf "works around" the problem by changing the integer length
  if test "${abi_fc_vendor}" != "ibm"; then
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
              integer*16 my_int
      ]])], [fc_has_int_quad="yes"])
    AC_LANG_POP()
  fi

  if test "${fc_has_int_quad}" = "yes"; then
    AC_DEFINE([HAVE_FC_INT_QUAD],1,
      [Define to 1 if your Fortran compiler accepts quadruple integers.])
  fi

  AC_MSG_RESULT(${fc_has_int_quad})
]) # _ABI_CHECK_FC_INT_QUAD


# _ABI_CHECK_FC_DTARRAYS()
# --------------------
#
# Checks whether the Fortran compiler supports allocatable arrays in Fortran datatypes.
#
AC_DEFUN([_ABI_CHECK_FC_DTARRAYS],[
  dnl Init
  fc_has_dtarrays="no"

  AC_MSG_CHECKING([whether the Fortran compiler supports allocatable arrays in datatypes])

  dnl Try to compile a type with an allocatable array
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], 
    [[

         integer, parameter :: dp=kind(1.0d0)
         integer, parameter :: dpc=kind((1.0_dp,1.0_dp))  

         type test_type
           integer,allocatable :: i(:) 
           real(dp),allocatable :: r(:,:) 
           complex(dpc),allocatable :: c(:,:,:) 
         end type test_type

    ]])], [fc_has_dtarrays="yes"])
  AC_LANG_POP()

  if test "${fc_has_dtarrays}" = "yes"; then
    AC_DEFINE([HAVE_FC_ALLOCATABLE_DTARRAYS],1, 
      [Define to 1 if your Fortran compiler supports allocatable arrays in datatypes.])
  fi

  AC_MSG_RESULT(${fc_has_dtarrays})
]) # _ABI_CHECK_FC_DTARRAYS


# _ABI_CHECK_FC_ISO_C_BINDING()
# -----------------------------
#
# Checks whether the Fortran compiler provides the intrinsic module ISO_C_BINDING.
#
AC_DEFUN([_ABI_CHECK_FC_ISO_C_BINDING],[
  dnl Init
  fc_has_iso_c_binding="no"

  AC_MSG_CHECKING([whether the Fortran compiler provides the iso_c_binding module])

  dnl Try to compile a simple piece of code using iso_c_binding
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], 
    [[
         use iso_c_binding
         implicit none
         integer(c_int) :: ii
         logical :: lbool
         type(c_ptr) :: ptr
         ptr = c_null_ptr
         lbool = c_associated(ptr) 

    ]])], [fc_has_iso_c_binding="yes"])
  AC_LANG_POP()

  if test "${fc_has_iso_c_binding}" = "yes"; then
    AC_DEFINE([HAVE_FC_ISO_C_BINDING],1, 
      [Define to 1 if your Fortran compiler provides the iso_c_binding module.])
  fi

  AC_MSG_RESULT(${fc_has_iso_c_binding})
]) # _ABI_CHECK_FC_ISO_C_BINDING



# _ABI_CHECK_FC_LONG_LINES()
# --------------------------
# 
# Checks whether the Fortran compiler supports long lines.
#
AC_DEFUN([_ABI_CHECK_FC_LONG_LINES],[
  dnl Init
  fc_has_long_lines="no"

  AC_MSG_CHECKING([whether the Fortran compiler accepts long lines])

  dnl Try to compile a single line exceeding 136 columns.
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], 
    [[
         write(*,*)'0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789' !142
    ]])], [fc_has_long_lines="yes"])
  AC_LANG_POP()

  if test "${fc_has_long_lines}" = "yes"; then
    AC_DEFINE([HAVE_FC_LONG_LINES],1, 
      [Define to 1 if your Fortran compiler supports long lines.])
  fi

  AC_MSG_RESULT(${fc_has_long_lines})
]) # _ABI_CHECK_FC_LONG_LINES



# _ABI_CHECK_FC_NULL()
# --------------------
#
# Checks whether the Fortran compiler supports the null() intrinsic (particularly in type declarations).
#
AC_DEFUN([_ABI_CHECK_FC_NULL],[
  dnl Init
  fc_has_null="no"

  AC_MSG_CHECKING([whether the Fortran compiler accepts the null() intrinsic])

  dnl Try to compile a type with a pointer set to null()
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], 
    [[

         type test_type
          integer,pointer :: i(:,:) => null() 
         end type test_type

    ]])], [fc_has_null="yes"])
  AC_LANG_POP()

  if test "${fc_has_null}" = "yes"; then
    AC_DEFINE([HAVE_FC_NULL],1, 
      [Define to 1 if your Fortran compiler supports the null() intrinsic.])
  fi

  AC_MSG_RESULT(${fc_has_null})
]) # _ABI_CHECK_FC_NULL


# _ABI_CHECK_FC_STREAM_IO()
# --------------------
#
# Checks whether the Fortran compiler supports stream IO.
#
AC_DEFUN([_ABI_CHECK_FC_STREAM_IO],[
  dnl Init
  fc_has_stream_io="no"

  AC_MSG_CHECKING([whether the Fortran compiler supports stream IO])

  dnl Try to compile a piece of code that opens a file using unformatted stream access.
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], 
    [[
         IMPLICIT NONE
         INTEGER :: myvalue = 12345, mypos
         OPEN(UNIT=11, FILE="ustream.demo", STATUS="NEW", ACCESS="STREAM")
         WRITE(11) "first"
         WRITE(11) "second"
         INQUIRE(UNIT=11, POS=mypos)
         PRINT *, "Myvalue will be written at position ", mypos
         WRITE(11) myvalue
         CLOSE(UNIT=11)

    ]])], [fc_has_stream_io="yes"])
  AC_LANG_POP()

  if test "${fc_has_stream_io}" = "yes"; then
    AC_DEFINE([HAVE_FC_STREAM_IO],1, 
      [Define to 1 if your Fortran compiler supports stream IO.])
  fi

  AC_MSG_RESULT(${fc_has_stream_io})
]) # _ABI_CHECK_FC_STREAM_IO


# _ABI_CHECK_FC_TIMING()
# ----------------------
#
# Tries to determine which Fortran timing routines are available.
#
AC_DEFUN([_ABI_CHECK_FC_TIMING],[
  dnl Init
  fc_timing="standard"
  fc_has_etime="no"

  dnl Look for etime() support
  if test "${fc_timing}" = "standard"; then
    AC_MSG_CHECKING([whether the Fortran compiler accepts etime()])

    dnl Try to compile a program calling etime()
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
            call etime(1)
      ]])], [fc_has_etime="yes"])
    AC_LANG_POP()

    if test "${fc_has_etime}" = "yes"; then
      AC_DEFINE([HAVE_FC_ETIME],1,
        [Define to 1 if your Fortran compiler supports etime().])
    fi

    AC_MSG_RESULT(${fc_has_etime})

  fi

  dnl Determine whether to use C clock for timings
  AC_MSG_CHECKING([whether to use C clock for timings])
  AC_MSG_RESULT([${enable_cclock}])
  if test "${enable_cclock}" = "yes"; then
    AC_DEFINE([HAVE_CCLOCK],1,[Use C clock for timings.])
    fc_timing="cclock"
  fi
  AM_CONDITIONAL(DO_BUILD_CCLOCK,[test "${enable_cclock}" = "yes"])

  dnl Schedule info for substitution
  AC_SUBST(fc_timing)
]) # _ABI_CHECK_FC_TIMING


# _ABI_CHECK_FC_CPUTIME()
# ----------------------
#
# Checks whether the Fortran compiler supports CPU_TIME 
# (Fortran 95 and later).
#
AC_DEFUN([_ABI_CHECK_FC_CPUTIME],[
  dnl Init
  fc_has_cputime="no"

  AC_MSG_CHECKING([whether the Fortran compiler accepts cpu_time()])

  dnl Try to compile a call to cpu_time
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], 
    [[
      real :: second
      call cpu_time(second)

    ]])], [fc_has_cputime="yes"])
  AC_LANG_POP()

  if test "${fc_has_cputime}" = "yes"; then
    AC_DEFINE([HAVE_FC_CPUTIME],1, 
      [Define to 1 if your Fortran compiler supports cpu_time().])
  fi

  AC_MSG_RESULT(${fc_has_cputime})
]) # _ABI_CHECK_FC_CPUTIME


# _ABI_CHECK_FC_GETPID()
# ----------------------
#
# Checks whether process IDs are available from Fortran.
#
AC_DEFUN([_ABI_CHECK_FC_GETPID],[
  dnl Init
  fc_has_getpid="no"

  dnl Look for getpid() support
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
          call getpid()
    ]])], [fc_has_getpid="yes"])
  AC_LANG_POP([Fortran])

  dnl Determine whether to use getpid()
  AC_MSG_CHECKING([whether the Fortran compiler accepts getpid()])
  AC_MSG_RESULT([${fc_has_getpid}])
  if test "${fc_has_getpid}" = "yes"; then
    AC_DEFINE([HAVE_FC_GETPID],1,
      [Define to 1 if your Fortran compiler supports getpid().])
  fi
]) # _ABI_CHECK_FC_GETPID



 #############################################################################



# ABI_FC_EXTENSIONS()
# -------------------
#
# Sets the default extensions of Fortran source files and modules,
# whenever possible.
#
AC_DEFUN([ABI_FC_EXTENSIONS],[
  dnl Set Fortran module extension
  AX_F90_MODULE_EXTENSION
  if test "${ax_cv_f90_modext}" != ""; then
    MODEXT="${ax_cv_f90_modext}"
  else
    MODEXT="mod"
    AC_MSG_NOTICE([setting Fortran module extension to ".${MODEXT}"])
  fi
  AC_SUBST(MODEXT)

  dnl Change the default Fortran extension for tests
  AC_FC_SRCEXT(F90,[abi_fc_src_ok="yes"],[abi_fc_src_ok="no"])
  if test "${abi_fc_src_ok}" != "yes"; then
    AC_MSG_WARN([Fortran file extension could not be changed])
    AC_MSG_WARN([some advanced Fortran tests may fail])
  fi
]) # ABI_FC_EXTENSIONS



# ABI_FC_FEATURES()
# -----------------
#
# Explores the capabilities of the Fortran compiler.
#
AC_DEFUN([ABI_FC_FEATURES],[
  dnl Explore compiler peculiarities
  _ABI_CHECK_FC_DTARRAYS
  _ABI_CHECK_FC_ISO_C_BINDING
  _ABI_CHECK_FC_EXIT
  _ABI_CHECK_FC_FLUSH
  _ABI_CHECK_FC_FLUSH_
  _ABI_CHECK_FC_GAMMA
  _ABI_CHECK_FC_GETENV
  _ABI_CHECK_FC_GETPID
  _ABI_CHECK_FC_NULL
  _ABI_CHECK_FC_INT_QUAD
  _ABI_CHECK_FC_LONG_LINES
  _ABI_CHECK_FC_STREAM_IO
  _ABI_CHECK_FC_CPUTIME
  _ABI_CHECK_FC_TIMING
]) # ABI_FC_FEATURES



# ABI_FC_MOD_CASE()
# -----------------
#
# Checks whether the Fortran compiler creates upper-case or lower-case
# module files.
#
AC_DEFUN([ABI_FC_MOD_CASE],[
  AC_REQUIRE([ABI_FC_EXTENSIONS])

  dnl Init
  fc_mod_lowercase="yes"
  fc_mod_uppercase="no"
  AC_MSG_NOTICE([determining Fortran module case])

  dnl Compile a dummy module
  AC_LANG_PUSH([Fortran])
  AC_COMPILE_IFELSE([[
    module conftest
    end module conftest
  ]],[],[AC_MSG_FAILURE([unable to compile a simple Fortran module])])
  AC_LANG_POP([Fortran])

  dnl Check module file existence
  if test -f "CONFTEST.${MODEXT}"; then
    fc_mod_lowercase="no"
    fc_mod_uppercase="yes"
  elif test ! -f "conftest.${MODEXT}"; then
    AC_MSG_WARN([conftest.${MODEXT} Fortran module could not be found])
  fi

  dnl Output final outcome
  AC_MSG_CHECKING([whether Fortran modules are upper-case])
  AC_MSG_RESULT([${fc_mod_uppercase}])
]) # ABI_FC_MOD_CASE



# ABI_FC_MOD_INCS(MODULE)
# -----------------------
#
# Checks whether the specified Fortran module is directly available, or
# if we need to add '-I/usr/include' to the compile flags. Returns the
# required includes.
#
AC_DEFUN([ABI_FC_MOD_INCS],[
  AC_MSG_CHECKING([for Fortran module includes])

  if test "${abi_fc_mod_incs_ok}" = "" -o \
          "${abi_fc_mod_incs_ok}" = "unknown"; then

    dnl Init
    fc_mod_incs=""

    dnl Prepare environment
    tmp_saved_FCFLAGS="${FCFLAGS}"
    AC_LANG_PUSH([Fortran])

    dnl Look for module without includes
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
      [[
        use $1
      ]])], [abi_fc_mod_incs_ok="none required"], [abi_fc_mod_incs_ok="unknown"])

    dnl Look for module with includes
    if test "${abi_fc_mod_incs_ok}" = "unknown"; then
      FCFLAGS="${FCFLAGS} -I/usr/include"
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
        [[
          use $1
        ]])],
        [abi_fc_mod_incs_ok="-I/usr/include"; fc_mod_incs="-I/usr/include"],
        [abi_fc_mod_incs_ok="unknown"])
    fi
    AC_MSG_RESULT([${abi_fc_mod_incs_ok}])

    dnl Restore environment
    AC_LANG_POP([Fortran])
    FCFLAGS="${tmp_saved_FCFLAGS}"

  else

    AC_MSG_RESULT([${abi_fc_mod_incs_ok} (cached)])

  fi

  dnl Substitute variables
  AC_SUBST(fc_mod_incs)
]) # ABI_FC_MOD_INCS



# ABI_PROG_FC()
# -------------
#
# Tries to determine which type of Fortran compiler is installed.
#
AC_DEFUN([ABI_PROG_FC],[
  dnl Init
  abi_fc_vendor="${with_fc_vendor}"
  abi_fc_version="${with_fc_version}"
  tmp_fc_info_file="${abinit_builddir}/config.fc_info.tmp"

  if test "${abi_fc_vendor}" = ""; then
    abi_fc_vendor="unknown"
  fi
  if test "${abi_fc_version}" = ""; then
    abi_fc_version="unknown"
  fi
  abi_fc_wrap="no"

  dnl Determine Fortran compiler type (the order is important)
  AC_MSG_CHECKING([which type of Fortran compiler we have])

  dnl Clear temporary info file
  rm -f "${tmp_fc_info_file}"

  dnl Get rid of that one as early as possible
  if test "${abi_fc_vendor}" = "unknown"; then
    _ABI_CHECK_FC_IBM(${FC})
  fi
  echo "${fc_info_string}" >>"${tmp_fc_info_file}"

  dnl Should be checked before gfortran because it mimics its behaviour
  if test "${abi_fc_vendor}" = "unknown"; then
    _ABI_CHECK_FC_INTEL(${FC})
  fi
  echo "${fc_info_string}" >>"${tmp_fc_info_file}"

  if test "${abi_fc_vendor}" = "unknown"; then
    _ABI_CHECK_FC_G95(${FC})
  fi
  echo "${fc_info_string}" >>"${tmp_fc_info_file}"

  if test "${abi_fc_vendor}" = "unknown"; then
    _ABI_CHECK_FC_GNU(${FC})
  fi
  echo "${fc_info_string}" >>"${tmp_fc_info_file}"

  if test "${abi_fc_vendor}" = "unknown"; then
    _ABI_CHECK_FC_PATHSCALE(${FC})
  fi
  echo "${fc_info_string}" >>"${tmp_fc_info_file}"

  #if test "${abi_fc_vendor}" = "unknown"; then
  #  _ABI_CHECK_FC_COMPAQ(${FC})
  #fi
  #echo "${fc_info_string}" >>"${tmp_fc_info_file}"

  if test "${abi_fc_vendor}" = "unknown"; then
    _ABI_CHECK_FC_ABSOFT(${FC})
  fi
  echo "${fc_info_string}" >>"${tmp_fc_info_file}"

  if test "${abi_fc_vendor}" = "unknown"; then
    _ABI_CHECK_FC_MIPSPRO(${FC})
  fi
  echo "${fc_info_string}" >>"${tmp_fc_info_file}"

  if test "${abi_fc_vendor}" = "unknown"; then
    _ABI_CHECK_FC_OPEN64(${FC})
  fi
  echo "${fc_info_string}" >>"${tmp_fc_info_file}"

  if test "${abi_fc_vendor}" = "unknown"; then
    _ABI_CHECK_FC_FUJITSU(${FC})
  fi
  echo "${fc_info_string}" >>"${tmp_fc_info_file}"

  if test "${abi_fc_vendor}" = "unknown"; then
    _ABI_CHECK_FC_SUN(${FC})
  fi
  echo "${fc_info_string}" >>"${tmp_fc_info_file}"

  if test "${abi_fc_vendor}" = "unknown"; then
    _ABI_CHECK_FC_HITACHI(${FC})
  fi
  echo "${fc_info_string}" >>"${tmp_fc_info_file}"

  if test "${abi_fc_vendor}" = "unknown"; then
    _ABI_CHECK_FC_NAG(${FC})
  fi
  echo "${fc_info_string}" >>"${tmp_fc_info_file}"

  if test "${abi_fc_vendor}" = "unknown"; then
    _ABI_CHECK_FC_PGI(${FC})
  fi
  echo "${fc_info_string}" >>"${tmp_fc_info_file}"

  dnl Fall back to generic when detection fails
  if test "${abi_fc_vendor}" = "unknown"; then
    abi_fc_vendor="generic"
  #else
  #  rm -f "${tmp_fc_info_file}"
  fi

  dnl Normalize Fortran compiler version
  if test "${abi_fc_version}" = "unknown"; then
    abi_fc_version="0.0"
  else
    abi_fc_version=`echo ${abi_fc_version} | cut -d. -f1-2`
  fi

  dnl Display final result
  AC_MSG_RESULT([${abi_fc_vendor} ${abi_fc_version}])

  dnl Schedule compiler info for substitution
  AC_SUBST(abi_fc_vendor)
  AC_SUBST(abi_fc_version)
  AC_SUBST(abi_fc_wrap)
  AC_SUBST(fc_info_string)
]) # ABI_PROG_FC
