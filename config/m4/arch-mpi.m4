# -*- Autoconf -*-
#
# Copyright (C) 2005-2012 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# MPI support for ABINIT
#



# _AX_MPIRUN_IFELSE(PROGRAM, [ACTION-IF-TRUE], [ACTION-IF-FALSE])
# ---------------------------------------------------------------
#
# Compile, link, and run in parallel.
# This macro can be used during the selection of a compiler.
# We also remove conftest.o as if the compilation fails, some compilers
# don't remove it.  We remove gmon.out and bb.out, which may be
# created during the run if the program is built with profiling support.
#
m4_define([_AX_MPIRUN_IFELSE],[
  m4_ifvaln([$1], [AC_LANG_CONFTEST([$1])])dnl
  rm -f conftest$ac_exeext
  test "${MPI_2PROCS}" = "" && MPI_2PROCS="-np 2"
  AS_IF([_AC_DO_VAR(ac_link) && _AC_DO_TOKENS(${MPI_RUNNER} ${MPI_2PROCS} conftest$ac_exeext) >&AS_MESSAGE_LOG_FD],
    [$2],
    [AS_ECHO(["$as_me: program exited with status $ac_status"]) >&AS_MESSAGE_LOG_FD
  _AC_MSG_LOG_CONFTEST
  m4_ifvaln([$3],
    [( exit $ac_status )
  $3])dnl])[]dnl
  rm -rf conftest.dSYM
  rm -f core *.core core.conftest.* gmon.out bb.out conftest$ac_exeext conftest.$ac_objext m4_ifval([$1],
    [conftest.$ac_ext])[]dnl
])# _AX_MPIRUN_IFELSE



# AX_MPIRUN_IFELSE(PROGRAM,
#                  [ACTION-IF-TRUE], [ACTION-IF-FALSE],
#                  [ACTION-IF-CROSS-COMPILING = RUNTIME-ERROR])
# -------------------------------------------------------------
#
# Compile, link, and run in parallel. Requires that the compiler for the
# current language was checked for, hence do not use this macro in macros
# looking for a compiler.
#
AC_DEFUN([AX_MPIRUN_IFELSE],
[AC_LANG_COMPILER_REQUIRE()dnl
m4_ifval([$4], [],
  [AC_DIAGNOSE([cross],
    [$0 called without default to allow cross compiling])])dnl
  AS_IF([test "$cross_compiling" = yes],
    [m4_default([$4],
      [AC_MSG_FAILURE([cannot run test program while cross compiling])])],
      [_AX_MPIRUN_IFELSE($@)])
]) # AX_MPIRUN_IFELSE



                    ########################################



# _ABI_MPI_CHECK_CC()
# -------------------
#
# Checks whether the C compiler is able to produce MPI binaries.
#
AC_DEFUN([_ABI_MPI_CHECK_CC],[
  dnl Set default values
  abi_mpi_cc_ok="no"

  dnl Try to compile a C MPI program
  AC_MSG_CHECKING([whether the C compiler supports MPI])

  dnl Back-up build environment
  ABI_ENV_BACKUP

  dnl Prepare build environment
  CPPFLAGS="${CPPFLAGS} ${lib_mpi_incs}"
  LDFLAGS="${CC_LDFLAGS}"
  LIBS="${CC_LIBS} ${lib_mpi_libs}"

  AC_LANG_PUSH([C])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [[#include <stdlib.h>
#include "mpi.h"]],
    [[
      int rc;

      MPI_Init(NULL,NULL);
      rc = MPI_Finalize();
    ]])], [abi_mpi_cc_ok="yes"], [abi_mpi_cc_ok="no"])
  AC_LANG_POP([C])

  dnl Restore build environment
  ABI_ENV_RESTORE

  AC_MSG_RESULT([${abi_mpi_cc_ok}])
]) # _ABI_MPI_CHECK_CC



                    ########################################



# _ABI_MPI_CHECK_CXX()
# --------------------
#
# Checks whether the C++ compiler is able to produce MPI binaries.
#
AC_DEFUN([_ABI_MPI_CHECK_CXX],[
  dnl Set default values
  abi_mpi_cxx_ok="no"

  dnl Try to compile a C++ MPI program
  AC_MSG_CHECKING([whether the C++ compiler supports MPI])

  dnl Back-up build environment
  ABI_ENV_BACKUP

  dnl Prepare build environment
  CPPFLAGS="${CPPFLAGS} ${lib_mpi_incs}"
  LDFLAGS="${CXX_LDFLAGS}"
  LIBS="${CXX_LIBS} ${lib_mpi_libs}"

  AC_LANG_PUSH([C++])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [[@%:@include "mpi.h"]],
    [[
      MPI::Init();
      MPI::Finalize();
    ]])], [abi_mpi_cxx_ok="yes"], [abi_mpi_cxx_ok="no"])
  AC_LANG_POP([C++])

  dnl Restore build environment
  ABI_ENV_RESTORE

  AC_MSG_RESULT([${abi_mpi_cxx_ok}])
]) # _ABI_MPI_CHECK_CXX



                    ########################################



# _ABI_MPI_CHECK_FC()
# -------------------
#
# Checks whether the Fortran compiler is able to produce MPI binaries.
#
AC_DEFUN([_ABI_MPI_CHECK_FC],[
  dnl Set default values
  abi_mpi_cxx_ok="no"

  dnl Try to compile a Fortran MPI program
  AC_MSG_CHECKING([whether the Fortran Compiler supports MPI])

  dnl Back-up build environment
  ABI_ENV_BACKUP

  dnl Prepare build environment
  CPPFLAGS="${CPPFLAGS} ${lib_mpi_incs}"
  LDFLAGS="${FC_LDFLAGS}"
  LIBS="${FC_LIBS} ${lib_mpi_libs}"

  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[  
      include "mpif.h"
      integer :: ierr
      call mpi_init(ierr)
      call mpi_finalize(ierr)
    ]])], [abi_mpi_fc_ok="yes"], [abi_mpi_fc_ok="no"])
  AC_LANG_POP([Fortran])

  dnl Restore build environment
  ABI_ENV_RESTORE

  AC_MSG_RESULT([${abi_mpi_fc_ok}])
]) # _ABI_MPI_CHECK_FC



                    ########################################



# _ABI_MPI_CHECK_FC_LEVEL()
# -------------------------
#
# Checks which MPI level is supported by the Fortran compiler.
#
AC_DEFUN([_ABI_MPI_CHECK_FC_LEVEL],[
  dnl Set default values
  abi_mpi_fc_level="none"

  dnl Try to compile a MPI-2 Fortran program
  AC_MSG_CHECKING([which level of MPI is supported by the Fortran compiler])

  if test "${abi_mpi_fc_ok}" = "yes"; then

    dnl Back-up build environment
    ABI_ENV_BACKUP

    dnl Prepare build environment
    CPPFLAGS="${CPPFLAGS} ${lib_mpi_incs}"
    LDFLAGS="${FC_LDFLAGS}"
    LIBS="${FC_LIBS} ${lib_mpi_libs}"

    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[  
              use mpi
              integer :: ierr
              call mpi_init(ierr)
              call mpi_finalize(ierr)
      ]])], [abi_mpi_fc_level="2"], [abi_mpi_fc_level="1"])
    AC_LANG_POP([Fortran])

    dnl Restore build environment
    ABI_ENV_RESTORE
  fi

  AC_MSG_RESULT([${abi_mpi_fc_level}])
]) # _ABI_MPI_CHECK_FC_LEVEL



                    ########################################



# _ABI_MPI_CHECK_CREATE_TYPE_STRUCT()
# -----------------------------------
#
# Checks whether the MPI library supports MPI_CREATE_TYPE_STRUCT (MPI2)
#
AC_DEFUN([_ABI_MPI_CHECK_CREATE_TYPE_STRUCT],[
  dnl Set default values
  abi_mpi_create_type_struct_ok="no"

  dnl Try to compile a Fortran MPI program
  AC_MSG_CHECKING([whether the MPI library supports MPI_CREATE_TYPE_STRUCT])

  if test "${abi_mpi_fc_level}" = "2"; then

    dnl No problem should appear for MPI2 but we test it anyway.

    dnl Back-up build environment
    ABI_ENV_BACKUP
                                                                                              
    dnl Prepare build environment
    CPPFLAGS="${CPPFLAGS} ${lib_mpi_incs}"
    LDFLAGS="${FC_LDFLAGS}"
    LIBS="${FC_LIBS} ${lib_mpi_libs}"
                                                                                              
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[  
        use mpi 
        integer,parameter :: ncount=10
        integer :: ierr,new_type
        integer :: block_length(ncount),block_type(ncount)
        integer(MPI_ADDRESS_KIND) :: block_displ(ncount)
                                                                                              
        call mpi_init(ierr)
        call MPI_TYPE_CREATE_STRUCT(ncount,block_length,block_displ,block_type,new_type,ierr)
        call mpi_finalize(ierr)
                                                                                              
      ]])], [abi_mpi_create_type_struct_ok="yes"], [abi_mpi_create_type_struct_ok="no"])
    AC_LANG_POP
                                                                                              
    dnl Restore build environment
    ABI_ENV_RESTORE

  else

    dnl Back-up build environment
    ABI_ENV_BACKUP

    dnl Prepare build environment
    CPPFLAGS="${CPPFLAGS} ${lib_mpi_incs}"
    LDFLAGS="${FC_LDFLAGS}"
    LIBS="${FC_LIBS} ${lib_mpi_libs}"

    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[  
        include "mpif.h"
        integer,parameter :: ncount=10
        integer :: ierr,new_type
        integer :: block_length(ncount),block_type(ncount)
        integer(MPI_ADDRESS_KIND) :: block_displ(ncount)

        call mpi_init(ierr)
        call MPI_TYPE_CREATE_STRUCT(ncount,block_length,block_displ,block_type,new_type,ierr)
        call mpi_finalize(ierr)

      ]])], [abi_mpi_create_type_struct_ok="yes"], [abi_mpi_create_type_struct_ok="no"])
    AC_LANG_POP

    dnl Restore build environment
    ABI_ENV_RESTORE

  fi

  AC_MSG_RESULT([${abi_mpi_create_type_struct_ok}])

  if test "${abi_mpi_create_type_struct_ok}" = "yes"; then
    AC_DEFINE([HAVE_MPI_TYPE_CREATE_STRUCT],1,
      [Define to 1 if your MPI library supports MPI_TYPE_CREATE_STRUCT.])
  fi

]) # _ABI_MPI_CHECK_CREATE_TYPE_STRUCT     



                    ########################################



# _ABI_MPI_CREATE_WRAPPER(COMPILER_TYPE, SERIAL_COMPILER, MPI_COMPILER)
# ---------------------------------------------------------------------
#
# Creates a wrapper for MPI compilers when they can be configured to
# accept different sequential compilers.
#
# Note: this is impossible with the Autotools, because they require CC
#       CXX, and FC to be set to the actual compilers.
#
AC_DEFUN([_ABI_MPI_CREATE_WRAPPER],[
  dnl Init
  tmp_name=`echo "$1" | sed -e 's/.*/\L&/'`
  ${MKDIR_P} config/wrappers

  dnl Create file
  cat >${abinit_builddir}/config/wrappers/wrap-mpi${tmp_name} <<EOF
#!/bin/sh

$1="$2"
export $1

$3 \[$]{*}
EOF

  dnl Fix permissions
  chmod u+x ${abinit_builddir}/config/wrappers/wrap-mpi${tmp_name}

  dnl Overwrite compiler setting
  eval $1="${abinit_builddir}/config/wrappers/wrap-mpi${tmp_name}"

  dnl Clean-up the mess
  unset tmp_name
]) # _ABI_MPI_CREATE_WRAPPER



                    ########################################



# ABI_MPI_INIT()
# --------------
#
# Looks for an implementation of MPI, using the provided prefix.
# Note 1: this is a convenience feature, purely for comfort.
# Note 2: should be run as early as possible
#
AC_DEFUN([ABI_MPI_INIT],[
  dnl Init
  abi_mpi_complete="unknown"
  abi_mpi_has_cc="no"
  abi_mpi_has_cxx="no"
  abi_mpi_has_fc="no"
  abi_mpi_has_runner="no"
  abi_mpi_level="${with_mpi_level}"
  abi_mpi_usable="no"
  lib_mpi_incs="${with_mpi_incs}"
  lib_mpi_libs="${with_mpi_libs}"
  test "${MPI_RUNNER}" = "" && MPI_RUNNER="${with_mpi_runner}"
  MPI_CC=""
  MPI_CXX=""
  MPI_FC=""

  if test "${enable_mpi}" = "yes"; then

    dnl Banner
    AC_MSG_NOTICE([Initializing MPI support])

    dnl Check whether to look for generic files
    if test "${with_mpi_prefix}" != ""; then
      AC_MSG_NOTICE([looking for MPI in ${with_mpi_prefix}])

      dnl Look for incompatibilities
      if test "${CPP}" != ""; then
        AC_MSG_WARN([${CPP} might not be fully compatible with MPI])
      fi
      if test "${MPI_RUNNER}" != ""; then
        AC_MSG_ERROR([use --with-mpi-prefix or set MPI_RUNNER, not both])
      fi
      if test "${with_mpi_incs}" != ""; then
        AC_MSG_ERROR([use --with-mpi-prefix or --with-mpi-includes, not both])
      fi
      if test "${with_mpi_level}" != ""; then
        AC_MSG_WARN([forcing MPI level to ${with_mpi_level} might make the build fail])
      fi
      if test "${with_mpi_libs}" != ""; then
        AC_MSG_ERROR([use --with-mpi-prefix or --with-mpi-libs, not both])
      fi
      if test "${with_mpi_runner}" != ""; then
        AC_MSG_ERROR([use --with-mpi-prefix or --with-mpi-runner, not both])
      fi

      dnl Look for a C compiler
      AC_MSG_CHECKING([for a MPI C compiler])
      if test -x "${with_mpi_prefix}/bin/mpicc"; then
        abi_mpi_has_cc="yes"
        MPI_CC="${with_mpi_prefix}/bin/mpicc"
      fi
      if test "${MPI_CC}" = ""; then
        AC_MSG_RESULT([none found])
      else
        AC_MSG_RESULT([${MPI_CC}])
      fi

      dnl Look for a C++ compiler
      AC_MSG_CHECKING([for a MPI C++ compiler])
      if test -x "${with_mpi_prefix}/bin/mpicxx"; then
        abi_mpi_has_cxx="yes"
        MPI_CXX="${with_mpi_prefix}/bin/mpicxx"
      elif test -x "${with_mpi_prefix}/bin/mpic++"; then
        abi_mpi_has_cxx="yes"
        MPI_CXX="${with_mpi_prefix}/bin/mpic++"
      fi
      if test "${MPI_CXX}" = ""; then
        AC_MSG_RESULT([none found])
      else
        AC_MSG_RESULT([${MPI_CXX}])
      fi

      dnl Look for a Fortran 90 compiler
      AC_MSG_CHECKING([for a MPI Fortran compiler])
      if test -x "${with_mpi_prefix}/bin/mpif90"; then
        abi_mpi_has_fc="yes"
        MPI_FC="${with_mpi_prefix}/bin/mpif90"
      fi
      if test "${MPI_FC}" = ""; then
        AC_MSG_RESULT([none found])
      else
        AC_MSG_RESULT([${MPI_FC}])
      fi

      dnl Look for a runner
      AC_MSG_CHECKING([for a MPI runner])
      if test -x "${with_mpi_prefix}/bin/mpirun"; then
        abi_mpi_has_runner="yes"
        MPI_RUNNER="${with_mpi_prefix}/bin/mpirun"
      elif test -x "${with_mpi_prefix}/bin/mpiexec"; then
        abi_mpi_has_runner="yes"
        MPI_RUNNER="${with_mpi_prefix}/bin/mpiexec"
      elif test -x "${with_mpi_prefix}/bin/dmpirun"; then
        abi_mpi_has_runner="yes"
        MPI_RUNNER="${with_mpi_prefix}/bin/dmpirun"
      fi
      if test "${MPI_RUNNER}" = ""; then
        AC_MSG_RESULT([none found])
      else
        AC_MSG_RESULT([${MPI_RUNNER}])
      fi

      dnl Report whether generic MPI implementation is sufficiently complete
      if test "${abi_mpi_has_cc}" = "yes" -a \
              "${abi_mpi_has_fc}" = "yes" -a \
              "${abi_mpi_has_runner}" = "yes"; then
        abi_mpi_complete="yes"

        dnl "Baise-couillon" for those who set compilers twice
        dnl FIXME: does not handle the case when full path is not
        dnl        provided, e.g. prefix=/usr, CC=mpicc.
        dnl Hint: TRUE_CC=`echo "${CC}" | cut -d' ' -f1`
        dnl       TRUE_CC=`type -t "${TRUE_CC}"`
        dnl       if "file" type -p ...
        tmp_chk_cc="no"
        tmp_chk_redundant="no"
        if test "${CC}" = "${MPI_CC}"; then
          tmp_chk_cc="yes"
          tmp_chk_redundant="yes"
          AC_MSG_WARN([redundant setting of MPI C compiler!
    Use --with-mpi-prefix preferably.])
        fi
        tmp_chk_cxx="no"
        if test "${CXX}" = "${MPI_CXX}"; then
          tmp_chk_cxx="yes"
          AC_MSG_WARN([redundant setting of MPI C++ compiler!
    Use --with-mpi-prefix preferably.])
        fi
        if test "${tmp_chk_cxx}" != "${tmp_chk_redundant}"; then
          AC_MSG_WARN([inconsistent compiler settings!
    Use --with-mpi-prefix or set (CC, CXX, FC), not both.])
        fi
        tmp_chk_fc="no"
        if test "${FC}" = "${MPI_FC}"; then
          tmp_chk_fc="yes"
          AC_MSG_WARN([redundant setting of MPI Fortran compiler
    Use --with-mpi-prefix preferably.])
        fi
        if test "${tmp_chk_fc}" != "${tmp_chk_redundant}"; then
          AC_MSG_ERROR([inconsistent compiler settings!
    Use --with-mpi-prefix or set (CC, CXX, FC), not both.])
        fi
        if test "${tmp_chk_redundant}" = "yes"; then
          CC=""
          CXX=""
          FC=""
          AC_MSG_NOTICE([ignoring CC, CXX, and FC settings])
        fi
        unset tmp_chk_cc tmp_chk_cxx tmp_chk_fc tmp_chk_redundant

        dnl Decide whether to wrap MPI compiler calls
        if test "${CC}" = ""; then
          CC="${MPI_CC}"
        else
          AC_MSG_NOTICE([creating wrapper for MPI C compiler])
          _ABI_MPI_CREATE_WRAPPER([CC],[${CC}],[${MPI_CC}])
        fi
        if test "${CXX}" = ""; then
          CXX="${MPI_CXX}"
        else
          AC_MSG_NOTICE([creating wrapper for MPI C++ compiler])
          _ABI_MPI_CREATE_WRAPPER([CXX],[${CXX}],[${MPI_CXX}])
        fi
        if test "${FC}" = ""; then
          FC="${MPI_FC}"
        else
          AC_MSG_NOTICE([creating wrapper for MPI Fortran compiler])
          _ABI_MPI_CREATE_WRAPPER([FC],[${FC}],[${MPI_FC}])
        fi
      else
        unset MPI_CC
        unset MPI_CXX
        unset MPI_FC
        unset MPI_RUNNER
        abi_mpi_complete="no"
      fi

    else

      dnl Look for a MPI runner if necessary
      if test "${MPI_RUNNER}" = ""; then
        AC_CHECK_PROGS(MPI_RUNNER,[mpirun mpiexec dmpirun])
        if test "${MPI_RUNNER}" != ""; then
          AC_MSG_WARN([MPI runner ${MPI_RUNNER} may be incompatible with MPI compilers])
        fi
      else
        abi_mpi_runner_works="no"
        if test -x "${MPI_RUNNER}"; then
          abi_mpi_runner_works="yes"
        else
          unset abi_mpi_runner_works
          AC_CHECK_PROG(abi_mpi_runner_works,${MPI_RUNNER},yes,no)
        fi
        if test "${abi_mpi_runner_works}" = "no"; then
          AC_MSG_ERROR([invalid MPI runner: ${MPI_RUNNER}])
        fi
      fi

      dnl Inform about compiler checks
      AC_MSG_NOTICE([compiler checks deferred])
  
      dnl Report whether MPI implementation is sufficiently complete
      if test "${MPI_RUNNER}" != ""; then
        abi_mpi_complete="yes"
      else
        abi_mpi_complete="no"
      fi

    fi dnl with_mpi_prefix

  else

    AC_MSG_NOTICE([MPI support disabled from command-line])
    enable_mpi_io="no"
    enable_mpi_trace="no"
    with_mpi_level=""
    with_mpi_prefix=""

  fi dnl enable_mpi

  dnl Enable substitution
  AC_SUBST(MPI_RUNNER)
  AC_SUBST(lib_mpi_fcflags)
  AC_SUBST(lib_mpi_ldflags)
  AC_SUBST(lib_mpi_incs)
  AC_SUBST(lib_mpi_libs)
]) # ABI_MPI_INIT



                    ########################################



# ABI_MPI_DETECT()
# ----------------
#
# Tries first to determine whether the MPI implementation is usable,
# then takes appropriate actions.
#
AC_DEFUN([ABI_MPI_DETECT],[
  dnl Init
  AC_REQUIRE([ABI_MPI_INIT])
  lib_mpi_fcflags=""
  lib_mpi_ldflags=""
  lib_mpi_incs=""
  lib_mpi_libs=""

  dnl Check whether MPI is usable
  if test "${abi_mpi_complete}" = "yes"; then
    _ABI_MPI_CHECK_CC
    _ABI_MPI_CHECK_CXX
    _ABI_MPI_CHECK_FC

    if test "${abi_mpi_cc_ok}" = "yes" -a \
                    "${abi_mpi_fc_ok}" = "yes"; then
      abi_mpi_usable="yes"
    fi
  fi
  AC_MSG_CHECKING([whether MPI is usable])
  AC_MSG_RESULT([${abi_mpi_usable}])

  dnl Make sure that main trigger is set
  if test "${enable_mpi}" = ""; then
    if test "${abi_mpi_usable}" = "yes"; then
      AC_MSG_NOTICE([enabling MPI support])
      enable_mpi="yes"
    else
      if test "${abi_mpi_complete}" = "no"; then
        AC_MSG_NOTICE([disabling MPI support])
      fi
      enable_mpi="no"
      enable_mpi_io="no"
      enable_mpi_trace="no"
    fi
  else
    if test "${enable_mpi}" = "yes" -a "${abi_mpi_usable}" = "no"; then
      AC_MSG_WARN([MPI support is broken!])
    fi
  fi

  dnl Set the I/O trigger accordingly
  if test "${enable_mpi_io}" = ""; then
    if test "${enable_mpi}" = "yes"; then
      AC_MSG_NOTICE([enabling MPI I/O support])
      enable_mpi_io="yes"
    else
      AC_MSG_NOTICE([disabling MPI I/O support])
      enable_mpi_io="no"
    fi
  else
    if test "${enable_mpi}" = "no" -a "${enable_mpi_io}" = "yes"; then
      AC_MSG_WARN([disabling MPI I/O support since MPI is disabled])
      enable_mpi_io="no"
    fi
    if test "${enable_mpi}" = "yes" -a "${enable_mpi_io}" = "no"; then
      AC_MSG_WARN([disabling MPI I/O is not recommended])
    fi
  fi

  dnl Report status
  AC_MSG_CHECKING([whether to build MPI code])
  AC_MSG_RESULT([${enable_mpi}])

  dnl Continue with advanced MPI checks
  if test "${enable_mpi}" = "yes"; then

    dnl Propagate main trigger
    AC_DEFINE([HAVE_MPI],1,[Define to 1 if you want to enable MPI support.])

    dnl Propagate MPI I/O trigger
    AC_MSG_CHECKING([whether to build MPI I/O code])
    AC_MSG_RESULT([${enable_mpi_io}])
    if test "${enable_mpi_io}" = "yes"; then
      AC_DEFINE([HAVE_MPI_IO],1,[Define to 1 if you want MPI I/O support.])
    fi

    dnl Set MPI time tracing support
    AC_MSG_CHECKING([whether to build MPI time tracing code])
    AC_MSG_RESULT([${enable_mpi_trace}])
    if test "${enable_mpi_trace}" = "yes"; then
      AC_DEFINE([HAVE_MPI_TRACE],1,[Define to 1 if you want MPI time tracing support.])
    fi

    dnl Check MPI level actually supported
    _ABI_MPI_CHECK_FC_LEVEL

    dnl Select MPI level
    if test "${abi_mpi_level}" = ""; then
      abi_mpi_level="${abi_mpi_fc_level}"
    else
      AC_MSG_NOTICE([forcing MPI-${abi_mpi_level} standard support])
      if test "${abi_mpi_level}" != "${abi_mpi_fc_level}"; then
        AC_MSG_WARN([detected MPI-${abi_mpi_fc_level} support but using MPI-${abi_mpi_level} instructions])
      fi
    fi

    dnl Propagate MPI level
    case "${abi_mpi_level}" in
      1)
        AC_DEFINE([HAVE_MPI1],1,[Define to 1 if you have a MPI-1 implementation.])
        ;;
      2)
        AC_DEFINE([HAVE_MPI2],1,[Define to 1 if you have a MPI-2 implementation.])
        ;;
    esac

  else

    lib_mpi_incs=""
    lib_mpi_libs=""
    MPI_RUNNER=""

  fi

  dnl Test the availability of problematic MPI primitives
  _ABI_MPI_CHECK_CREATE_TYPE_STRUCT()

  AM_CONDITIONAL(DO_TEST_MPI,[test "${enable_mpi}" = "yes"])
]) # ABI_MPI_DETECT
