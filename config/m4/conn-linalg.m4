# -*- Autoconf -*-
#
# Copyright (C) 2005-2012 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for external linear algebra libraries
#



# _ABI_LINALG_CHECK_LIBS()
# ------------------------
#
# Check whether the specified libraries are BLAS and LAPACK
# implementations.
#
AC_DEFUN([_ABI_LINALG_CHECK_LIBS],[
  dnl Init
  abi_linalg_has_blas="no"
  abi_linalg_has_lapack="no"
  abi_linalg_has_blacs="no"
  abi_linalg_has_scalapack="no"
  abi_linalg_has_magma="no"
  abi_linalg_has_plasma="no"

  dnl Prepare environment
  tmp_saved_LIBS="${LIBS}"
  LIBS="${LIBS} ${lib_gpu_libs} ${lib_mpi_libs}"

  dnl BLAS?
  AC_MSG_CHECKING([for BLAS support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [
      call zgemm
    ])], [abi_linalg_has_blas="yes"], [abi_linalg_has_blas="no"])
  AC_MSG_RESULT([${abi_linalg_has_blas}])

  dnl LAPACK?
  AC_MSG_CHECKING([for LAPACK support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [
      call zhpev
    ])], [abi_linalg_has_lapack="yes"], [abi_linalg_has_lapack="no"])
  AC_MSG_RESULT([${abi_linalg_has_lapack}])

  dnl BLACS?
  AC_MSG_CHECKING([for BLACS support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [
      call blacs_gridinit
    ])], [abi_linalg_has_blacs="yes"], [abi_linalg_has_blacs="no"])
  AC_MSG_RESULT([${abi_linalg_has_blacs}])

  dnl ScaLAPACK?
  AC_MSG_CHECKING([for ScaLAPACK support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [
      call pzheevx
    ])], [abi_linalg_has_scalapack="yes"], [abi_linalg_has_scalapack="no"])
  AC_MSG_RESULT([${abi_linalg_has_scalapack}])

  dnl MAGMA?
  AC_MSG_CHECKING([for MAGMA support in specified libraries])
  AC_MSG_CHECKING([TEST])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [
      call magmaf_zheevd
    ])], [abi_linalg_has_magma="yes"], [abi_linalg_has_magma="no"])
  AC_MSG_RESULT([${abi_linalg_has_magma}])

  dnl PLASMA?
  dnl FIXME: implement test

  dnl Restore environment
  LIBS="${tmp_saved_LIBS}"
]) # _ABI_LINALG_CHECK_LIBS



# _ABI_LINALG_SEARCH_BLAS(BLAS,EXTRA_LIBS)
# ----------------------------------------
#
# Look for a BLAS implementation.
#
AC_DEFUN([_ABI_LINALG_SEARCH_BLAS],[
  dnl Init
  abi_linalg_has_blas="no"

  dnl Look for libraries and routines
  AC_SEARCH_LIBS([zgemm],$1,
    [abi_linalg_has_blas="yes"],[abi_linalg_has_blas="no"],
    [$2 ${abi_linalg_libs}])
  if test "${abi_linalg_has_blas}" = "yes"; then
    if test "${ac_cv_search_zgemm}" != "none required"; then
      abi_linalg_libs="${ac_cv_search_zgemm} $2 ${abi_linalg_libs}"
    fi
  fi
]) # _ABI_LINALG_SEARCH_BLAS



# _ABI_LINALG_SEARCH_LAPACK(LAPACK,EXTRA_LIBS)
# --------------------------------------------
#
# Look for a LAPACK implementation.
#
AC_DEFUN([_ABI_LINALG_SEARCH_LAPACK],[
  dnl Init
  abi_linalg_has_lapack="no"

  dnl Look for libraries and routines
  AC_SEARCH_LIBS([zhpev],$1,
    [abi_linalg_has_lapack="yes"],[abi_linalg_has_lapack="no"],
    [$2 ${abi_linalg_libs}])
  if test "${abi_linalg_has_lapack}" = "yes"; then
    if test "${ac_cv_search_zhpev}" != "none required"; then
      abi_linalg_libs="${ac_cv_search_zhpev} $2 ${abi_linalg_libs}"
    fi
  fi
]) # _ABI_LINALG_SEARCH_LAPACK



# _ABI_LINALG_SEARCH_BLACS(BLACS,EXTRA_LIBS)
# --------------------------------------------
#
# Look for a BLACS implementation.
#
AC_DEFUN([_ABI_LINALG_SEARCH_BLACS],[
  dnl Init
  abi_linalg_has_blacs="no"

  dnl Look for libraries and routines
  AC_SEARCH_LIBS([blacs_gridinit],$1,
    [abi_linalg_has_blacs="yes"],[abi_linalg_has_blacs="no"],
    [$2 ${abi_linalg_libs}])
  if test "${abi_linalg_has_blacs}" = "yes"; then
    if test "${ac_cv_search_blacs_gridinit}" != "none required"; then
      abi_linalg_libs="${ac_cv_search_blacs_gridinit} $2 ${abi_linalg_libs}"
    fi
  fi
]) # _ABI_LINALG_SEARCH_BLACS



# _ABI_LINALG_SEARCH_SCALAPACK(SCALAPACK,EXTRA_LIBS)
# --------------------------------------------
#
# Look for a ScaLAPACK implementation.
#
AC_DEFUN([_ABI_LINALG_SEARCH_SCALAPACK],[
  dnl Init
  abi_linalg_has_scalapack="no"

  dnl Look for libraries and routines
  AC_SEARCH_LIBS([pzheevx],$1,
    [abi_linalg_has_scalapack="yes"],[abi_linalg_has_scalapack="no"],
    [$2 ${abi_linalg_libs}])
  if test "${abi_linalg_has_scalapack}" = "yes"; then
    if test "${ac_cv_search_pzheevx}" != "none required"; then
      abi_linalg_libs="${ac_cv_search_pzheevx} $2 ${abi_linalg_libs}"
    fi
  fi
]) # _ABI_LINALG_SEARCH_SCALAPACK



# _ABI_LINALG_SEARCH_MAGMA(MAGMA,EXTRA_LIBS)
# ------------------------------------------
#
# Look for a MAGMA implementation.
#
AC_DEFUN([_ABI_LINALG_SEARCH_MAGMA],[
  dnl Init
  abi_linalg_has_magma="no"

  dnl Look for libraries and routines
  AC_SEARCH_LIBS([magmaf_zheevd],$1,
    [abi_linalg_has_magma="yes"],[abi_linalg_has_magma="no"],
    [$2 ${abi_linalg_libs}])
  if test "${abi_linalg_has_magma}" = "yes"; then
    if test "${ac_cv_search_magmaf_zheevd}" != "none required"; then
      abi_linalg_libs="${ac_cv_search_magmaf_zheevd} $2 ${abi_linalg_libs}"
    fi
  fi
]) # _ABI_LINALG_SEARCH_MAGMA



# ABI_CONNECT_LINALG()
# --------------------
#
# Sets all variables needed to handle the optimized linear algebra
# libraries.
#
AC_DEFUN([ABI_CONNECT_LINALG],[
  dnl Initial setup
  abi_linalg_chk_gpu=""
  abi_linalg_chk_mpi=""
  abi_linalg_chk_serial=""
  abi_linalg_gpu="no"
  abi_linalg_mpi="no"
  abi_linalg_serial="no"
  abi_linalg_has_blas="no"
  abi_linalg_has_lapack="no"
  abi_linalg_has_blacs="no"
  abi_linalg_has_scalapack="no"
  abi_linalg_has_magma="no"
  abi_linalg_has_plasma="no"
  abi_linalg_incs="${with_linalg_incs}"
  abi_linalg_libs="${with_linalg_libs}"
  lib_linalg_flavor="${with_linalg_flavor}"
  lib_linalg_fcflags=""
  lib_linalg_incs=""
  lib_linalg_ldflags=""
  lib_linalg_libs=""

  dnl Prepare environment
  ABI_ENV_BACKUP
  LDFLAGS="${FC_LDFLAGS}"
  abi_saved_FCFLAGS="${FCFLAGS}"
  abi_saved_LDFLAGS="${LDFLAGS}"
  abi_saved_LIBS="${LIBS}"
  CPPFLAGS="${with_linalg_incs} ${CPPFLAGS}"
  LIBS="${with_linalg_libs} ${LIBS}"
  AC_LANG_PUSH([Fortran])

  dnl Make sure the 'none' flavor is not overriden
  if test "${with_linalg_flavor}" = "none"; then
    if test "${with_linalg_incs}" != "" -o \
            "${with_linalg_libs}" != ""; then
      AC_MSG_ERROR([user-defined linear algebra includes and libraries
                  are not allowed when the flavor is set to 'none'
           solution: use consistent linear algebra options])
    fi
  fi

  dnl Display requested flavor
  AC_MSG_CHECKING([for the requested linear algebra support])
  AC_MSG_RESULT([${lib_linalg_flavor}])

  dnl Reformat flavor
  abi_linalg_iter=`echo "${lib_linalg_flavor}" | tr '+' '\n' | sort -u | awk '{printf " %s",[$]1}'`

  dnl Check serial and parallel flavor unicity
  for abi_linalg_flavor in ${abi_linalg_iter}; do
    case "${abi_linalg_flavor}" in
      magma)
        if test "${abi_linalg_chk_gpu}" != ""; then
          AC_MSG_ERROR([only one GPU linear algebra flavor is permitted])
        fi
        abi_linalg_chk_gpu="${abi_linalg_flavor}"
        ;;
      plasma|scalapack)
        if test "${abi_linalg_chk_mpi}" != ""; then
          AC_MSG_ERROR([only one MPI linear algebra flavor is permitted])
        fi
        abi_linalg_chk_mpi="${abi_linalg_flavor}"
        ;;
      *)
        if test "${abi_linalg_chk_serial}" != ""; then
          AC_MSG_ERROR([only one serial linear algebra flavor is permitted])
        fi
        abi_linalg_chk_serial="${abi_linalg_flavor}"
        ;;
    esac
  done
  if test "${abi_linalg_chk_serial}" = ""; then
    AC_MSG_ERROR([you must choose a serial linear algebra flavor])
  fi

  dnl Check if the user has requested a fallback
  AC_MSG_CHECKING([whether to select a fallback for linear algebra])
  abi_linalg_fallback=`echo "${abi_linalg_chk_serial}" | cut -s -d- -f2`
  if test "${abi_linalg_fallback}" = "fallback"; then
    abi_linalg_fallback="yes"
  else
    abi_linalg_fallback="no"
  fi
  AC_MSG_RESULT([${abi_linalg_fallback}])
  if test "${abi_linalg_fallback}" = "yes"; then
    if test "${enable_fallbacks}" = "no"; then
      AC_MSG_ERROR([fallback requested while fallbacks have been globally disabled])
    fi
    if test "${with_linalg_incs}" != "" -o "${with_linalg_libs}" != ""; then
      AC_MSG_ERROR([you may not specify include or link flags when requesting
                  a fallback (--with-linalg-incs and --with-linalg-libs)])
    fi
  fi

  dnl Look for linear algebra libraries
  if test "${with_linalg_libs}" != "" -o \
          "${lib_linalg_flavor}" = "custom"; then
    _ABI_LINALG_CHECK_LIBS
  elif test "${lib_linalg_flavor}" != "none"; then
    case "${abi_linalg_chk_serial}" in

      acml)
        abi_linalg_fcflags=""
        abi_linalg_ldflags=""
        abi_linalg_blas_libs="acml"
        abi_linalg_blas_prqs="-lacml_mv"
        abi_linalg_lapack_libs="acml"
        abi_linalg_lapack_prqs=""
        abi_linalg_blacs_libs="acml"
        abi_linalg_blacs_prqs=""
        abi_linalg_scalapack_libs="acml"
        abi_linalg_scalapack_prqs=""
        ;;

      asl)
        abi_linalg_fcflags=""
        abi_linalg_ldflags=""
        abi_linalg_blas_libs="asl"
        abi_linalg_blas_prqs=""
        abi_linalg_lapack_libs="asl"
        abi_linalg_lapack_prqs=""
        abi_linalg_blacs_libs="asl"
        abi_linalg_blacs_prqs=""
        abi_linalg_scalapack_libs="asl"
        abi_linalg_scalapack_prqs=""
        ;;

      atlas)
        abi_linalg_fcflags=""
        abi_linalg_ldflags=""
        abi_linalg_blas_libs="f77blas"
        abi_linalg_blas_prqs="-lcblas -latlas"
        abi_linalg_lapack_libs="lapack"
        abi_linalg_lapack_prqs=""
        abi_linalg_blacs_libs="atlas"
        abi_linalg_blacs_prqs=""
        abi_linalg_scalapack_libs="atlas"
        abi_linalg_scalapack_prqs=""
        ;;

      essl)
        abi_linalg_fcflags="-qessl"
        abi_linalg_ldflags="-qessl"
        abi_linalg_blas_libs="essl"
        abi_linalg_blas_prqs=""
        abi_linalg_lapack_libs="essl"
        abi_linalg_lapack_prqs=""
        abi_linalg_blacs_libs="essl"
        abi_linalg_blacs_prqs=""
        abi_linalg_scalapack_libs="essl"
        abi_linalg_scalapack_prqs=""
        ;;

      mkl)
        abi_linalg_fcflags=""
        abi_linalg_ldflags=""
        abi_linalg_blas_libs="mkl_intel_lp64"
        abi_linalg_blas_prqs="-lmkl_sequential -lmkl_core"
        abi_linalg_lapack_libs="mkl_intel_lp64"
        abi_linalg_lapack_prqs=""
        abi_linalg_blacs_libs="mkl_blacs_lp64"
        abi_linalg_blacs_prqs=""
        abi_linalg_scalapack_libs="mkl_scalapack_lp64"
        abi_linalg_scalapack_prqs=""
        ;;

      mlib)
        abi_linalg_fcflags=""
        abi_linalg_ldflags=""
        abi_linalg_blas_libs="mlib"
        abi_linalg_blas_prqs=""
        abi_linalg_lapack_libs="mlib"
        abi_linalg_lapack_prqs=""
        abi_linalg_blacs_libs="mlib"
        abi_linalg_blacs_prqs=""
        abi_linalg_scalapack_libs="mlib"
        abi_linalg_scalapack_prqs=""
        ;;

      netlib|goto)
        abi_linalg_fcflags=""
        abi_linalg_ldflags=""
        if test "${abi_linalg_chk_serial}" = "goto"; then
          abi_linalg_blas_libs="goto"
          abi_linalg_blas_prqs=""
        else
          abi_linalg_blas_libs="blas"
          abi_linalg_blas_prqs=""
        fi
        abi_linalg_lapack_libs="lapack"
        abi_linalg_lapack_prqs=""
        abi_linalg_blacs_libs="blacs"
        abi_linalg_blacs_prqs="-lblacsCinit -lblacsF77init"
        abi_linalg_scalapack_libs="scalapack"
        abi_linalg_scalapack_prqs=""
        ;;

      sgimath)
        abi_linalg_fcflags=""
        abi_linalg_ldflags=""
        abi_linalg_blas_libs="complib.sgimath"
        abi_linalg_blas_prqs=""
        abi_linalg_lapack_libs="complib.sgimath"
        abi_linalg_lapack_prqs=""
        abi_linalg_blacs_libs="complib.sgimath"
        abi_linalg_blacs_prqs=""
        abi_linalg_scalapack_libs="complib.sgimath"
        abi_linalg_scalapack_prqs=""
        ;;

      sunperf)
        abi_linalg_fcflags=""
        abi_linalg_ldflags=""
        abi_linalg_blas_libs="sunperf"
        abi_linalg_blas_prqs=""
        abi_linalg_lapack_libs="sunperf"
        abi_linalg_lapack_prqs=""
        abi_linalg_blacs_libs="sunperf"
        abi_linalg_blacs_prqs=""
        abi_linalg_scalapack_libs="sunperf"
        abi_linalg_scalapack_prqs=""
        ;;

      *)
        if test "${abi_linalg_fallback}" = "no"; then
          AC_MSG_ERROR([unknown linear algebra flavor '${lib_linalg_flavor}'])
        fi
        ;;

    esac

    dnl MAGMA support is always separate
    abi_linalg_magma_libs="magma"
    abi_linalg_magma_prqs="${lib_gpu_libs}"

    dnl Look for the selected libraries
    if test "${abi_linalg_fallback}" = "no"; then
      FCFLAGS="${abi_saved_FCFLAGS} ${abi_linalg_fcflags}"
      LDFLAGS="${abi_saved_LDFLAGS} ${abi_linalg_ldflags}"
      _ABI_LINALG_SEARCH_BLAS([${abi_linalg_blas_libs}],[${abi_linalg_blas_prqs}])
      _ABI_LINALG_SEARCH_LAPACK([${abi_linalg_lapack_libs}],[${abi_linalg_lapack_prqs}])
      case "${abi_linalg_chk_mpi}" in
        plasma)
          AC_MSG_ERROR([not implemented yet - thank you for your patience])
          ;;
        scalapack)
          if test "${enable_mpi}" != "yes"; then
            AC_MSG_ERROR([ScaLAPACK support requires MPI])
          fi
          _ABI_LINALG_SEARCH_BLACS([${abi_linalg_blacs_libs}],[${abi_linalg_blacs_prqs}])
          _ABI_LINALG_SEARCH_SCALAPACK([${abi_linalg_scalapack_libs}],[${abi_linalg_scalapack_prqs}])
          ;;
        *)
          if test "${abi_linalg_chk_mpi}" != ""; then
            AC_MSG_ERROR([library search for ${abi_linalg_chk_mpi} not implemented])
          fi
          ;;
      esac
      case "${abi_linalg_chk_gpu}" in
        magma)
          if test "${enable_gpu}" != "yes"; then
            AC_MSG_ERROR([MAGMA requires GPU support])
          fi
          _ABI_LINALG_SEARCH_MAGMA([${abi_linalg_magma_libs}],[${abi_linalg_magma_prqs}])
          ;;
        *)
          if test "${abi_linalg_chk_gpu}" != ""; then
            AC_MSG_ERROR([library search for ${abi_linalg_chk_gpu} not implemented])
          fi
          ;;
      esac
    fi
  fi

  dnl Set serial, MPI and GPU status
  if test "${abi_linalg_has_blas}" = "yes" -a \
          "${abi_linalg_has_lapack}" = "yes"; then
    abi_linalg_serial="yes"
    if test "${abi_linalg_has_blacs}" = "yes" -a \
            "${abi_linalg_has_scalapack}" = "yes"; then
      abi_linalg_mpi="yes"
    fi
    if test "${abi_linalg_has_magma}" = "yes"; then
      abi_linalg_gpu="yes"
    fi
  fi

  dnl Transmit serial status to the source code
  AC_MSG_CHECKING([whether we have a serial linear algebra support])
  AC_MSG_RESULT([${abi_linalg_serial}])
  if test "${abi_linalg_serial}" = "yes"; then
    AC_DEFINE([HAVE_LINALG],1,[Define to 1 if you have an optimized linear algebra library.])
    AC_DEFINE([HAVE_LINALG_SERIAL],1,[Define to 1 if you have an optimized serial linear algebra library.])

    case "${abi_linalg_chk_serial}" in
      asl)
        AC_DEFINE([HAVE_LINALG_ASL],1,[Define to 1 if you have the ASL linear algebra library.])
        ;;
      essl)
        AC_DEFINE([HAVE_LINALG_ESSL],1,[Define to 1 if you have the ESSL linear algebra library.])
        ;;
    esac

    lib_linalg_fcflags="${abi_linalg_fcflags}"
    lib_linalg_ldflags="${abi_linalg_ldflags}"
    lib_linalg_incs="${abi_linalg_incs}"
    lib_linalg_libs="${abi_linalg_libs}"
  else
    lib_linalg_flavor="broken"
    AC_MSG_WARN([falling back to internal linear algebra libraries])
    abi_fallbacks="${abi_fallbacks} linalg"
    lib_linalg_flavor="netlib-fallback"
    abi_dft_linalg_fallback="yes"
  fi

  dnl Transmit MPI status to the source code
  AC_MSG_CHECKING([whether we have a MPI linear algebra support])
  AC_MSG_RESULT([${abi_linalg_mpi}])
  if test "${abi_linalg_mpi}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_MPI],1,[Define to 1 if you have an optimized MPI-parallel linear algebra library.])
    case "${abi_linalg_chk_mpi}" in
      #plasma)
        #AC_DEFINE([HAVE_LINALG_PLASMA],1,[Define to 1 if you have an optimized PLASMA linear algebra library.])
        #;;
      scalapack)
        AC_DEFINE([HAVE_LINALG_SCALAPACK],1,[Define to 1 if you have an optimized ScaLAPACK linear algebra library.])
        ;;
    esac
  elif test "${abi_linalg_chk_mpi}" != ""; then
    lib_linalg_flavor="broken"
  fi

  dnl Transmit GPU status to the source code
  AC_MSG_CHECKING([whether we have a GPU linear algebra support])
  AC_MSG_RESULT([${abi_linalg_gpu}])
  if test "${abi_linalg_gpu}" = "yes"; then
    AC_DEFINE([HAVE_LINALG_GPU],1,[Define to 1 if you have an optimized GPU-compatible linear algebra library.])
    case "${abi_linalg_chk_gpu}" in
      magma)
        AC_DEFINE([HAVE_LINALG_MAGMA],1,[Define to 1 if you have the MAGMA linear algebra library.])
        ;;
    esac
  elif test "${abi_linalg_chk_gpu}" != ""; then
    lib_linalg_flavor="broken"
  fi

  dnl Restore build environment
  AC_LANG_POP([Fortran])
  FCFLAGS="${abi_saved_FCFLAGS}"
  LDFLAGS="${abi_saved_LDFLAGS}"
  LIBS="${abi_saved_LIBS}"
  ABI_ENV_RESTORE

  dnl Output final flavor
  AC_MSG_CHECKING([for the actual linear algebra support])
  AC_MSG_RESULT([${lib_linalg_flavor}])
  if test "${lib_linalg_flavor}" = "broken"; then
    ABI_MSG_NOTICE([connectors-failure],[Connector detection failure])
    AC_MSG_ERROR([the requested ${with_linalg_flavor} linear algebra flavor is not supported on this architecture])
  fi

  dnl Substitute variables needed for the use of the library
  AC_SUBST(lib_linalg_flavor)
  AC_SUBST(lib_linalg_fcflags)
  AC_SUBST(lib_linalg_ldflags)
  AC_SUBST(lib_linalg_incs)
  AC_SUBST(lib_linalg_libs)
]) # ABI_CONNECT_LINALG
