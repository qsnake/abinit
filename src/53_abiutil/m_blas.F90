!{\src2tex{textfont=tt}}
!!****f* ABINIT/m_blas
!! NAME
!!  m_blas
!!
!! FUNCTION
!! This module defines interfaces for overloading BLAS routines. 
!! whose goal is twofold. On one hand, using generic interfaces renders 
!! the code more readable, especially when the routine can be compiled with 
!! different precision type (single-precision or double precision as done for example in the GW code)
!! On the other hand, the generic interfaces defined here introduce a programming 
!! layer that can be exploited for interfacing non-standard libraries such as for 
!! example CUBLAS routines for GPU computations.
!!
!! COPYRIGHT
!! Copyright (C) 1992-2012 ABINIT group (M.Giantomassi)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!!
!! The convention about names of interfaced routine is: x<name>, 
!! where <name> is usually equal to the name of the standard routine
!! without the first character specifying the type and kind.
!! The full list of names is reported below.
!! BLAS procedures interfaced in this module are marked with an asterisk.
!! A complete list of possible overloaded interfaces is provided as guide for future additions.
!!
!! ================
!! ==== BLAS 1 ====
!! ================
!! FUNCTION idamax isamax icamax izamax  ---> XIAMAX(n,dx,incx)
!! * FUNCTION  snrm2  dnrm2 scnrm2 dznmr2  ---> XNRM2(n,x,incx)
!! FUNCTION  sasum  dasum scasum dzasum  ---> XASUM(n,x,incx)
!! * FUNCTION               cdotu  zdotu ---> XDOTU(n,x,incx,y,incy)
!! * FUNCTION               cdotc  zdotc ---> XDOTC(n,x,incx,y,incy)
!! FUNCTION  sdot   ddot                 ---> XDOT(n,x,incx,y,incy)
!! FUNCTION  sdsdot sdot                 ---> XDSDOT(n,x,incx,y,incy)
!! SUBROUTINE saxpy daxpy caxpy  zaxpy   ---> XAXPY(n,ca,cx,incx,cy,incy)
!! * SUBROUTINE scopy dcopy ccopy  zcopy   ---> XCOPY(n,cx,incx,cy,incy)
!! SUBROUTINE srotg drotg crotg  zrotg   ---> XROTG(a,b,c,s)
!! SUBROUTINE srot  drot  csrot  zdrot   ---> XROT(n,x,incx,y,incy,c,s)
!! SUBROUTINE sscal dscal cscal  zscal  
!!                        csscal zdscal  ---> XSCAL(n,a,x,incx)
!! SUBROUTINE sswap dswap cswap  zswap   ---> XSWAP(n,x,incx,y,incy)
!!
!! ================
!! ==== BLAS 2 ====
!! ================
!! SUBROUTINE sgbmv dgbmv cgbmv zgbmv    ---> XGBMV(trans,m,kl,ku,n,alpha,A,lda,X,incx,beta,Y,incy)
!! * SUBROUTINE sgemv dgemv cgemv zgemv    ---> XGEMV(trans,m,n,alpha,A,lda,X,incx,beta,Y,incy)
!! * SUBROUTINE             cgerc zgerc    ---> XGERC(m,n,alpha,x,incx,y,incy,A,lda)
!! SUBROUTINE             cgeru zgeru    ---> XGERU(m,n,alpha,x,incx,y,incy,A,lda)
!! SUBROUTINE             chbmv zhbmv    ---> XHBMV(uplo,n,k,alpha,A,lda,X,incx,beta,Y,incy)
!! SUBROUTINE             chemv zhemv    ---> XHEMV(uplo,n,alpha,A,lda,X,incx,beta,Y,incy)
!! SUBROUTINE             cher  zher     ---> XHER(uplo,n,alpha,x,incx,A,lda)
!! SUBROUTINE             cher2 zher2    ---> XHER2(uplo,n,alpha,x,incx,y,incy,A,lda)
!! SUBROUTINE             chpr  zhpr     ---> XHPR(uplo,n,alpha,x,incx,AP)
!! SUBROUTINE             chpr2 zhpr2    ---> XHPR2(uplo,n,alpha,x,incx,y,incy,AP)
!! SUBROUTINE             chpmv zhpmv    ---> XHPMV(uplo,n,alpha,AP,X,incx,beta,Y,incy)
!! SUBROUTINE stbmv dtbmv ctbmv ztbmv    ---> XTBMV(uplo,trans,diag,n,k,A,lda,X,incx)
!! SUBROUTINE stpmv dtpmv ctpmv ztpmv    ---> XTPMV(uplo,trans,diag,n,AP,X,incx)
!! SUBROUTINE strmv dtrmv ctrmv ztrmv    ---> XTRMV(uplo,trans,diag,n,A,lda,X,incx)
!! SUBROUTINE ssymv dsymv                ---> XSYMV(uplo,n,alpha,A,lda,X,incx,beta,Y,incy)
!! SUBROUTINE ssbmv dsbmv                ---> XSBMV(uplo,n,k,alpha,A,lda,X,incx,beta,Y,incy)
!! SUBROUTINE sspmv dspmv                ---> XSPMV(uplo,n,alpha,AP,X,incx,beta,Y,incy)
!! SUBROUTINE stbsv dtbsv ctbsv ztbsv    ---> XTBSV(uplo,trans,diag,n,k,A,lda,X,incx)
!! SUBROUTINE stpsv dtpsv ctpsv ztpsv    ---> XTPSV(uplo,trans,diag,n,AP,X,incx)
!! SUBROUTINE strsv dtrsv ctrsv ztrsv    ---> XTRSV(uplo,trans,diag,n,A,lda,X,incx)
!! SUBROUTINE  sger  dger                ---> XGER(m,n,alpha,x,incx,y,incy,A,lda)
!! SUBROUTINE  sspr  dspr                ---> XSPR(uplo,n,alpha,x,incx,AP)
!! SUBROUTINE sspr2 dspr2                ---> XSPR2(uplo,n,alpha,x,incx,y,incy,AP)
!! SUBROUTINE  ssyr  dsyr                ---> XSYR(uplo,n,alpha,x,incx,A,lda)
!! SUBROUTINE ssyr2 dsyr2                ---> XSYR2(uplo,n,alpha,x,incx,y,incy,A,lda)
!!
!! ================
!! ==== BLAS 3 ====
!! ================
!! * SUBROUTINE sgemm dgemm cgemm zgemm      ---> XGEMM(transA,transB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc)
!! SUBROUTINE             chemm zhemm      ---> XHEMM(side,uplo,m,n,alpha,A,lda,B,ldb,beta,C,ldc)
!! SUBROUTINE            cher2k zher2k     ---> XHER2K(uplo,trans,n,k,alpha,A,lda,B,ldb,beta,C,ldc)
!! SUBROUTINE             cherk zherk      ---> XHERK(uplo,trans,n,k,alpha,A,lda,beta,C,ldc)
!! SUBROUTINE ssymm dsymm csymm zsymm      ---> XSYMM(side,uplo,m,n,alpha,A,lda,B,ldb,beta,C,ldc)
!! SUBROUTINE ssyr2k dsyr2k csyr2k zsyr2k  ---> XSYR2K(uplo,trans,n,k,alpha,A,lda,B,ldb,beta,C,ldc)
!! SUBROUTINE ssyrk dsyrk csyrk zsyrk      ---> XSYRK(uplo,trans,n,k,alpha,A,lda,beta,C,ldc)
!! SUBROUTINE strmm dtrmm ctrmm ztrmm      ---> XTRMM(side,uplo,transa,diag,m,n,alpha,A,lda,B,ldb)
!! SUBROUTINE strsm dtrsm ctrsm ztrsm      ---> XTRSM(side,uplo,transa,diag,m,n,alpha,A,lda,B,ldb)
!!-------------------------------------------------------------------------------
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h" 

#define HAVE_LINALG_ZDOTC_BUG
#define HAVE_LINALG_ZDOTU_BUG

MODULE m_blas

 use m_profiling

 use defs_basis
 use m_errors

 implicit none

 private

 public :: xnrm2
 public :: xdotu
 public :: xdotc
 public :: xcopy
 public :: xgemv
 public :: xgerc
 public :: xgemm
 public :: blas_cholesky_ortho

!----------------------------------------------------------------------

interface xnrm2  
  !
  function snrm2 ( n, x, incx )
    use defs_basis
    real(sp) ::  snrm2
    integer,intent(in) :: incx, n
    real(sp),intent(in) ::  x( * )
  end function snrm2
  !
  function dnrm2 ( n, x, incx )
    use defs_basis
    real(dp) :: dnrm2
    integer,intent(in) :: incx, n
    real(dp),intent(in) ::  x( * )
  end function dnrm2
  !
  function scnrm2( n, x, incx )
    use defs_basis
    real(sp) :: scnrm2
    integer,intent(in) :: incx, n
    complex(spc),intent(in) :: x( * )
  end function scnrm2
  !
  function dznrm2( n, x, incx )
    use defs_basis
    real(dp) :: dznrm2
    integer,intent(in) :: incx, n
    complex(dpc),intent(in) :: x( * )
  end function dznrm2
  !
end interface xnrm2
!***

!-------------------------------------------------------------------------------

interface xdotu
  !
#ifdef HAVE_LINALG_ZDOTU_BUG 
  module procedure cdotu
  module procedure zdotu
#else
  function cdotu(n,cx,incx,cy,incy)
    use defs_basis
    complex(spc) :: cdotu
    complex(spc),intent(in) :: cx(*),cy(*)
    integer,intent(in) :: incx,incy,n
  end function cdotu
  !
  function zdotu(n,zx,incx,zy,incy)
    use defs_basis
    complex(dpc) :: zdotu
    complex(dpc),intent(in) :: zx(*),zy(*)
    integer,intent(in) :: incx,incy,n
  end function zdotu
#endif
  !
end interface xdotu
!***

!-------------------------------------------------------------------------------


! CDOTC, CDOTU, ZDOTC, and ZDOTU are problematic if Mac OS X's Vec lib is used.
! See http://developer.apple.com/hardwaredrivers/ve/errata.html.
! If needed, we replace them with plain Fortran code.

interface xdotc
  !
#ifdef HAVE_LINALG_ZDOTC_BUG
   module procedure cdotc
   module procedure zdotc
#else
  function cdotc(n,cx,incx,cy,incy)
    use defs_basis
    complex(spc) :: cdotc
    complex(spc),intent(in) :: cx(*),cy(*)
    integer,intent(in) :: incx,incy,n
  end function cdotc
  !
  function zdotc(n,zx,incx,zy,incy)
    use defs_basis
    complex(dpc) :: zdotc
    complex(dpc),intent(in) :: zx(*),zy(*)
    integer,intent(in) :: incx,incy,n
  end function zdotc
#endif
  !
end interface xdotc
!***


!-------------------------------------------------------------------------------

interface xcopy
 !
 subroutine scopy(n,sx,incx,sy,incy)
   use defs_basis
   implicit none
   integer,intent(in) :: incx
   integer,intent(in) :: incy
   integer,intent(in) :: n
   real(sp),intent(in) ::  sx(*)
   real(sp),intent(inout) :: sy(*)
 end subroutine scopy
 !
 subroutine  dcopy(n,dx,incx,dy,incy)
   use defs_basis
   implicit none
   integer,intent(in) :: incx
   integer,intent(in) :: incy
   integer,intent(in) :: n
   real(dp),intent(in) :: dx(*)
   real(dp),intent(inout) :: dy(*)
 end subroutine dcopy
 !
 subroutine  ccopy(n,cx,incx,cy,incy)
   use defs_basis
   implicit none
   integer,intent(in) :: incx
   integer,intent(in) :: incy
   integer,intent(in) :: n
   complex(spc),intent(in) :: cx(*)
   complex(spc),intent(inout) :: cy(*)
 end subroutine ccopy
 !
 subroutine  zcopy(n,cx,incx,cy,incy)
   use defs_basis
   implicit none
   integer,intent(in) :: incx
   integer,intent(in) :: incy
   integer,intent(in) :: n
   complex(dpc),intent(in) :: cx(*)
   complex(dpc),intent(inout) :: cy(*)
 end subroutine zcopy
 !
end interface xcopy
!***

!-------------------------------------------------------------------------------

interface xgemv
  !
  subroutine sgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )
    use defs_basis
    real(sp),intent(in) :: alpha, beta
    integer,intent(in) :: incx, incy, lda, m, n
    character(len=1),intent(in) :: trans
    real(sp),intent(in) :: a( lda, * ), x( * )
    real(sp),intent(inout) :: y( * )
  end subroutine sgemv
  !
  subroutine dgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )
    use defs_basis
    real(dp),intent(in) :: alpha, beta
    integer,intent(in) :: incx, incy, lda, m, n
    character(len=1),intent(in) :: trans
    real(dp),intent(in) :: a( lda, * ), x( * )
    real(dp),intent(inout) :: y( * )
  end subroutine dgemv
  !
  subroutine cgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )
    use defs_basis
    complex(spc),intent(in) :: alpha, beta
    integer,intent(in) :: incx, incy, lda, m, n
    character(len=1),intent(in) :: trans
    complex(spc),intent(in) :: a( lda, * ), x( * )
    complex(spc),intent(inout) :: y( * )
  end subroutine cgemv
  !
  subroutine zgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )
    use defs_basis
    complex(dpc),intent(in) :: alpha, beta
    integer,intent(in) :: incx, incy, lda, m, n
    character(len=1),intent(in) :: trans
    complex(dpc),intent(in) :: a( lda, * ), x( * )
    complex(dpc),intent(inout) :: y( * )
  end subroutine zgemv
  !
end interface xgemv
!***

!-------------------------------------------------------------------------------

interface xgerc
  !
  subroutine cgerc ( m, n, alpha, x, incx, y, incy, a, lda )
    use defs_basis
    complex(spc),intent(in) :: alpha
    integer,intent(in) :: incx, incy, lda, m, n
    complex(spc),intent(inout) ::  a( lda, * )
    complex(spc),intent(in) :: x( * ), y( * )
  end subroutine cgerc
  !
  subroutine zgerc ( m, n, alpha, x, incx, y, incy, a, lda )
    use defs_basis
    complex(dpc),intent(in) :: alpha
    integer,intent(in) :: incx, incy, lda, m, n
    complex(dpc),intent(inout) :: a( lda, * )
    complex(dpc),intent(in) :: x( * ), y( * )
  end subroutine zgerc
  !
end interface xgerc
!***

!-------------------------------------------------------------------------------

interface xgemm
  !
  subroutine sgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
    use defs_basis
    character(len=1),intent(in) :: transa, transb
    integer,intent(in) :: m, n, k, lda, ldb, ldc
    real(sp),intent(in) :: alpha, beta
    real(sp),intent(in) :: a( lda, * ), b( ldb, * )
    real(sp),intent(inout) :: c( ldc, * )
  end subroutine sgemm
  !
  subroutine dgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
    use defs_basis
    character(len=1),intent(in) :: transa, transb
    integer,intent(in) :: m, n, k, lda, ldb, ldc
    real(dp),intent(in) :: alpha, beta
    real(dp),intent(in) :: a( lda, * ), b( ldb, * )
    real(dp),intent(inout) :: c( ldc, * )
  end subroutine dgemm
  !
  subroutine cgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
    use defs_basis
    character(len=1),intent(in) :: transa, transb
    integer,intent(in) :: m, n, k, lda, ldb, ldc
    complex(spc),intent(in) :: alpha, beta
    complex(spc),intent(in) :: a( lda, * ), b( ldb, * )
    complex(spc),intent(inout) :: c( ldc, * )
  end subroutine cgemm
  !
  subroutine zgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
    use defs_basis
    character(len=1),intent(in) :: transa, transb
    integer,intent(in) :: m, n, k, lda, ldb, ldc
    complex(dpc),intent(in) :: alpha, beta
    complex(dpc),intent(in) :: a( lda, * ), b( ldb, * )
    complex(dpc),intent(inout) :: c( ldc, * )
  end subroutine zgemm
  !
end interface xgemm
!!***

!-------------------------------------------------------------------------------

interface blas_cholesky_ortho      
  module procedure blas_cholesky_ortho_spc
  module procedure blas_cholesky_ortho_dpc
end interface blas_cholesky_ortho
!***

 complex(spc),private,parameter :: czero_spc =(0._sp,0._sp)
 complex(spc),private,parameter :: cone_spc  =(1._sp,0._sp)

 complex(dpc),private,parameter :: czero_dpc =(0._dp,0._dp)
 complex(dpc),private,parameter :: cone_dpc  =(1._dp,0._dp)

CONTAINS  !========================================================================================
!!***

! CDOTC, CDOTU, ZDOTC, and ZDOTU are problematic if Mac OS X's Vec lib is used.
! See http://developer.apple.com/hardwaredrivers/ve/errata.html.
! Here we replace them with plain Fortran code.

#ifdef HAVE_LINALG_ZDOTC_BUG
!#warning "Using internal replacement for zdotc. External library cannot be used"
#include "replacements/cdotc.f"
#include "replacements/zdotc.f"
#endif

#ifdef HAVE_LINALG_ZDOTU_BUG
!#warning "Using internal replacement for zdotu. External library cannot be used"
#include "replacements/cdotu.f"
#include "replacements/zdotu.f"
#endif

!----------------------------------------------------------------------

!!****f* m_blas/blas_cholesky_ortho_spc  
!! NAME
!!  blas_cholesky_ortho_spc
!!
!! FUNCTION
!!  Performs the Cholesky orthonormalization of the vectors stored in iomat.
!!
!! INPUTS
!!  vec_size=Size of each vector.
!!  nvec=Number of vectors in iomat
!!
!! OUTPUT
!!  cf_ovlp=Cholesky factorization of the overlap matrix. ovlp = U^H U with U upper triangle matrix returned in cf_ovlp
!!
!! SIDE EFFECTS
!!  iomat(vec_size,nvec)
!!    input: Input set of vectors.
!!    output: Orthonormalized set.
!!
!! PARENTS
!!
!! CHILDREN
!!      xgemm,zpotrf,ztrsm
!!
!! SOURCE

subroutine blas_cholesky_ortho_spc(vec_size,nvec,iomat,cf_ovlp)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'blas_cholesky_ortho_spc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: vec_size,nvec
 complex(spc),intent(inout) :: iomat(vec_size,nvec)
 complex(spc),intent(out) :: cf_ovlp(nvec,nvec)

!Local variables ------------------------------
!scalars
 integer :: ierr
 character(len=500) :: msg

! *************************************************************************

 ! 1) Calculate overlap_ij =  <phi_i|phi_j>
 ! TODO: use dsyrk
 call xgemm("Conjugate","Normal",nvec,nvec,vec_size,cone_spc,iomat,vec_size,iomat,vec_size,czero_spc,cf_ovlp,nvec)
 !
 ! 2) Cholesky factorization: ovlp = U^H U with U upper triangle matrix.
 call CPOTRF('U',nvec,cf_ovlp,nvec,ierr)
 if (ierr/=0)  then
   write(msg,'(a,i0)')' ZPOTRF returned info= ',ierr
   MSG_ERROR(msg)
 end if 
 !
 ! 3) Solve X U = io_mat. On exit io_mat is orthonormalized.
 call CTRSM('Right','Upper','Normal','Normal',vec_size,nvec,cone_spc,cf_ovlp,nvec,iomat,vec_size)

end subroutine blas_cholesky_ortho_spc
!!***

!----------------------------------------------------------------------

!!****f* m_blas/blas_cholesky_ortho_dpc  
!! NAME
!!  blas_cholesky_ortho_dpc
!!
!! FUNCTION
!!  Performs the Cholesky orthonormalization of the vectors stored in iomat.
!!
!! INPUTS
!!  vec_size=Size of each vector.
!!  nvec=Number of vectors in iomat
!!
!! OUTPUT
!!  cf_ovlp=Cholesky factorization of the overlap matrix. ovlp = U^H U with U upper triangle matrix returned in cf_ovlp
!!
!! SIDE EFFECTS
!!  iomat(vec_size,nvec)
!!    input: Input set of vectors.
!!    output: Orthonormalized set.
!!
!! PARENTS
!!
!! CHILDREN
!!      xgemm,zpotrf,ztrsm
!!
!! SOURCE

subroutine blas_cholesky_ortho_dpc(vec_size,nvec,iomat,cf_ovlp)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'blas_cholesky_ortho_dpc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: vec_size,nvec
 complex(dpc),intent(inout) :: iomat(vec_size,nvec)
 complex(dpc),intent(out) :: cf_ovlp(nvec,nvec)

!Local variables ------------------------------
!scalars
 integer :: ierr
 character(len=500) :: msg

! *************************************************************************

 ! 1) Calculate overlap_ij =  <phi_i|phi_j>
 call xgemm("Conjugate","Normal",nvec,nvec,vec_size,cone_dpc,iomat,vec_size,iomat,vec_size,czero_dpc,cf_ovlp,nvec)
 !
 ! 2) Cholesky factorization: ovlp = U^H U with U upper triangle matrix.
 call ZPOTRF('U',nvec,cf_ovlp,nvec,ierr)
 if (ierr/=0)  then
   write(msg,'(a,i0)')' ZPOTRF returned info= ',ierr
   MSG_ERROR(msg)
 end if 
 !
 ! 3) Solve X U = io_mat. On exit io_mat is orthonormalized.
 call ZTRSM('Right','Upper','Normal','Normal',vec_size,nvec,cone_dpc,cf_ovlp,nvec,iomat,vec_size)

end subroutine blas_cholesky_ortho_dpc
!!***

!----------------------------------------------------------------------

END MODULE m_blas
!!***
