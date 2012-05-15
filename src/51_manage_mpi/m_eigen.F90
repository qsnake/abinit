!{\src2tex{textfont=tt}}
!!****f* ABINIT/m_eigen
!! NAME
!!  m_eigen
!!
!! FUNCTION
!!  management of eigen problems computation
!!  with or without ScaLAPACK
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2012 ABINIT group (FBottin,CS)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~ABINIT/Infos/copyright
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_eigen

 use m_profiling

 use defs_basis
 use defs_scalapack
 use m_wfutils

 implicit none

 private

!Constants ------------------------------------
 logical,public :: x_withsclapack

#if defined HAVE_LINALG_SCALAPACK
!Scalapack variables ---------------------------
 type(processor_scalapack),public :: x_processor
 integer,public                   :: x_communicator
#endif

!Work arrays for eigen problem -----------------
 real(dp), dimension(:), allocatable :: work,rwork
 integer, dimension(:), allocatable :: iwork

!Procedures ------------------------------------
 public :: xsetscalapack ! Define use of scalapack or not
 public :: xeigeninit    ! Initialisation of eigen computing
 public :: xeigenend     ! Terminate eigen computations
 public :: xeigen1       ! Computes all eigenvalues & eigenvectors
 public :: xeigen2       ! Computes all eigenvalues & eigenvectors
 public :: xev           ! Computes all eigenvalues & eigenvectors
 public :: xgv           ! Computes all eigenvalues & eigenvectors

contains

!! NAME
!!
!! xsetscalapack
!!
!! FUNCTION
!! Define use of scalapack or not
!!
!! COPYRIGHT
!!   Copyright (C) 1998-2012 ABINIT group (FBottin,CS)
!!   this file is distributed under the terms of the
!!   gnu general public license, see ~abinit/COPYING
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   withscalapack true if using scalapack or false  without using scalapack
!!
!! PARENTS
!!      lobpcgwf
!!
!! CHILDREN
!!
!! SOURCE
!!
 subroutine xsetscalapack(withscalapack)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xsetscalapack'
!End of the abilint section

 logical, intent(in) :: withscalapack

! *********************************************************************

 x_withsclapack=withscalapack
#if !(defined HAVE_LINALG_SCALAPACK)
 x_withsclapack=.false.
#endif


 end subroutine xsetscalapack
!!***

!! NAME
!!   xeigeninit
!!
!! FUNCTION
!!   Initialisation of eigen computing
!!   Do nothing if not using scalapack
!!
!! COPYRIGHT
!!   Copyright (C) 1998-2012 ABINIT group (FBottin,CS,FDahm)
!!   this file is distributed under the terms of the
!!   gnu general public license, see ~abinit/COPYING
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   paralbd
!!   commcart:         communicator for the full cartesian array
!!   num_group:        number of group of my processor. 0 if my processor is not in a group
!!
!! PARENTS
!!      lobpcgwf
!!
!! CHILDREN
!!      dcopy,dpotrf,dsygst,dtrmm,dtrsm,magmaf_dsyevd,magmaf_zheevd,xerbla
!!      zhegst,zpotrf,ztrmm,ztrsm
!!
!! SOURCE
!!
 subroutine xeigeninit(paralbd,commcart,num_group,cplx,maxsize,usegpu)

#if defined HAVE_MPI2
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xeigeninit'
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer, intent(in) :: commcart,cplx,maxsize,paralbd
 integer, intent(in) :: num_group,usegpu
!Local variables ------------------------------
 integer :: lwork,lrwork,liwork

! *********************************************************************

 if (usegpu==1) then
   lwork=33*maxsize + (maxsize**2); !work only if  lwork=n**2 + 33*n?
   if(cplx==1) lwork= 1+maxsize*(6 + 2*maxsize)
   rwork= 1 + 5*maxsize + 2*(maxsize**2)
   liwork= 3 + 5*maxsize
   ABI_ALLOCATE(work,(cplx*lwork))
   ABI_ALLOCATE(rwork,(lrwork))
   ABI_ALLOCATE(iwork,(liwork))
 else
#if defined HAVE_LINALG_SCALAPACK
   if ( x_withsclapack ) then
     if (paralbd <=1 ) then
       x_communicator = commcart
     else
       x_communicator = num_group
     endif
     call init_scalapack(x_processor,x_communicator)
   else
#endif
     lwork=3*maxsize - (cplx-1)*2
     lrwork=lwork
     ABI_ALLOCATE(work,(cplx*lwork))
     ABI_ALLOCATE(rwork,(lrwork))
#if defined HAVE_LINALG_SCALAPACK
   end if
#endif
 end if

#if ! defined HAVE_LINALG_SCALAPACK
 if (.false.) write(std_out,*) paralbd,commcart,num_group
#endif


 end subroutine xeigeninit
!!***

!! NAME
!!   xeigenend
!!
!! FUNCTION
!!   Terminate eigen computations
!!   Do nothing if not using scalapack
!!
!! COPYRIGHT
!!   Copyright (C) 1998-2012 ABINIT group (FBottin,CS,FDahm)
!!   this file is distributed under the terms of the
!!   gnu general public license, see ~abinit/COPYING
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! PARENTS
!!      lobpcgwf
!!
!! CHILDREN
!!      dcopy,dpotrf,dsygst,dtrmm,dtrsm,magmaf_dsyevd,magmaf_zheevd,xerbla
!!      zhegst,zpotrf,ztrmm,ztrsm
!!
!! SOURCE
!!
 subroutine xeigenend()

!Arguments ------------------------------------

! *********************************************************************

# if defined HAVE_LINALG_SCALAPACK

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xeigenend'
!End of the abilint section

 if ( x_withsclapack ) then
   call end_scalapack(x_processor)
 else
#endif
   ABI_DEALLOCATE(work)
   ABI_DEALLOCATE(rwork)
   if(allocated(iwork))  then
     ABI_DEALLOCATE(iwork)
   end if
#if defined HAVE_LINALG_SCALAPACK
 end if
#endif


 end subroutine xeigenend
!!***

!! NAME
!!   xeigen1
!!
!! FUNCTION
!!   Computes all eigenvalues and, optionally, eigenvectors
!!   with or without scalapack
!!
!! COPYRIGHT
!!   Copyright (C) 1998-2012 ABINIT group (FBottin,CS,FDahm)
!!   this file is distributed under the terms of the
!!   gnu general public license, see ~abinit/COPYING
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!!   cplx
!!   blocksize
!!   gramxax
!!   eigen
!!   istwf_k
!!   info
!!
!! PARENTS
!!      lobpcgwf
!!
!! CHILDREN
!!      dcopy,dpotrf,dsygst,dtrmm,dtrsm,magmaf_dsyevd,magmaf_zheevd,xerbla
!!      zhegst,zpotrf,ztrmm,ztrsm
!!
!! SOURCE
!!
 subroutine xeigen1(cplx,blocksize,gramxax,eigen,istwf_k,info,&
 &                  timopt,tim_xeigen,usegpu) ! optional arguments

 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xeigen1'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) ::cplx
 integer, intent(in) ::blocksize
 integer, intent(in), optional :: timopt,tim_xeigen,usegpu
 real(dp), intent(inout) :: gramxax(cplx*blocksize,blocksize)
 real(dp), intent(inout) ::eigen(blocksize)
 integer, intent(in) ::istwf_k
 integer, intent(inout) ::info
!Local variables-------------------------------
 integer :: lwork,usegpu_
 real(dp) :: tsec(2)
 character(len=500) :: msg
#if defined HAVE_LINALG_MAGMA
 integer :: lrwork,liwork
#endif

! *********************************************************************

 usegpu_=0;if (present(usegpu)) usegpu_=usegpu
 if (present(tim_xeigen).and.present(timopt)) then
   if(abs(timopt)==3) call timab(tim_xeigen,1,tsec)
 end if

#if defined HAVE_LINALG_MAGMA
 if (usegpu_==1) then
   lwork=33*blocksize + (blocksize**2); !work only if  lwork=n**2 + 33*n?
   if(cplx==1) lwork= 1+blocksize*(6 + 2*blocksize)
   lrwork= 1 + 5*blocksize + 2*(blocksize**2)
   liwork= 3 + 5*blocksize
   call xev_gpu('v','u',blocksize,gramxax,blocksize,eigen,work(1:cplx*lwork),&
&               lwork,rwork(1:lrwork),lrwork,iwork(1:liwork),liwork,info)
   if (info/=0)  then
     write(msg,'(a,i3)')' Problem in xev_gpu, info= ',info
     MSG_WARNING(msg)
   endif
 else
#endif
#if defined HAVE_LINALG_SCALAPACK
   if ( x_withsclapack ) then
     call compute_eigen1(x_communicator,x_processor,blocksize,blocksize,gramxax,eigen,istwf_k)
   else
#endif
     lwork=3*blocksize - (cplx-1)*2
     call xev('v','u',blocksize,gramxax,blocksize,eigen,work(1:cplx*lwork),&
&             lwork,rwork(1:lwork),info)
     if (info/=0)  then
       write(msg,'(a,i3)')' Problem in xev, info= ',info
       MSG_WARNING(msg)
     endif
#if defined HAVE_LINALG_SCALAPACK
   endif
# endif
#if defined HAVE_LINALG_MAGMA
 endif
#endif

 if (present(tim_xeigen).and.present(timopt)) then
   if(abs(timopt)==3) call timab(tim_xeigen,2,tsec)
 end if

#if !defined HAVE_LINALG_SCALAPACK
 if (.false.) write(std_out,*) istwf_k
#endif

 end subroutine xeigen1


!!***

!! NAME
!!   xeigen2
!!
!! FUNCTION
!!   Computes all eigenvalues and, optionally, eigenvectors
!!   with or without scalapack
!!
!! COPYRIGHT
!!   Copyright (C) 1998-2012 ABINIT group (FBottin,CS,FDahm)
!!   this file is distributed under the terms of the
!!   gnu general public license, see ~abinit/COPYING
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   cplx
!!   blocksize
!!   gramxax
!!   eigen
!!   istwf_k
!!   info
!!
!! PARENTS
!!      lobpcgwf
!!
!! CHILDREN
!!      dcopy,dpotrf,dsygst,dtrmm,dtrsm,magmaf_dsyevd,magmaf_zheevd,xerbla
!!      zhegst,zpotrf,ztrmm,ztrsm
!!
!! SOURCE
!!
 subroutine xeigen2(cplx,blocksize,grama,gramb,eigen,istwf_k,info,&
     &             timopt,tim_xeigen,usegpu) ! optional arguments

 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xeigen2'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) ::cplx
 integer, intent(in) ::blocksize
 integer, intent(in), optional :: timopt,tim_xeigen,usegpu
 real(dp), intent(inout) :: grama(cplx*blocksize,blocksize)
 real(dp), intent(inout) :: gramb(cplx*blocksize,blocksize)
 real(dp), intent(inout) ::eigen(blocksize)
 integer, intent(in) ::istwf_k
 integer, intent(inout) ::info
!Local variables-------------------------------
 integer ::lwork,usegpu_
 real(dp) :: tsec(2)
 character(len=500) :: msg
#if defined HAVE_LINALG_MAGMA
 integer :: lrwork,liwork
#endif

! *********************************************************************

 usegpu_=0;if (present(usegpu)) usegpu_=usegpu
 if (present(tim_xeigen).and.present(timopt)) then
   if(abs(timopt)==3) call timab(tim_xeigen,1,tsec)
 end if

#if defined HAVE_LINALG_MAGMA
 if (usegpu_==1) then
   lwork=33*blocksize + (blocksize**2); !work only if  lwork=n**2 + 33*n?
   if(cplx==1) lwork= 1+blocksize*(6 + 2*blocksize)
   lrwork= 1 + 5*blocksize + 2*(blocksize**2)
   liwork= 3 + 5*blocksize
   call xgv_gpu(1,'v','u',blocksize,grama,blocksize,gramb,blocksize,eigen,work(1:cplx*lwork),&
&               lwork,rwork(1:lrwork),lrwork,iwork(1:liwork),liwork,info)
   if (info/=0)  then
     write(msg,'(a,i3)')' Problem in xgv_gpu, info= ',info
     MSG_WARNING(msg)
   endif
 else
#endif
#if defined HAVE_LINALG_SCALAPACK
   if ( x_withsclapack ) then
     call compute_eigen2(x_communicator,x_processor,blocksize,&
&                        blocksize,grama,gramb,eigen,istwf_k)
   else
#endif
     lwork=3*blocksize - (cplx-1)*2
     call xgv(1,'v','u',blocksize,grama,blocksize,gramb,blocksize,eigen,work(1:cplx*lwork),&
&             lwork,rwork(1:lwork),info)
     if (info/=0)  then
       write(msg,'(a,i3)')' Problem in xgv, info= ',info
       MSG_WARNING(msg)
     endif
#if defined HAVE_LINALG_SCALAPACK
   endif
# endif
#if defined HAVE_LINALG_MAGMA
 endif
#endif

 if (present(tim_xeigen).and.present(timopt)) then
   if(abs(timopt)==3) call timab(tim_xeigen,2,tsec)
 end if

#if !defined HAVE_LINALG_SCALAPACK
 if (.false.) write(std_out,*) istwf_k
#endif

 end subroutine xeigen2


!!***

!! NAME
!! xev
!!
!! FUNCTION
!! Computes all eigenvalues and, optionally, eigenvectors
!! of a real symmetric matrix.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FBottin,CS,FDahm)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! jobz CHARACTER*1. Must be 'N' or 'V'.
!!     If jobz ='N', then only eigenvalues are computed.
!!     If jobz ='V', then eigenvalues and eigenvectors are computed.
!! uplo CHARACTER*1. Must be 'U' or 'L'.
!!     If uplo = 'U', a stores the upper triangular part of A.
!!     If uplo = 'L', a stores the lower triangular part of A.
!!     n INTEGER. The order of the matrix A (n \u2265 0).
!! a, work REAL for ssyev Arrays:
!!     a(lda,*) is an array containing either upper or lower triangular part of the
!!     symmetric matrix A, as specified by uplo.
!!     The second dimension of a must be at least max(1, n).
!!     work(lwork) is a workspace array.
!! lda INTEGER. The first dimension of the array a.
!!     Must be at least max(1, n).
!! lwork INTEGER. The dimension of the array work.
!!     Constraint: lwork \u2265 max(1, 3n-1).
!!     If lwork = -1, then a workspace query is assumed; the routine only calculates
!!     the optimal size of the work array, returns this value as the first entry of the
!!     work array, and no error message related to lwork is issued by xerbla.
!!
!! OUTPUT
!! a, w, work, info
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_eigen
!!
!! CHILDREN
!!      dcopy,dpotrf,dsygst,dtrmm,dtrsm,magmaf_dsyevd,magmaf_zheevd,xerbla
!!      zhegst,zpotrf,ztrmm,ztrsm
!!
!! SOURCE

subroutine xev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)

 use defs_basis
 use m_wfutils
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xev'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 real(dp), intent(inout) :: a(x_cplx*n,n)
 integer, intent(in) :: lda
 real(dp), intent(inout) :: work(x_cplx*lwork)
 integer, intent(in) :: lwork
 real(dp), intent(inout) :: rwork(lwork)
 real(dp), intent(out) :: w(n)
 integer, intent(out) :: info

! *********************************************************************

 if ( x_cplx == 1 ) then
   call dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
 else
   call zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
 endif
 return

end subroutine xev


!!***

!! NAME
!! xgv
!!
!! FUNCTION
!! Computes all eigenvalues and, optionally, eigenvectors
!! of a real generalized symmetric definite eigenproblem
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FBottin,CS,FDahm)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! itype INTEGER. Must be 1 or 2 or 3.
!!     Specifies the problem type to be solved:
!!     if itype = 1, the problem type is Ax = \u03bb Bx;
!!     if itype = 2, the problem type is ABx = \u03bb x;
!!     if itype = 3, the problem type is B Ax = \u03bb x.
!! jobz CHARACTER*1. Must be 'N' or 'V'.
!!     If jobz ='N', then only eigenvalues are computed.
!!     If jobz ='V', then eigenvalues and eigenvectors are computed.
!! uplo CHARACTER*1. Must be 'U' or 'L'.
!!     If uplo = 'U', arrays a and b store the upper triangles of A and B;
!!     If uplo = 'L', arrays a and b store the lower triangles of A and B.
!! n INTEGER. The order of the matrix A (n \u2265 0).
!!a, b, work
!!     a(lda,*) contains the upper or lower triangle of the symmetric matrix A, as
!!     specified by uplo.
!!     The second dimension of a must be at least max(1, n).
!!     b(ldb,*) contains the upper or lower triangle of the symmetric positive
!!     definite matrix B, as specified by uplo.
!!     The second dimension of b must be at least max(1, n).
!!     work(lwork) is a workspace array.
!!lda INTEGER. The first dimension of a; at least max(1, n).
!!ldb INTEGER. The first dimension of b; at least max(1, n).
!!lwork INTEGER. The dimension of the array work;
!!     lwork >= max(1, 3n-1).
!!     If lwork = -1, then a workspace query is assumed; the routine only calculates
!!     the optimal size of the work array, returns this value as the first entry of the
!!     work array, and no error message related to lwork is issued by xerbla.
!!
!! OUTPUT
!! a, b work, info
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_eigen
!!
!! CHILDREN
!!      dcopy,dpotrf,dsygst,dtrmm,dtrsm,magmaf_dsyevd,magmaf_zheevd,xerbla
!!      zhegst,zpotrf,ztrmm,ztrsm
!!
!! SOURCE

subroutine xgv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,rwork,info)

 use defs_basis
 use m_wfutils
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgv'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: itype
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 real(dp), intent(inout) :: a(x_cplx*n,n)
 integer, intent(in) :: lda
 real(dp), intent(inout) :: b(x_cplx*n,n)
 integer, intent(in) :: ldb
 real(dp), intent(inout) :: work(x_cplx*lwork)
 integer, intent(in) :: lwork
 real(dp), intent(inout) :: rwork(lwork)
 real(dp), intent(out) :: w(n)
 integer, intent(out) :: info

! *********************************************************************

 if ( x_cplx == 1 ) then
   call dsygv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,info)
 else
   call zhegv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,rwork,info)
 endif
 return

end subroutine xgv


!!***

#if defined HAVE_LINALG_MAGMA
!! NAME
!!  xev_gpu
!!
!! FUNCTION
!!   Computes all eigenvalues and, optionally, eigenvectors
!!   of a real symmetric matrix on gpu.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2011 abinit group (FDahm)
!!  this file is distributed under the terms of the
!!  gnu general public license, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  for the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! jobz character*1. must be 'n' or 'v'.
!!     if jobz ='n', then only eigenvalues are computed.
!!     if jobz ='v', then eigenvalues and eigenvectors are computed.
!! uplo character*1. must be 'u' or 'l'.
!!     if uplo = 'u', a stores the upper triangular part of a.
!!     if uplo = 'l', a stores the lower triangular part of a.
!!     n integer. the order of the matrix a (n \u2265 0).
!! a, work real for ssyev arrays:
!!     a(lda,*) is an array containing either upper or lower triangular part of the
!!     symmetric matrix a, as specified by uplo.
!!     the second dimension of a must be at least max(1, n).
!!     work(lwork) is a workspace array.
!! lda integer. the first dimension of the array a.
!!     must be at least max(1, n).
!! lwork integer. the dimension of the array work.
!!     constraint: lwork \u2265 max(1, 3n-1).
!!     if lwork = -1, then a workspace query is assumed; the routine only calculates
!!     the optimal size of the work array, returns this value as the first entry of the
!!     work array, and no error message related to lwork is issued by xerbla.
!!
!! OUTPUT
!! a, w, work, info
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_eigen
!!
!! CHILDREN
!!      dcopy,dpotrf,dsygst,dtrmm,dtrsm,magmaf_dsyevd,magmaf_zheevd,xerbla
!!      zhegst,zpotrf,ztrmm,ztrsm
!!
!! SOURCE

subroutine xev_gpu(jobz,uplo,n,a,lda,w,work,lwork,rwork,lrwork,iwork,liwork,info)

 use defs_basis
 use m_wfutils

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xev_gpu'
!End of the abilint section

 implicit none

!arguments ------------------------------------
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 real(dp), intent(inout) :: a(x_cplx*n,n)
 integer, intent(in) :: lda
 real(dp), intent(inout) :: work(x_cplx*lwork)
 integer, intent(in) :: lwork,lrwork,liwork
 real(dp), intent(inout) :: rwork(lrwork)
 integer, intent(inout) :: iwork(liwork)
 real(dp), intent(inout) :: w(n)
 integer, intent(inout) :: info

! *********************************************************************

 if ( x_cplx == 1 ) then
   !write(std_out,*) 'sizes:',n,lda,lwork,liwork
   !call magmaf_dsyevd(jobz,uplo,n,a,lda,w,work,-1,iwork,liwork,info)
   !write(std_out,*) 'optimal sizes:',work(1),iwork(1)
   call magmaf_dsyevd(jobz,uplo,n,a,lda,w,work,lwork,iwork,liwork,info)
 else
   !write(std_out,*) 'sizes:',n,lda,lwork,lrwork,liwork
   !call magmaf_zheevd(jobz,uplo,n,a,lda,w,work,-1,rwork,lrwork,iwork,liwork,info)
   !write(std_out,*) 'optimal sizes:',work(1),rwork(1),iwork(1)
   call magmaf_zheevd(jobz,uplo,n,a,lda,w,work,lwork,rwork,lrwork,iwork,liwork,info)
 endif
 end subroutine xev_gpu
!!***
#endif


#if defined HAVE_LINALG_MAGMA
!! NAME
!!   xgv_gpu
!!
!! FUNCTION
!!   Computes all eigenvalues and, optionally, eigenvectors
!!   of a real generalized symmetric definite eigenproblem using gpu
!!
!! COPYRIGHT
!! copyright (c) 1998-2011 abinit group (FDahm)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/copying
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! itype integer. must be 1 or 2 or 3.
!!     specifies the problem type to be solved:
!!     if itype = 1, the problem type is ax = \u03bb bx;
!!     if itype = 2, the problem type is abx = \u03bb x;
!!     if itype = 3, the problem type is b ax = \u03bb x.
!! jobz character*1. must be 'n' or 'v'.
!!     if jobz ='n', then only eigenvalues are computed.
!!     if jobz ='v', then eigenvalues and eigenvectors are computed.
!! uplo character*1. must be 'u' or 'l'.
!!     if uplo = 'u', arrays a and b store the upper triangles of a and b;
!!     if uplo = 'l', arrays a and b store the lower triangles of a and b.
!! n integer. the order of the matrix a (n \u2265 0).
!! a, b, work
!!     a(lda,*) contains the upper or lower triangle of the symmetric matrix a, as
!!     specified by uplo.
!!     the second dimension of a must be at least max(1, n).
!!     b(ldb,*) contains the upper or lower triangle of the symmetric positive
!!     definite matrix b, as specified by uplo.
!!     the second dimension of b must be at least max(1, n).
!!     work(lwork) is a workspace array.
!! lda integer. the first dimension of a; at least max(1, n).
!! ldb integer. the first dimension of b; at least max(1, n).
!! lwork integer. the dimension of the array work;
!!     lwork >= max(1, 3n-1).
!!     if lwork = -1, then a workspace query is assumed; the routine only calculates
!!     the optimal size of the work array, returns this value as the first entry of the
!!     work array, and no error message related to lwork is issued by xerbla.
!!
!! OUTPUT
!! a, b work, info
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_eigen
!!
!! CHILDREN
!!      dcopy,dpotrf,dsygst,dtrmm,dtrsm,magmaf_dsyevd,magmaf_zheevd,xerbla
!!      zhegst,zpotrf,ztrmm,ztrsm
!!
!! SOURCE

 subroutine xgv_gpu(itype,jobz,uplo,n,a,lda,b,ldb,w,work,&
&                   lwork,rwork,lrwork,iwork,liwork,info)

 use defs_basis
 use m_wfutils

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgv_gpu'
!End of the abilint section

 implicit none

!arguments ------------------------------------
 integer, intent(in) :: itype
 character(len=1), intent(in) :: jobz
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n
 real(dp), intent(inout) :: a(x_cplx*n,n)
 integer, intent(in) :: lda
 real(dp), intent(inout) :: b(x_cplx*n,n)
 integer, intent(in) :: ldb
 real(dp), intent(inout) :: work(x_cplx*lwork)
 integer, intent(in) :: lwork,lrwork,liwork
 real(dp), intent(inout) :: rwork(lwork)
 integer, intent(inout) :: iwork(liwork)
 real(dp), intent(inout) :: w(n)
 integer, intent(inout) :: info
!local variables  ----------------------------------
!scalars
 logical :: lquery, upper, wantz
 character :: trans
 integer :: liopt, liwmin, lopt, lwmin,lropt, lrwmin
!arrays
 complex(dpc), dimension(:,:),allocatable :: z_a,z_b
!external functions ------------------------
 logical :: lsame
 external :: lsame

! *********************************************************************

! !IMPORTANT!
! WHEN MAGMA GVD FUNCTIONS WILL BE COMPLETLY FUNCTIONAL
! THE BODY OF THIS FUNCTION SHOULDE BE :
!     if ( x_cplx == 1 ) then
!        call magmaf_dsygvd(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,iwork,liwork,info)
!     else
!        call magmaf_zhegvd(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,rwork,lrwork,iwork,liwork,info)
!     endif
! NOWADAYS (SEPTEMBER 2011) MAGMA IS IN V1.0.0 AND MAGMAGF_ZHEEVD FAILS FOR SOME SMALL CASE SO WE USE
! THIS HAND MADE FUNCTION INSTEAD, WHICH IS NOT OPTIMAL
 if ( x_cplx == 1 ) then
   !DSYGVD BOBY
   wantz = lsame( jobz, 'v' )
   upper = lsame( uplo, 'u' )
   lquery = ( lwork.eq.-1 .or. liwork.eq.-1 )
   !
   info = 0
   if( n.le.1 ) then
     liwmin = 1
     lwmin = 1
   else if( wantz ) then
     liwmin = 3 + 5*n
     lwmin = 1 + 6*n + 2*n**2
   else
     liwmin = 1
     lwmin = 2*n + 1
   end if
   lopt = lwmin
   liopt = liwmin
   if( itype.lt.1 .or. itype.gt.3 ) then
     info = -1
   else if( .not.( wantz .or. lsame( jobz, 'n' ) ) ) then
     info = -2
   else if( .not.( upper .or. lsame( uplo, 'l' ) ) ) then
     info = -3
   else if( n.lt.0 ) then
     info = -4
   else if( lda.lt.max( 1, n ) ) then
     info = -6
   else if( ldb.lt.max( 1, n ) ) then
     info = -8
   end if

   if( info.eq.0 ) then
     work( 1 ) = lopt
     iwork( 1 ) = liopt
     if( lwork.lt.lwmin .and. .not.lquery ) then
       info = -11
     else if( liwork.lt.liwmin .and. .not.lquery ) then
       info = -13
     end if
   end if

   if( info.ne.0 ) then
     call xerbla( 'dsygvd', -info )
     return
   else if( lquery ) then
     return
   end if
   !     quick return if possible
   if( n.eq.0 )  return
   !*     form a cholesky factorization of b.
   call dpotrf( uplo, n, b, ldb, info )
   if( info.ne.0 ) then
     info = n + info
     return
   end if
!  Transform problem to standard eigenvalue problem and solve.
   call dsygst( itype, uplo, n, a, lda, b, ldb, info )
   call magmaf_dsyevd( jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork,info)
   lopt = max( dble( lopt ), dble( work( 1 ) ) )
   liopt = max( dble( liopt ), dble( iwork( 1 ) ) )
   if( wantz .and. info.eq.0 ) then
     !*        backtransform eigenvectors to the original problem.
     if( itype.eq.1 .or. itype.eq.2 ) then
       !*           for a*x=(lambda)*b*x and a*b*x=(lambda)*x;
       !*           backtransform eigenvectors: x = inv(l)**t*y or inv(u)*y
       if( upper ) then
         trans = 'n'
       else
         trans = 't'
       end if
       call dtrsm( 'left', uplo, trans, 'non-unit', n, n, one,b, ldb, a, lda )
     else if( itype.eq.3 ) then
       !           for b*a*x=(lambda)*x;
       !           backtransform eigenvectors: x = l*y or u**t*y
       if( upper ) then
         trans = 't'
       else
         trans = 'n'
       end if
       call dtrmm( 'left', uplo, trans, 'non-unit', n, n, one,b, ldb, a, lda )
     end if
   end if
   work( 1 ) = lopt
   iwork( 1 ) = liopt
 else
   !ZHEGVD BODY
   wantz = lsame( jobz, 'v' )
   upper = lsame( uplo, 'u' )
   lquery = ( lwork.eq.-1 .or. lrwork.eq.-1 .or. liwork.eq.-1 )
   info = 0
   if( n.le.1 ) then
     lwmin = 1
     lrwmin = 1
     liwmin = 1
   else if( wantz ) then
     lwmin = 2*n + n*n
     lrwmin = 1 + 5*n + 2*n*n
     liwmin = 3 + 5*n
   else
     lwmin = n + 1
     lrwmin = n
     liwmin = 1
   end if
   lopt = lwmin
   lropt = lrwmin
   liopt = liwmin
   if( itype.lt.1 .or. itype.gt.3 ) then
     info = -1
   else if( .not.( wantz .or. lsame( jobz, 'n' ) ) ) then
     info = -2
   else if( .not.( upper .or. lsame( uplo, 'l' ) ) ) then
     info = -3
   else if( n.lt.0 ) then
     info = -4
   else if( lda.lt.max( 1, n ) ) then
     info = -6
   else if( ldb.lt.max( 1, n ) ) then
     info = -8
   end if

   if( info.eq.0 ) then
     work( 1 ) = lopt
     rwork( 1 ) = lropt
     iwork( 1 ) = liopt

     if( lwork.lt.lwmin .and. .not.lquery ) then
       info = -11
     else if( lrwork.lt.lrwmin .and. .not.lquery ) then
       info = -13
     else if( liwork.lt.liwmin .and. .not.lquery ) then
       info = -15
     end if
   end if

   if( info.ne.0 ) then
     call xerbla( 'zhegvd', -info )
     return
   else if( lquery ) then
     return
   end if
!*
!* Quick return if possible
!*
   if( n.eq.0 )  return
!*
!* Form a cholesky factorization of b.
!*
   call zpotrf( uplo, n, b, ldb, info )
   if( info.ne.0 ) then
     info = n + info
     return
   end if
!*
!* Transform problem to standard eigenvalue problem and solve.
!*
   call zhegst( itype, uplo, n, a, lda, b, ldb, info )
   call magmaf_zheevd( jobz, uplo, n, a, lda, w, work, lwork, rwork, &
&                      lrwork,iwork, liwork, info )
   lopt = max( dble( lopt ), dble( work( 1 ) ) )
   lropt = max( dble( lropt ), dble( rwork( 1 ) ) )
   liopt = max( dble( liopt ), dble( iwork( 1 ) ) )
   if( wantz .and. info.eq.0 ) then
     !*        backtransform eigenvectors to the original problem.
     if( itype.eq.1 .or. itype.eq.2 ) then
       !*           for a*x=(lambda)*b*x and a*b*x=(lambda)*x;
       !*           backtransform eigenvectors: x = inv(l)**h *y or inv(u)*y
       if( upper ) then
         trans = 'n'
       else
         trans = 'c'
       end if
       call ztrsm( 'left', uplo, trans, 'non-unit', n, n, cone,b, ldb, a, lda )
     else if( itype.eq.3 ) then
       !*           for b*a*x=(lambda)*x;
       !*           backtransform eigenvectors: x = l*y or u**h *y
       if( upper ) then
         trans = 'c'
       else
         trans = 'n'
       end if
       ABI_ALLOCATE(z_a,(n,n))
       ABI_ALLOCATE(z_b,(n,n))
       call dcopy(x_cplx*n*n,a,1,z_a)
       call dcopy(x_cplx*n*n,b,1,z_b)
       call ztrmm( 'left', uplo, trans, 'non-unit', n, n, cone, z_b, ldb, z_a, lda )
       call dcopy(x_cplx*n*n,z_a,1,a)
       call dcopy(x_cplx*n*n,z_b,1,b)
       ABI_DEALLOCATE(z_a)
       ABI_DEALLOCATE(z_b)
     end if
   end if
   work( 1 ) = lopt
   rwork( 1 ) = lropt
   iwork( 1 ) = liopt
   !*     end of zhegvd
 endif

 end subroutine xgv_gpu
!!***
#endif

end module m_eigen
!!***
