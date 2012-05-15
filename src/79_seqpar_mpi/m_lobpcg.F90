!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_lobpcg
!! NAME
!!  m_lobpcg
!!
!! FUNCTION
!!  This module provides the procedures used in the LOBPCGWF routine.
!!  They permit to hide the complex/real form of the WFs.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (FBottin,CS,FDahm,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_lobpcg

 use m_profiling

 use defs_basis
 use defs_datatypes
 use m_errors

#if defined HAVE_FC_ISO_C_BINDING
 use iso_c_binding, only: c_ptr
#endif

 implicit none

 private

#ifndef HAVE_FC_ISO_C_BINDING
 type, public :: c_ptr
   integer :: val
 end type
#endif

!public procedures.
 public :: my_xcopy
 public :: my_xgemm
 public :: my_xtrsm
 public :: xorthonormalize
 public :: gpu_xorthonormalize
 public :: xprecon
#ifndef HAVE_GPU_CUDA
!dummy routines replace gpu helper routines
 public :: alloc_on_gpu
 public :: copy_from_gpu
 public :: copy_on_gpu
 public :: dealloc_on_gpu
 public :: gpu_linalg_init
 public :: gpu_linalg_shutdown
 public :: gpu_xgemm
 public :: gpu_xtrsm
#endif
!!***

CONTAINS

!----------------------------------------------------------------------

!!****f* m_lobpcg/my_xcopy
!! NAME
!!  my_xcopy
!!
!! FUNCTION
!!  The xcopy routines perform a vector-vector operation defined as y = x
!!  where x and y are vectors.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FBottin,CS)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  size= Specifies the order of vectors x and y.
!!  tsrc= for dcopy or zcopy Array, DIMENSION at least (1 + (n-1)*abs(incx)).
!!  incsrc= Specifies the increment for the elements of x.
!!  tdest= for dcopy or zcopy Array, DIMENSION at least (1 + (n-1)*abs(incy)).
!!  incdest= Specifies the increment for the elements of y.
!!
!! OUTPUT
!!  tdest
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      lobpcgwf
!!
!! CHILDREN
!!
!! SOURCE

subroutine my_xcopy(size,tsrc,incsrc,tdest,incdest,&
&                timopt,tim_xcopy) ! optional arguments
 use defs_basis
 use m_wfutils

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'my_xcopy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: size
 integer, intent(in) :: incsrc
 integer, intent(in) :: incdest
 integer, intent(in), optional :: timopt,tim_xcopy
 real(dp), dimension(*) :: tsrc
 real(dp), dimension(*) :: tdest
!Local variables-------------------------------
 real(dp) :: tsec(2)

! *********************************************************************

 if (present(tim_xcopy).and.present(timopt)) then
   if(abs(timopt)==3) call timab(tim_xcopy,1,tsec)
 end if
 if ( x_cplx == 1 ) then
   call dcopy(size,tsrc,incsrc,tdest,incdest)
 else
   call zcopy(size,tsrc,incsrc,tdest,incdest)
 endif
 if (present(tim_xcopy).and.present(timopt)) then
   if(abs(timopt)==3) call timab(tim_xcopy,2,tsec)
 end if
 return

end subroutine my_xcopy
!!***

!----------------------------------------------------------------------

!!****f* m_lobpcg/my_xgemm
!! NAME
!!  my_xgemm
!!
!! FUNCTION
!!  Compute a scalar-matrix-matrix product and return a scalar-matrix product
!!  c = alpha * op(a) * op(b) + beta * c
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FBottin,CS)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  transa= from of op(a) to be used in the matrix multiplication
!!  transb= from of op(b) to be used in the matrix multiplication
!!  m     = number of rows of the matrix op(a) and of the matrix c
!!  n     = number of rows of the matrix op(b) and the number of columns of the matrix c
!!  k     = number of columns of the matrix op(a) and the number of rows of the matrix op(b)
!!  alpha = alpha scalar coefficient for matrix op(a)
!!  a     = a matrix
!!  lda   = first dimension of a
!!  b     = b matrix
!!  ldb   = first dimension of b
!!  beta  = beta scalar coefficient for matrix c
!!  c     = c matrix
!!  ldc   = first dimension of c
!!
!! OUTPUT
!!  c     = c matrix
!!
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      lobpcgwf
!!
!! CHILDREN
!!
!! SOURCE
subroutine my_xgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc,&
&                timopt,tim_xgemm) ! optional arguments
 use defs_basis
 use m_wfutils
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'my_xgemm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: lda
 integer, intent(in) :: ldb
 integer, intent(in) :: ldc
 integer, intent(in) :: m
 integer, intent(in) :: n
 integer, intent(in) :: k
 integer, intent(in), optional :: timopt,tim_xgemm
 real(dp), dimension(*), intent(in) :: a
 real(dp), dimension(*), intent(in) :: b
 real(dp), dimension(*), intent(inout) :: c
 complex(dpc), intent(in) :: alpha
 complex(dpc), intent(in) :: beta
 character(len=1), intent(in) :: transa
 character(len=1), intent(in) :: transb
!Local variables-------------------------------
 real(dp) :: tsec(2)

! *********************************************************************

 if (present(tim_xgemm).and.present(timopt)) then
   if(abs(timopt)==3) call timab(tim_xgemm,1,tsec)
 end if
 if ( x_cplx == 1 ) then
   call dgemm(transa,transb,m,n,k,real(alpha,dp),a,lda,b,ldb,real(beta,dp),c,ldc)
 else
   call zgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
 endif
 if (present(tim_xgemm).and.present(timopt)) then
   if(abs(timopt)==3) call timab(tim_xgemm,2,tsec)
 end if
 return

end subroutine my_xgemm
!!***

!----------------------------------------------------------------------

!!****f* m_lobpcg/my_xtrsm
!! NAME
!!  my_xtrsm
!!
!! FUNCTION
!! Solves a matrix equation (one matrix operand is triangular).
!! The xtrsm routines solve one of the following matrix equations
!! op(a)*x = alpha*b
!! or
!! x*op(a) = alpha*b,
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FBottin,CS)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! side= Specifies whether op(a) appears on the left or right of x for
!!      the operation to be performed as follows:
!!      L or l op(a)*x = alpha*b
!!      R or r x*op(a) = alpha*b
!! uplo= Specifies whether the matrix a is an upper or lower triangular
!!      matrix as follows:
!!      U or u Matrix a is an upper triangular matrix.
!!      L or l Matrix a is a lower triangular matrix
!! transa= Specifies the form of op(a) to be used in the matrix
!!      multiplication as follows:
!!      N or n op(a) = a
!!      T or t op(a) = a'
!!      C or c op(a) = conjg(a')
!! diag= Specifies whether or not a is unit triangular as follows:
!!      U or u Matrix a is assumed to be unit triangular.
!!      N or n Matrix a is not assumed to be unit triangular.
!! m= Specifies the number of rows of b. The value of m must be at least zero
!! n= Specifies the number of columns of b. The value of n must be at least zero
!! alpha= Specifies the scalar alpha. When alpha is zero, then a is not referenced and b
!!      need not be set before entry.
!! a= Array, DIMENSION (lda, k), where k is m when side = 'L' or 'l' and is n
!!      when side = 'R' or 'r'.
!! lda= Specifies the first dimension of a as declared in the calling
!!     (sub)program. When side = 'L' or 'l', then lda must be at least max(1,
!!      m), when side = 'R' or 'r', then lda must be at least max(1, n).
!! b= Array, DIMENSION (ldb,n). Before entry, the leading m-by-n part of the array
!!     b must contain the right-hand side matrix b.
!! ldb= Specifies the first dimension of b as declared in the calling
!!     (sub)program. The value of ldb must be at least max(1, m).
!!
!! OUTPUT
!!  b
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      lobpcgwf
!!
!! CHILDREN
!!
!! SOURCE
subroutine my_xtrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb,&
&                timopt,tim_xtrsm) ! optional arguments
 use defs_basis
 use m_wfutils
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'my_xtrsm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=1), intent(in) :: side
 character(len=1), intent(in) :: uplo
 character(len=1), intent(in) :: transa
 character(len=1), intent(in) :: diag
 integer, intent(in) :: m
 integer, intent(in) :: n
 complex(dpc), intent(in) :: alpha
 real(dp), dimension(*), intent(in) :: a
 integer, intent(in) :: lda
 real(dp), dimension(*), intent(inout) :: b
 integer, intent(in) :: ldb
 integer, intent(in), optional :: timopt,tim_xtrsm

!Local variables-------------------------------
 real(dp) :: tsec(2)

! *********************************************************************

 if (present(tim_xtrsm).and.present(timopt)) then
   if(abs(timopt)==3) call timab(tim_xtrsm,1,tsec)
 end if
 if ( x_cplx == 1 ) then
   call dtrsm(side,uplo,transa,diag,m,n,real(alpha,dp),a,lda,b,ldb)
 else
   call ztrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
 endif
 if (present(tim_xtrsm).and.present(timopt)) then
   if(abs(timopt)==3) call timab(tim_xtrsm,2,tsec)
 end if
 return

end subroutine my_xtrsm
!!***

!----------------------------------------------------------------------

!! NAME
!!  xorthonormalize
!!
!! FUNCTION
!!  This routine computes the overlap of two complex wavefunctions (for a given number of bands)
!!  and orthonormalizes it:
!!      - Computes the products of two rectangular matrices
!!         containing the wavefunctions psi and S.psi (where S is the
!!         overlap (with the PAW terms if necessary)).
!!      - Does a Cholesky decomposition of this overlap
!!      - rotates the initial matrix blockvectorx by the triangular matrix to
!!         have an orthonormal set of wavefunctions
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FBottin,CS)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  blockvectorbx = matrix of dimension (blocksize,vectsize)
!!                  (e.g. block of overlap*wavefunction)
!!  blocksize     = dimension of matrices (e.g number of bands)
!!  mpi_enreg     = informations about MPI parallelization
!!  vectsize      = dimension of matrices (e.g number of G vector)
!!
!! OUTPUT
!!  sqgram        = Choleski decomposition of transpose(blockvector)*blockvectorx
!!
!! SIDE EFFECTS
!!  blockvectorx  = on input, matrix of dimension (vectsize,blocksize)
!!                  (e.g block of wavefunction)
!!  blockvectorx  = on output, orthonormalized wavefunction.
!!
!!
!! PARENTS
!!      lobpcgwf
!!
!! CHILDREN
!!
!! SOURCE

subroutine xorthonormalize(blockvectorx,blockvectorbx,blocksize,mpi_enreg,sqgram,vectsize,&
&                          timopt,tim_xortho) ! optional arguments

 use defs_basis
 use defs_abitypes
 use m_wfutils

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xorthonormalize'
 use interfaces_53_abiutil
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: blocksize,vectsize
 integer, intent(in), optional :: timopt,tim_xortho
 type(mpi_type) :: mpi_enreg
!arrays
 real(dp) :: blockvectorbx(vectsize,blocksize),blockvectorx(vectsize,blocksize)
 real(dp) :: sqgram(blocksize,blocksize)

!Local variables-------------------------------
 complex(dpc),dimension(:,:), allocatable :: z_blockvectorbx,z_blockvectorx,z_sqgram
 real(dp) :: tsec(2)

! *********************************************************************

 if (present(tim_xortho).and.present(timopt)) then
   if(abs(timopt)==3) call timab(tim_xortho,1,tsec)
 end if
 if ( x_cplx == 1 ) then
   call orthonormalize(blockvectorx,blockvectorbx,blocksize,mpi_enreg,sqgram,vectsize)
 else
   ABI_ALLOCATE(z_blockvectorbx,(vectsize,blocksize))
   ABI_ALLOCATE(z_blockvectorx,(vectsize,blocksize))
   ABI_ALLOCATE(z_sqgram,(blocksize,blocksize))

   call dcopy(x_cplx*vectsize*blocksize,blockvectorbx,1,z_blockvectorbx,1)
   call dcopy(x_cplx*vectsize*blocksize,blockvectorx,1,z_blockvectorx,1)
   call dcopy(x_cplx*blocksize*blocksize,sqgram,1,z_sqgram,1)

   call zorthonormalize(z_blockvectorx,z_blockvectorbx,blocksize,mpi_enreg,z_sqgram,vectsize)

   call dcopy(x_cplx*vectsize*blocksize,z_blockvectorbx,1,blockvectorbx,1)
   call dcopy(x_cplx*vectsize*blocksize,z_blockvectorx,1,blockvectorx,1)
   call dcopy(x_cplx*blocksize*blocksize,z_sqgram,1,sqgram,1)

   ABI_DEALLOCATE(z_blockvectorbx)
   ABI_DEALLOCATE(z_blockvectorx)
   ABI_DEALLOCATE(z_sqgram)

 endif
 if (present(tim_xortho).and.present(timopt)) then
   if(abs(timopt)==3) call timab(tim_xortho,2,tsec)
 end if
 return

end subroutine xorthonormalize
!!***
!----------------------------------------------------------------------

!! NAME
!!  xprecon
!!
!! FUNCTION
!!  precondition $<G|(H-e_{n,k})|C_{n,k}>$
!!  for a block of band (band-FFT parallelisation)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FBottin,CS)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  blocksize= size of blocks of bands
!!  $cg(vectsize,blocksize)=<G|C_{n,k}> for a block of bands$.
!!  $eval(blocksize,blocksize)=current block of bands eigenvalues=<C_{n,k}|H|C_{n,k}>$.
!!  $ghc(vectsize,blocksize)=<G|H|C_{n,k}> for a block of bands$.
!!  iterationnumber=number of iterative minimizations in LOBPCG
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (Hartree)
!!  mpi_enreg=informations about MPI parallelization
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  $vect(vectsize,blocksize)=<G|H|C_{n,k}> for a block of bands$.
!!  npw=number of planewaves at this k point.
!!  optekin= 1 if the kinetic energy used in preconditionning is modified
!!             according to Kresse, Furthmuller, PRB 54, 11169 (1996)
!!           0 otherwise
!!  optpcon= 0 the TPA preconditionning matrix does not depend on band
!!           1 the TPA preconditionning matrix (not modified)
!!           2 the TPA preconditionning matrix is independant of iterationnumber
!!  vectsize= size of vectors
!!
!! OUTPUT
!!  vect(2,npw)=<g|(h-eval)|c_{n,k}>*(polynomial ratio)
!!
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      lobpcgwf
!!
!! CHILDREN
!!
!! SOURCE

subroutine xprecon(cg,eval,blocksize,iterationnumber,kinpw,&
& mpi_enreg,npw,nspinor,optekin,optpcon,pcon,ghc,vect,vectsize,&
&                timopt,tim_xprecon) ! optional arguments
 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_wfutils
!Arguments ------------------------------------
!scalars

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xprecon'
 use interfaces_66_wfs
!End of the abilint section

 integer,intent(in) :: blocksize,iterationnumber,npw,nspinor,optekin
 integer,intent(in) :: optpcon,vectsize
 integer, intent(in), optional :: timopt,tim_xprecon
 type(mpi_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: cg(vectsize,blocksize),eval(blocksize,blocksize)
 real(dp),intent(in) :: ghc(vectsize,blocksize),kinpw(npw)
 real(dp),intent(inout) :: pcon(npw,blocksize),vect(vectsize,blocksize)

!Local variables-------------------------------
 complex(dpc),dimension(:,:),allocatable :: z_cg,z_eval,z_ghc,z_vect
 real(dp) :: tsec(2)

! *********************************************************************

 if (present(tim_xprecon).and.present(timopt)) then
   if(abs(timopt)==3) call timab(tim_xprecon,1,tsec)
 end if
 if ( x_cplx == 1 ) then
  call precon2(cg,eval,blocksize,iterationnumber,kinpw,&
& mpi_enreg,npw,nspinor,optekin,optpcon,pcon,ghc,vect,vectsize)
 else
  ABI_ALLOCATE(z_cg,(vectsize,blocksize))
  ABI_ALLOCATE(z_eval,(blocksize,blocksize))
  ABI_ALLOCATE(z_ghc,(vectsize,blocksize))
  ABI_ALLOCATE(z_vect,(vectsize,blocksize))

  call dcopy(x_cplx*vectsize*blocksize,cg,1,z_cg,1)
  call dcopy(x_cplx*vectsize*blocksize,ghc,1,z_ghc,1)
  call dcopy(x_cplx*vectsize*blocksize,vect,1,z_vect,1)
  call dcopy(x_cplx*blocksize*blocksize,eval,1,z_eval,1)

  call zprecon3(z_cg,z_eval,blocksize,iterationnumber,kinpw,&
& mpi_enreg,npw,nspinor,optekin,optpcon,pcon,z_ghc,z_vect,vectsize)

  call dcopy(x_cplx*vectsize*blocksize,z_cg,1,cg,1)
  call dcopy(x_cplx*vectsize*blocksize,z_ghc,1,ghc,1)
  call dcopy(x_cplx*vectsize*blocksize,z_vect,1,vect,1)
  call dcopy(x_cplx*blocksize*blocksize,z_eval,1,eval,1)

  ABI_DEALLOCATE(z_cg)
  ABI_DEALLOCATE(z_eval)
  ABI_DEALLOCATE(z_ghc)
  ABI_DEALLOCATE(z_vect)
 endif
 if (present(tim_xprecon).and.present(timopt)) then
   if(abs(timopt)==3) call timab(tim_xprecon,2,tsec)
 end if
 return

end subroutine xprecon
!!***

!----------------------------------------------------------------------

!!****f* m_lobpcg/gpu_xorthonormalize
!! NAME
!!  gpu_xorthonormalize
!!
!! FUNCTION
!!  This routine computes the overlap of two complex wavefunctions (for a given number of bands)
!!  and orthonormalizes it using gpu:
!!      - Computes the products of two rectangular matrices
!!         containing the wavefunctions psi and S.psi (where S is the
!!         overlap (with the PAW terms if necessary)).
!!      - Does a Cholesky decomposition of this overlap
!!      - rotates the initial matrix blockvectorx by the triangular matrix to
!!         have an orthonormal set of wavefunctions
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (FDahm)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  blockvectorbx = matrix of dimension (blocksize,vectsize) as a GPU ptr
!!                  (e.g. block of overlap*wavefunction)
!!  blocksize     = dimension of matrices (e.g number of bands)
!!  mpi_enreg     = informations about MPI parallelization
!!  vectsize      = dimension of matrices (e.g number of G vector)
!!
!! OUTPUT
!!  sqgram        = Choleski decomposition of transpose(blockvector)*blockvectorx as a GPU ptr
!!
!! SIDE EFFECTS
!!  blockvectorx  = on input, matrix of dimension (vectsize,blocksize) as a GPU ptr
!!                  (e.g block of wavefunction)
!!  blockvectorx  = on output, orthonormalized wavefunction. as a GPU ptr
!!
!!
!! PARENTS
!!      lobpcgwf
!!
!! CHILDREN
!!
!! SOURCE

subroutine gpu_xorthonormalize(blockvectorx_gpu,blockvectorbx_gpu,blocksize,mpi_enreg,&
&                              sqgram_gpu,vectsize,&
     &                         timopt,tim_xortho) ! optional arguments

 use defs_basis
 use defs_abitypes
 use m_wfutils
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gpu_xorthonormalize'
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: blocksize,vectsize
 integer, intent(in), optional :: timopt,tim_xortho
 type(mpi_type) :: mpi_enreg
!arrays
 type(c_ptr),intent(inout) :: blockvectorbx_gpu, blockvectorx_gpu,sqgram_gpu
!Local variables-------------------------------
#if defined HAVE_GPU_CUDA
 integer :: ierr,info,old_paral_level,spaceComm
 real(dp),dimension(:,:),allocatable :: d_sqgram
 complex(dpc),dimension(:,:),allocatable :: z_sqgram
 character :: tr
 real(dp) :: tsec(2)
#else
 type(c_ptr) :: cptr_a
#endif
 character(len=500) :: message

! *********************************************************************
#if defined HAVE_GPU_CUDA
 if (present(tim_xortho).and.present(timopt)) then
   if(abs(timopt)==3) call timab(tim_xortho,1,tsec)
 end if

 if ( x_cplx == 1 ) then
   tr='t'
   ABI_ALLOCATE(d_sqgram,(blocksize,blocksize))
 else
   tr='c'
   ABI_ALLOCATE(z_sqgram,(blocksize,blocksize))
 end if

 call gpu_xgemm(x_cplx,tr,'n',blocksize,blocksize,vectsize, &
& cone,blockvectorx_gpu,vectsize,blockvectorbx_gpu,vectsize,czero,sqgram_gpu,blocksize)

 old_paral_level= mpi_enreg%paral_level
 mpi_enreg%paral_level=3

 call xcomm_init(mpi_enreg,spaceComm,spaceComm_bandfft=mpi_enreg%commcart)
 if ( x_cplx == 1 ) then
   call copy_from_gpu(d_sqgram,sqgram_gpu,x_cplx*dp*blocksize*blocksize)
   call xsum_mpi(d_sqgram,spaceComm,ierr)
   call dpotrf('u',blocksize,d_sqgram,blocksize,info)
   call copy_on_gpu(d_sqgram,sqgram_gpu,x_cplx*dp*blocksize*blocksize)
 else
   call copy_from_gpu(z_sqgram,sqgram_gpu,x_cplx*dp*blocksize*blocksize)
   call xsum_mpi(z_sqgram,spaceComm,ierr)
   call zpotrf('u',blocksize,z_sqgram,blocksize,info)
   call copy_on_gpu(z_sqgram,sqgram_gpu,x_cplx*dp*blocksize*blocksize)
 end if
 mpi_enreg%paral_level= old_paral_level

 if (info /= 0 )  then
   write(message,'(a,i3)') '  xpotrf, info=',info
   MSG_WARNING(message)
 end if

 call gpu_xtrsm(x_cplx,'r','u','n','n',vectsize,blocksize,cone,sqgram_gpu,blocksize,&
&               blockvectorx_gpu,vectsize)

 if(x_cplx==1) then
   ABI_DEALLOCATE(d_sqgram)
 else
   ABI_DEALLOCATE(z_sqgram)
 end if
 if (present(tim_xortho).and.present(timopt)) then
   if(abs(timopt)==3) call timab(tim_xortho,2,tsec)
 end if
 return

#else
 message='  This routine is not allowed when Cuda is disabled !'
 MSG_BUG(message)
 if (.false.) then
   write(std_out,*) mpi_enreg%nproc,blocksize,vectsize
   if (present(timopt))  write(std_out,*) timopt
   if (present(tim_xortho))  write(std_out,*) tim_xortho
   cptr_a=blockvectorbx_gpu;cptr_a=blockvectorx_gpu;cptr_a=sqgram_gpu
 end if
#endif

end subroutine gpu_xorthonormalize
!!***

!---- Fortran interfaces to GPU helper routines  ----
#ifndef HAVE_GPU_CUDA

!!****f* m_lobpcg/alloc_on_gpu
!! NAME
!!  alloc_on_gpu
!!
!! FUNCTION
!!  Allocate size byte in gpu memory and returns in gpu_ptr this location
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (FDahm)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  size= size in byte to allocate
!!
!! OUTPUT
!!  gpu_ptr= C_PTR on gpu memory location that has been allocated
!!
!! SIDE EFFECTS
!!   WARNING! : this routine is a dummy one when HAVE_GPU_CUDA is not enabled
!!   the correct one is in 15_gpu_toolbox/dev_spec.cu
!!
!! PARENTS
!!      lobpcgwf
!!
!! CHILDREN
!!
!! SOURCE

subroutine alloc_on_gpu(gpu_ptr,size)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'alloc_on_gpu'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer :: size ! size in bytes to allocate
 type(c_ptr) :: gpu_ptr
!Local variables ------------------------------
 type(c_ptr) :: cptr

! *********************************************************************

 if (.false.) then
   cptr=gpu_ptr;write(std_out,*) size
 end if

end subroutine alloc_on_gpu
!!***

!!****f* m_lobpcg/copy_from_gpu
!! NAME
!!  copy_from_gpu
!!
!! FUNCTION
!!  copy size byte from gpu memory pointed by gpu_ptr to dtab
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (FDahm)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  size= size in byte to allocate
!!  gpu_ptr= C_PTR on gpu memory location that has been allocated
!!
!! OUTPUT
!!  dtab = fortran tab which will contains data
!!
!! SIDE EFFECTS
!!   WARNING! : this routine is a dummy one when HAVE_GPU_CUDA is not enabled
!!   the correct one is in 15_gpu_toolbox/dev_spec.cu
!!
!! PARENTS
!!      lobpcgwf,m_lobpcg
!!
!! CHILDREN
!!
!! SOURCE

subroutine copy_from_gpu(dtab,gpu_ptr,size)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'copy_from_gpu'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer :: size ! taille en octet a transferer
 real(dp),dimension(*),optional :: dtab
 type(c_ptr) :: gpu_ptr
!Local variables ------------------------------
 type(c_ptr) :: cptr

! *********************************************************************

 if (.false.) then
   cptr=gpu_ptr;write(std_out,*) size,dtab(1)
 end if

end subroutine copy_from_gpu
!!***

!!****f* m_lobpcg/copy_on_gpu
!! NAME
!!  copy_on_gpu
!!
!! FUNCTION
!!  copy size byte from  dtab to gpu memory pointed by gpu_ptr
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (FDahm)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  size= size in byte to allocate
!!  dtab = fortran tab to copy
!!
!! OUTPUT
!!  gpu_ptr= C_PTR on gpu memory location
!!
!! SIDE EFFECTS
!!   WARNING! : this routine is a dummy one when HAVE_GPU_CUDA is not enabled
!!   the correct one is in 15_gpu_toolbox/dev_spec.cu
!!
!! PARENTS
!!      lobpcgwf,m_lobpcg
!!
!! CHILDREN
!!
!! SOURCE

subroutine copy_on_gpu(dtab,gpu_ptr,size)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'copy_on_gpu'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer :: size ! size in byte (to be transfered)
 real(dp), dimension(*),optional :: dtab
 type(c_ptr) :: gpu_ptr
!Local variables ------------------------------
 type(c_ptr) :: cptr

! *********************************************************************

 if (.false.) then
   cptr=gpu_ptr;write(std_out,*) size,dtab(1)
 end if

end subroutine copy_on_gpu
!!***

!!****f* m_lobpcg/dealloc_on_gpu
!! NAME
!!  dealloc_on_gpu
!!
!! FUNCTION
!!  free memory location pointed by gpu_ptr
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (FDahm)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!  gpu_ptr= C_PTR on gpu memory location that has been allocated
!!
!! SIDE EFFECTS
!!   WARNING! : this routine is a dummy one when HAVE_GPU_CUDA is not enabled
!!   the correct one is in 15_gpu_toolbox/dev_spec.cu
!!
!! PARENTS
!!      lobpcgwf
!!
!! CHILDREN
!!
!! SOURCE

subroutine dealloc_on_gpu(gpu_ptr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dealloc_on_gpu'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(c_ptr) :: gpu_ptr
!Local variables ------------------------------
 type(c_ptr) :: cptr

! *********************************************************************

 if (.false.) cptr=gpu_ptr

end subroutine dealloc_on_gpu
!!***

!!****f* m_lobpcg/gpu_linalg_init
!! NAME
!!  gpu_linalg_init
!!
!! FUNCTION
!!  initialisation of linalg environnement on GPU
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (FDahm)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!   WARNING! : this routine is a dummy one when HAVE_GPU_CUDA is not enabled
!!   the correct one is in 15_gpu_toolbox/gpu_linalg.cu
!!
!! PARENTS
!!      lobpcgwf
!!
!! CHILDREN
!!
!! SOURCE

subroutine gpu_linalg_init()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gpu_linalg_init'
!End of the abilint section

 implicit none

end subroutine gpu_linalg_init
!!***

!!****f* m_lobpcg/gpu_linalg_shutdown
!! NAME
!!  gpu_linalg_shutdown
!!
!! FUNCTION
!!  close linalg environnement on GPU
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (FDahm)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!   WARNING! : this routine is a dummy one when HAVE_GPU_CUDA is not enabled
!!   the correct one is in 15_gpu_toolbox/gpu_linalg.cu
!!
!! PARENTS
!!      lobpcgwf
!!
!! CHILDREN
!!
!! SOURCE
subroutine gpu_linalg_shutdown()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gpu_linalg_shutdown'
!End of the abilint section

 implicit none

end subroutine gpu_linalg_shutdown
!!***

!!****f* m_lobpcg/gpu_xgemm
!! NAME
!!  gpu_xgemm
!!
!! FUNCTION
!!  Compute a scalar-matrix-matrix product and return a scalar-matrix product on GPU
!!  c = alpha * op(a) * op(b) + beta * c
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (FDahm)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cplx  = 1 if real 2 if complex
!!  transa= from of op(a) to be used in the matrix multiplication
!!  transb= from of op(b) to be used in the matrix multiplication
!!  m     = number of rows of the matrix op(a) and of the matrix c
!!  n     = number of rows of the matrix op(b) and the number of columns of the matrix c
!!  k     = number of columns of the matrix op(a) and the number of rows of the matrix op(b)
!!  alpha = alpha scalar coefficient for matrix op(a)
!!  a_gpu = pointer to gpu memory location of  matrix a
!!  lda   = first dimension of a
!!  b_gpu = pointer to gpu memory location of  matrix b
!!  ldb   = first dimension of b
!!  beta  = beta scalar coefficient for matrix c
!!  c_gpu = pointer to gpu memory location of  matrix c
!!  ldc   = first dimension of c
!!
!! OUTPUT
!!  c     = c matrix
!!
!!
!! SIDE EFFECTS
!!   WARNING! : this routine is a dummy one when HAVE_GPU_CUDA is not enabled
!!   the correct one is in 15_gpu_toolbox/gpu_linalg.cu
!!
!! PARENTS
!!      lobpcgwf,m_lobpcg
!!
!! CHILDREN
!!
!! SOURCE

subroutine gpu_xgemm(cplx,transa,transb,m,n,k,alpha,a_gpu,lda,b_gpu,ldb,beta,c_gpu,ldc)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gpu_xgemm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: cplx,lda,ldb,ldc,m,n,k
 complex(dpc),intent(in) :: alpha,beta
 character(len=1),intent(in) :: transa,transb
 type(c_ptr),intent(in) :: a_gpu,b_gpu
 type(c_ptr),intent(inout) :: c_gpu
!Local variables ------------------------------
 type(c_ptr) :: cptr

! *********************************************************************

 if (.false.) then
   cptr=a_gpu;cptr=b_gpu;cptr=c_gpu
   write(std_out,*) transa,transb,cplx,lda,ldb,ldc,m,n,k,alpha,beta
 end if

end subroutine gpu_xgemm
!!***

!!****f* m_lobpcg/gpu_xtrsm
!! NAME
!!  gpu_xtrsm
!!
!! FUNCTION
!! Solves a matrix equation (one matrix operand is triangular) on GPU.
!! The xtrsm routines solve one of the following matrix equations
!! op(a)*x = alpha*b
!! or
!! x*op(a) = alpha*b,
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (FDahm)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cplx= 1 if real 2 if complex
!! side= Specifies whether op(a) appears on the left or right of x for
!!      the operation to be performed as follows:
!!      L or l op(a)*x = alpha*b
!!      R or r x*op(a) = alpha*b
!! uplo= Specifies whether the matrix a is an upper or lower triangular
!!      matrix as follows:
!!      U or u Matrix a is an upper triangular matrix.
!!      L or l Matrix a is a lower triangular matrix
!! transa= Specifies the form of op(a) to be used in the matrix
!!      multiplication as follows:
!!      N or n op(a) = a
!!      T or t op(a) = a'
!!      C or c op(a) = conjg(a')
!! diag= Specifies whether or not a is unit triangular as follows:
!!      U or u Matrix a is assumed to be unit triangular.
!!      N or n Matrix a is not assumed to be unit triangular.
!! m= Specifies the number of rows of b. The value of m must be at least zero
!! n= Specifies the number of columns of b. The value of n must be at least zero
!! alpha= Specifies the scalar alpha. When alpha is zero, then a is not referenced and b
!!      need not be set before entry.
!!  a_gpu = pointer to gpu memory location of  array a, DIMENSION (lda, k), where k is m when side = 'L' or 'l' and is n
!!      when side = 'R' or 'r'.
!! lda= Specifies the first dimension of a as declared in the calling
!!     (sub)program. When side = 'L' or 'l', then lda must be at least max(1,
!!      m), when side = 'R' or 'r', then lda must be at least max(1, n).
!!  b_gpu = pointer to gpu memory location of  b Array, DIMENSION (ldb,n). Before entry, the leading m-by-n part of the array
!!     b must contain the right-hand side matrix b.
!! ldb= Specifies the first dimension of b as declared in the calling
!!     (sub)program. The value of ldb must be at least max(1, m).
!!
!! OUTPUT
!!  b_gpu
!!
!! SIDE EFFECTS
!!   WARNING! : this routine is a dummy one when HAVE_GPU_CUDA is not enabled
!!   the correct one is in 15_gpu_toolbox/gpu_linalg.cu
!!
!! PARENTS
!!      lobpcgwf,m_lobpcg
!!
!! CHILDREN
!!
!! SOURCE
subroutine gpu_xtrsm(cplx,side,uplo,transa,diag,m,n,alpha,a_gpu,lda,b_gpu,ldb)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gpu_xtrsm'
!End of the abilint section

 implicit none

! !Arguments ------------------------------------
 integer, intent(in) :: cplx,lda,ldb,m,n
 complex(dpc), intent(in) :: alpha
 character(len=1), intent(in) :: side,uplo,transa,diag
 type(c_ptr),intent(in) :: a_gpu
 type(c_ptr),intent(inout) :: b_gpu
!Local variables ------------------------------
 type(c_ptr) :: cptr

! *********************************************************************

 if (.false.) then
   cptr=a_gpu;cptr=b_gpu
   write(std_out,*) side,uplo,transa,diag,cplx,lda,ldb,m,n,alpha
 end if

end subroutine gpu_xtrsm
!!***

#endif

END MODULE m_lobpcg
!!***
