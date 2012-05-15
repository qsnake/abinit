!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ddb_blk
!! NAME
!!  m_ddb_blk
!!
!! FUNCTION
!!  This module contains the declaration of data types and methods 
!!  used to handle the blocks of data in DDB files:
!!  blkval, nrm, qpt, flg, and associated dimensions
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!! 
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_ddb_blk

 use m_profiling

 use defs_basis
 use m_errors

 implicit none

 private 
!!***

!!****t* m_ddb_blk/ddb_blk_type
!! NAME
!! ddb_blk_type
!!
!! FUNCTION
!!  Stores data related to derivative database contents
!!
!! SOURCE

 type,public :: ddb_blk_type

  integer :: msize
  ! maximum size of dynamical matrices and other perturbations (ddk, dde...)

  integer :: nblok
  ! number of 2dte blocks in present object

  integer,pointer :: flg(:,:)   SET2NULL
  ! flg(msize,nblok)
  ! flag to indicate presence of a given block
 
  integer,pointer :: typ(:)   SET2NULL
  ! typ(nblok)
  ! type of each block - ddk dde, phonon etc...

  real(dp),pointer :: nrm(:,:)   SET2NULL
  ! nrm(3,nblok)
  ! norm of the q-points for each block - can be 0 to indicate a direction of approach to gamma

  real(dp),pointer :: qpt(:,:)   SET2NULL
  ! qpt(9,nblok)
  ! q-point vector in reciprocal space (reduced lattice coordinates) for each block

  real(dp),pointer :: val(:,:,:)   SET2NULL
  ! val(2,msize,nblok)
  ! values of the second energy derivatives in each block

 end type ddb_blk_type
!!***
 
 ! Bound methods:
 public :: create_ddb_blk
 public :: nullify_ddb_blk
 public :: destroy_ddb_blk

!----------------------------------------------------------------------

CONTAINS  !===========================================================
!!***

!!****f* m_ddb_blk/nullify_ddb_blk
!! NAME
!!
!! FUNCTION
!!  Initialize pointer types for the ddb_blk_type structure, by nullifying them
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine nullify_ddb_blk(ddb_blk)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_ddb_blk'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(ddb_blk_type),intent(out) :: ddb_blk

! ************************************************************************

 nullify(ddb_blk%val)
 nullify(ddb_blk%nrm)
 nullify(ddb_blk%flg)
 nullify(ddb_blk%qpt)
 nullify(ddb_blk%typ)

end subroutine nullify_ddb_blk
!!***

!----------------------------------------------------------------------

!!****f* m_ddb_blk/destroy_ddb_blk
!! NAME
!!
!! FUNCTION
!!  Clean and deallocate pointer types for the ddb_blk_type structure
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_ddb_blk(ddb_blk)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_ddb_blk'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(ddb_blk_type),intent(inout) :: ddb_blk

! ************************************************************************
 ! real
 if (associated(ddb_blk%val))  then
   ABI_DEALLOCATE(ddb_blk%val)
 end if
 if (associated(ddb_blk%nrm))  then
   ABI_DEALLOCATE(ddb_blk%nrm)
 end if
 if (associated(ddb_blk%qpt))  then
   ABI_DEALLOCATE(ddb_blk%qpt)
 end if
 if (associated(ddb_blk%typ))  then
   ABI_DEALLOCATE(ddb_blk%typ)
 end if
 if (associated(ddb_blk%flg))  then
   ABI_DEALLOCATE(ddb_blk%flg)
 end if

end subroutine destroy_ddb_blk
!!***

!----------------------------------------------------------------------

!!****f* m_ddb_blk/create_ddb_blk
!! NAME
!!
!! FUNCTION
!!  Allocate pointer types for the ddb_blk_type structure
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine create_ddb_blk(msize, nblok, ddb_blk)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'create_ddb_blk'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 integer, intent(in) :: msize
 integer, intent(in) :: nblok
 type(ddb_blk_type),intent(inout) :: ddb_blk

! ************************************************************************

 ! in case someone has forgotten to deallocate stuff
 call destroy_ddb_blk(ddb_blk)
 
 ddb_blk%msize = msize
 ddb_blk%nblok = nblok
 ABI_ALLOCATE(ddb_blk%val,(2,msize,nblok))
 ABI_ALLOCATE(ddb_blk%nrm,(3,nblok))
 ABI_ALLOCATE(ddb_blk%qpt,(9,nblok))
 ABI_ALLOCATE(ddb_blk%typ,(nblok))
 ABI_ALLOCATE(ddb_blk%flg,(msize,nblok))

end subroutine create_ddb_blk
!!***

!----------------------------------------------------------------------

!!****f* m_ddb_blk/copy_ddb_blk
!! NAME
!!
!! FUNCTION
!!  create object and copy all types for the ddb_blk_type structure
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine copy_ddb_blk(blk_in, blk_out)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'copy_ddb_blk'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(ddb_blk_type),intent(in) :: blk_in
 type(ddb_blk_type),intent(out) :: blk_out

! ************************************************************************

 call create_ddb_blk(blk_in%msize, blk_in%nblok, blk_out)

 blk_out%val = blk_in%val
 blk_out%nrm = blk_in%nrm
 blk_out%qpt = blk_in%qpt
 blk_out%typ = blk_in%typ
 blk_out%flg = blk_in%flg

end subroutine copy_ddb_blk
!!***


END MODULE m_ddb_blk
!!***
