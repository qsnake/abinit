!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ifc
!! NAME
!!  m_ifc
!!
!! FUNCTION
!!  This module contains the declaration of data types and methods 
!!  used to handle interatomic force constant sets 
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

MODULE m_ifc

 use m_profiling

 use defs_basis
 use m_errors

 implicit none

 private 
!!***

!!****t* m_ifc/ifc_type
!! NAME
!! ifc_type
!!
!! FUNCTION
!!  Stores data related to interatomic force constants, their real space lattice points
!!
!! SOURCE

 type,public :: ifc_type

  integer :: natom
  ! number of atoms

  integer :: nrpt
  ! Number of r space points in Wigner Seitz cell

  real(dp),pointer :: atmfrc(:,:,:,:,:,:)   SET2NULL
  ! inter atomic forces in real space (2,3,natom,3,natom,nrpt)

  real(dp), pointer :: rpt(:,:) SET2NULL
  ! real space points in canonical type coordinates (3,nrpt)

  real(dp),pointer :: wghatm(:,:,:)   SET2NULL
  ! Weights for each point and atom in the Wigner Seitz supercell in real space (natom,natom,nrpt)

 end type ifc_type
!!***
 
 ! Bound methods:
 public :: nullify_ifc
 public :: destroy_ifc

! should add mkifc9 into this module eventually, as an initializer

!----------------------------------------------------------------------

CONTAINS  !===========================================================
!!***

!!****f* m_ifc/nullify_ifc
!! NAME
!!
!! FUNCTION
!!  Initialize pointer types for the ifc_type structure, by nullifying them
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine nullify_ifc(ifc)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_ifc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(ifc_type),intent(out) :: ifc

! ************************************************************************

 nullify(ifc%atmfrc)
 nullify(ifc%rpt)
 nullify(ifc%wghatm)

end subroutine nullify_ifc
!!***

!----------------------------------------------------------------------

!!****f* m_ifc/destroy_ifc
!! NAME
!!
!! FUNCTION
!!  Clean and deallocate pointer types for the ifc_type structure
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_ifc(ifc)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_ifc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(ifc_type),intent(inout) :: ifc

! ************************************************************************
 ! real
 if (associated(ifc%atmfrc   ))  then
   ABI_DEALLOCATE(ifc%atmfrc)
 end if
 if (associated(ifc%rpt      ))  then
   ABI_DEALLOCATE(ifc%rpt)
 end if
 if (associated(ifc%wghatm   ))  then
   ABI_DEALLOCATE(ifc%wghatm)
 end if

end subroutine destroy_ifc
!!***

END MODULE m_ifc
!!***
