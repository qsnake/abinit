!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_delocint
!!
!! NAME
!! m_delocint
!!
!! FUNCTION
!! Module for delocalized internal coordinates: container type, and nullify+destruction routines
!!
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUT
!!
!! OUTPUT
!!
!! NOTES
!! 
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_delocint

 use m_profiling

 use defs_basis
 use m_errors

 implicit none
!!***

!!****t* m_delocint/ab_delocint
!! NAME
!! ab_delocint
!!
!! FUNCTION
!! Datatype with the important variables in pred_delocint
!!
!! NOTES
!!
!!   deloc <type(ab_delocint)>=Important variables for
!!   |                           pred_delocint
!!   |
!!   | nang     = Number of angles
!!   | nbond    = Number of bonds
!!   | ncart    = Number of cartesian directions
!!   |             (used for constraints)
!!   | ndihed   = Number of dihedrals
!!   | nrshift  = Dimension of rshift
!!   | ninternal= Number of internal coordinates
!!   |            ninternal=nbond+nang+ndihed+ncart
!!   |
!!   | angs(2,3,nang)  = Indexes to characterize angles
!!   | bonds(2,2,nbond)= For a bond between iatom and jatom
!!   |                   bonds(1,1,nbond) = iatom
!!   |                   bonds(2,1,nbond) = icenter
!!   |                   bonds(1,2,nbond) = jatom
!!   |                   bonds(2,2,nbond) = irshift
!!   | carts(2,ncart)  = Index of total primitive internal,
!!   |                   and atom (carts(2,:))
!!   | dihedrals(2,4,ndihed)= Indexes to characterize dihedrals
!!   |
!!   | rshift(3,nrshift)= Shift in xred that must be done to find
!!   |                    all neighbors of a given atom within a
!!   |                    given number of neighboring shells
!!
!! SOURCE

type ab_delocint

! scalars
 integer :: nang
 integer :: nbond
 integer :: ncart
 integer :: ndihed
 integer :: nrshift
 integer :: ninternal

! arrays
 integer,pointer :: angs(:,:,:)
 integer,pointer :: bonds(:,:,:)
 integer,pointer :: carts(:,:)
 integer,pointer :: dihedrals(:,:,:)
 real(dp),pointer :: rshift(:,:)

end type ab_delocint
!!***


contains

!----------------------------------------------------------------------

!!****f* m_delocint/nullify_delocint
!!
!! NAME
!! nullify_delocint
!!
!! FUNCTION
!! nullify function for delocalized internals object
!!
!! INPUT
!! deloc= container object for delocalized internal coordinates
!!
!! OUTPUT
!!
!! PARENTS
!!      defs_mover
!!
!! CHILDREN
!!
!! SOURCE

subroutine nullify_delocint(deloc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_delocint'
!End of the abilint section

  type(ab_delocint), intent(out) :: deloc
 
  nullify(deloc%angs)
  nullify(deloc%bonds)
  nullify(deloc%carts)
  nullify(deloc%dihedrals)
  nullify(deloc%rshift)

end subroutine nullify_delocint
!!***

!----------------------------------------------------------------------

!!****f* m_delocint/destroy_delocint
!!
!! NAME
!! destroy_delocint
!!
!! FUNCTION
!! destructor function for delocint object
!!
!! INPUT
!! deloc= container object for delocalized internal coordinates
!!
!! OUTPUT
!!
!! PARENTS
!!      defs_mover
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_delocint(deloc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_delocint'
!End of the abilint section

  type(ab_delocint), intent(out) :: deloc
 
  if(associated(deloc%angs))  then
    ABI_DEALLOCATE(deloc%angs)
  end if
  if(associated(deloc%bonds))  then
    ABI_DEALLOCATE(deloc%bonds)
  end if
  if(associated(deloc%carts))  then
    ABI_DEALLOCATE(deloc%carts)
  end if
  if(associated(deloc%dihedrals))  then
    ABI_DEALLOCATE(deloc%dihedrals)
  end if
  if(associated(deloc%rshift))  then
    ABI_DEALLOCATE(deloc%rshift)
  end if

end subroutine destroy_delocint
!!***

end module m_delocint
!!***
