!{\src2tex{textfont=tt}}
!!****m* ABINIT/mod_prc_memory
!! NAME
!! mod_prc_memory
!!
!! FUNCTION
!! This modules defines arrays and data used for the real-space kerker
!! preconditionning of potential residuals.
!! 
!! COPYRIGHT
!! Copyright (C) 2009-2012 ABINIT group (PMA).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  FIXME: this is highly non-kosher. Should be a datastructure which is declared dynamically
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module mod_prc_memory

 use m_profiling

 use defs_basis

 implicit none

!Arguments -------------------------------
!Local variables -------------------------

private

  real(dp),public, allocatable :: rdiemac(:)

  integer, public :: cycle=0
  real(dp), public :: energy_min

public :: prc_mem_init
public :: prc_mem_free

! *********************************************************************

 contains
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/prc_mem_init
!! NAME
!! prc_mem_init
!!
!! FUNCTION
!! This subroutine allocates the module's main component
!! 
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (MJV).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      prcrskerker1
!!
!! CHILDREN
!!
!! SOURCE

subroutine prc_mem_init(nfft)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prc_mem_init'
!End of the abilint section

implicit none

!Arguments -------------------------------
integer, intent(in) :: nfft
!Local variables -------------------------
! *********************************************************************

   if (.not. allocated(rdiemac))  then
     ABI_ALLOCATE(rdiemac,(nfft))
   end if
   if(nfft.ne.size(rdiemac)) then ! This steps should be done over "istep" instead
     ABI_DEALLOCATE(rdiemac)
     ABI_ALLOCATE(rdiemac,(nfft))
     cycle=0
   end if

 end subroutine prc_mem_init
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/prc_mem_free
!! NAME
!! prc_mem_free
!!
!! FUNCTION
!! This subroutine deallocates the module's main component
!! 
!! COPYRIGHT
!! Copyright (C) 2011-2012 ABINIT group (MJV).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!
!! SOURCE

subroutine prc_mem_free()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prc_mem_free'
!End of the abilint section

implicit none

!Arguments -------------------------------

!Local variables -------------------------

! *********************************************************************

   if (allocated(rdiemac))  then
     ABI_DEALLOCATE(rdiemac)
   end if

 end subroutine prc_mem_free
!!***



end module mod_prc_memory
!!***
