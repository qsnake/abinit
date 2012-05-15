!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp2params_init
!! NAME
!! psp2params_init
!!
!! FUNCTION
!! Allocate and initialise the data structure holding parameters for the GTH
!! pseudo-potentials.
!!
!!  MJV note: this should be renamed: psp2 suggests it relates to pspcod 2,
!!     whereas it is actually 3 
!!    the parameters would also be better off separated into C and h arrays
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  npsp=number of true pseudo used (not alchemy).
!!
!! OUTPUT
!!  gth_params <type (pseudopotential_gth_type)>=the values to allocate and initialise.
!!
!! PARENTS
!!      psps_init_global
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp2params_init(gth_params, npsp)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp2params_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npsp
 type(pseudopotential_gth_type),intent(out) :: gth_params

!Local variables-------------------------------

! *********************************************************************
!Check array, no params are currently set.
 ABI_ALLOCATE(gth_params%set,(npsp))
 gth_params%set(:) = .false.

!Check array, have geometric informations been filled?
 ABI_ALLOCATE(gth_params%hasGeometry,(npsp))
 gth_params%hasGeometry(:) = .false.

!Coefficients for local part and projectors
 ABI_ALLOCATE(gth_params%psppar,(0:4, 0:6, npsp))
 gth_params%psppar = real(0, dp)

!Coefficients for spin orbit part
 ABI_ALLOCATE(gth_params%psp_k_par,(1:4, 1:3, npsp))
 gth_params%psp_k_par = zero

!Different radii
 ABI_ALLOCATE(gth_params%radii_cov,(npsp))
 ABI_ALLOCATE(gth_params%radii_cf,(npsp, 3))

!Number of semicore electrons
 ABI_ALLOCATE(gth_params%semicore,(npsp))
end subroutine psp2params_init
!!***

!!****f* ABINIT/psp2params_free
!! NAME
!! psp2params_free
!!
!! FUNCTION
!! Deallocate a previously allocated data structure for storage of GTH parameters.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!  gth_params <type (pseudopotential_gth_type)>=the values to deallocate.
!!
!! PARENTS
!!      psps_free
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine psp2params_free(gth_params)

 use m_profiling

 use defs_basis
 use defs_datatypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp2params_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(pseudopotential_gth_type),intent(inout) :: gth_params

!Local variables-------------------------------

! *********************************************************************

!Check arrays.
 ABI_DEALLOCATE(gth_params%set)
 ABI_DEALLOCATE(gth_params%hasGeometry)

!Coefficients for local part and projectors
 ABI_DEALLOCATE(gth_params%psppar)

!Coefficients for spin orbit part
 ABI_DEALLOCATE(gth_params%psp_k_par)

!Different radii
 ABI_DEALLOCATE(gth_params%radii_cov)
 ABI_DEALLOCATE(gth_params%radii_cf)

!Number of semicore electrons
 ABI_DEALLOCATE(gth_params%semicore)

end subroutine psp2params_free
!!***
