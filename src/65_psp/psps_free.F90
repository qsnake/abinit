!{\src2tex{textfont=tt}}
!!****f* ABINIT/psps_free
!! NAME
!! psps_free
!!
!! FUNCTION
!! Deallocate all memory of psps structure.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (DC)
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
!! psps=<type pseudopotential_type>the pseudopotentials description
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      psp2params_free
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psps_free(psps)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psps_free'
 use interfaces_65_psp, except_this_one => psps_free
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(pseudopotential_type),intent(inout) :: psps
!arrays

!Local variables-------------------------------
!scalars

! *************************************************************************

!Allocation of some arrays independent of the dataset
 ABI_DEALLOCATE(psps%filpsp)
 ABI_DEALLOCATE(psps%pspcod)
 ABI_DEALLOCATE(psps%pspdat)
 ABI_DEALLOCATE(psps%pspso)
 ABI_DEALLOCATE(psps%pspxc)
 ABI_DEALLOCATE(psps%title)
 ABI_DEALLOCATE(psps%zionpsp)
 ABI_DEALLOCATE(psps%znuclpsp)
 call psp2params_free(psps%gth_params)

 ABI_DEALLOCATE(psps%algalch)
 ABI_DEALLOCATE(psps%mixalch)

 ABI_DEALLOCATE(psps%ekb)
 ABI_DEALLOCATE(psps%indlmn)
 ABI_DEALLOCATE(psps%ffspl)
 ABI_DEALLOCATE(psps%qgrid_ff)
 ABI_DEALLOCATE(psps%qgrid_vl)
 ABI_DEALLOCATE(psps%vlspl)
 if (.not.psps%vlspl_recipSpace) then
   ABI_DEALLOCATE(psps%dvlspl)
 end if
 ABI_DEALLOCATE(psps%xccc1d)
 ABI_DEALLOCATE(psps%xcccrc)
 ABI_DEALLOCATE(psps%ziontypat)
 ABI_DEALLOCATE(psps%znucltypat)

end subroutine psps_free
!!***
