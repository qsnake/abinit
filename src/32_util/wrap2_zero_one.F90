!{\src2tex{textfont=tt}}
!!****f* ABINIT/wrap2_zero_one
!! NAME
!! wrap2_zero_one
!!
!! FUNCTION
!! Transforms a real number (num) in its corresponding reduced number
!! (red) in the interval [0,1[ where 1 is not included (tol12)
!! num=red+shift
!!
!! COPYRIGHT
!!  Copyright (C) 2006-2012 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  num=real number
!!
!! OUTPUT
!! red=reduced number of num in the interval [0,1[ where 1 is not included
!! shift=num-red
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      exc_plot,k_neighbors,lin_interpq_gam,m_paw_toolbox,mknesting,outnesting
!!      pawmkaewf
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine wrap2_zero_one(num,red,shift)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wrap2_zero_one'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: num
 real(dp),intent(out) :: red,shift

! *************************************************************************

 if (num>zero) then
   red=mod((num+tol12),one)-tol12
 else
   red=-mod(-(num-one+tol12),one)+one-tol12
 end if
 if(abs(red)<tol12)red=0.0_dp
 shift=num-red

end subroutine wrap2_zero_one
!!***
