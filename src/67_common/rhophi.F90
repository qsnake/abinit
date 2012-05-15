!{\src2tex{textfont=tt}}
!!****f* ABINIT/rhophi
!! NAME
!! rhophi
!!
!! FUNCTION
!! Compute the phase and the module of a complex number.
!! The phase angle is fold into the interval [-pi,pi]
!!
!! COPYRIGHT
!! Copyright (C) 2000-2012 ABINIT  group (MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cx(2) = complex number
!!
!! OUTPUT
!!  phi = phase of cx fold into [-pi,pi]
!!  rho = modul of cx
!!
!! PARENTS
!!      berryphase_new,etheta,linemin
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine rhophi(cx,phi,rho)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhophi'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(out) :: phi,rho
!arrays
 real(dp),intent(in) :: cx(2)

!Local variables-------------------------------

! ***********************************************************************


 rho = sqrt(cx(1)*cx(1) + cx(2)*cx(2))

 if (abs(cx(1)) > tol8) then

   phi = atan(cx(2)/cx(1))

!  phi is an element of [-pi,pi]
   if (cx(1) < zero) then
     if (phi < zero) then
       phi = phi + pi
     else
       phi = phi - pi
     end if
   end if

 else

   if (cx(2) > tol8) then
     phi = pi*half
   else if (cx(2) < tol8) then
     phi = -0.5_dp*pi
   else
     phi = 0
   end if

 end if

end subroutine rhophi
!!***
