!{\src2tex{textfont=tt}}
!!****f* ABINIT/initang
!! NAME
!! initang
!!
!! FUNCTION
!! Initialize angular mesh for PAW calculations
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!       pawang%angl_size  - Total number of sample points in the angular mesh
!!       pawang%ntheta     - Number of sample points in the theta dir
!!       pawang%nphi       - Number of sample points in the phi dir
!!
!! OUTPUT
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!       pawang%anginit    - (3 x angl_size) array, the ntheta*nphi
!!                           dimensional arrays ax, ay, and az
!!       pawang%angwgth    - (angl_size) array, the weight factor of the
!!                           point (ax, ay, az)
!!
!! PARENTS
!!      pawinit
!!
!! CHILDREN
!!      coeffs_gausslegint
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine initang(pawang)

 use m_profiling

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initang'
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(pawang_type),intent(inout) :: pawang

!Local variables-------------------------------
!scalars
 integer :: ip,it,npoints
 real(dp) :: ang,con,cos_phi,cos_theta,sin_phi,sin_theta
 character(len=500) :: msg
!arrays
 real(dp) :: th(pawang%ntheta),wth(pawang%ntheta)

! ***********************************************************************

 DBG_ENTER("COLL")

 if (pawang%angl_size==0) return

!Initializations
 npoints=0
 con=two_pi / pawang%nphi
 call coeffs_gausslegint(-one,one,th,wth,pawang%ntheta)

!We now open two nested do-loops. The first loops through the number
!of theta angles, the second through the number of phi angles (?).
!The two together initialize anginit.

 do it = 1, pawang%ntheta

   cos_theta = th(it)
   sin_theta = sqrt(one - cos_theta*cos_theta)

   do ip = 1, pawang%nphi

     ang = con * (ip-1)
     cos_phi = cos(ang)
     sin_phi = sin(ang)

     npoints = npoints + 1

     pawang%anginit(1, npoints) = sin_theta * cos_phi
     pawang%anginit(2, npoints) = sin_theta * sin_phi
     pawang%anginit(3, npoints) = cos_theta

!    Normalization required
     pawang%angwgth(npoints) = wth(it) / (2 * pawang%nphi)

   end do
 end do

!The following is an error statement that will be generated
!if npoints exceeds nang...
 if (npoints > pawang%angl_size) then
   write(msg, '(a,i4,a,a,i4)' ) &
&   '  anginit%npoints =',npoints,ch10,&
&   '        angl_size =',pawang%angl_size
   MSG_BUG(msg)
 end if

 DBG_EXIT("COLL")

end subroutine initang
!!***
