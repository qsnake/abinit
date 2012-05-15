!{\src2tex{textfont=tt}}
!!****f* ABINIT/evspln
!! NAME
!! evspln
!!
!! FUNCTION
!! Evaluation of the value of a function, given its spline coefficients,
!! in an interval between two points.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2012 ABINIT group (PCasek,FF,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  aa,bb = coefficients of the local cubic spline
!!  dxx   = difference between the end points of the interval
!!  fld(ndim)=value of the function at all points
!!  indx  = index of the interval
!!  kod   = option for actual numerical evaluation of val, der and dder
!!         (kod is expected between 0 and 7)
!!  ndim  = dimension of indx and sdfd
!!  sdfd(ndim)=value of the second derivative of the function at all points
!!
!! OUTPUT
!!  val= value of the function (if kod is odd)
!!  der=derivative of the function (if kod is 2, 3, 6 or 7)
!!  dder=second derivative of the function (if kod is >=4)
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


subroutine evspln(aa,bb,indx,dxx,ndim,fld,sdfd,kod,val,der,dder)

 use m_profiling

 use defs_basis
 use defs_aimfields

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'evspln'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: indx,kod,ndim
 real(dp),intent(in) :: aa,bb,dxx
 real(dp),intent(out) :: dder,der,val
!arrays
 real(dp),intent(in) :: fld(ndim),sdfd(ndim)

!Local variables-------------------------------
!scalars

! *************************************************************************

 if (indx /= ndim) then
   if (btest(kod,0)) then
     val=aa*fld(indx)+bb*fld(indx+1)+(aa**3-aa)*dxx**2*sdfd(indx)/6._dp &
&     +(bb**3-bb)*dxx**2*sdfd(indx+1)/6._dp
   end if
   if (btest(kod,1)) then
     der=(fld(indx+1)-fld(indx))/dxx-(3*aa**2-1._dp)*dxx*sdfd(indx)/6._dp &
&     +(3*bb**3-1)*dxx*sdfd(indx+1)/6._dp
   end if
   if (btest(kod,2)) then
     dder=aa*sdfd(indx)+bb*sdfd(indx+1)
   end if
 else
   if (btest(kod,0)) then
     val=aa*fld(indx)+bb*fld(1)+(aa**3-aa)*dxx**2*sdfd(indx)/6._dp &
&     +(bb**3-bb)*dxx**2*sdfd(1)/6._dp
   end if
   if (btest(kod,1)) then
     der=(fld(1)-fld(indx))/dxx-(3*aa**2-1._dp)*dxx*sdfd(indx)/6._dp &
&     +(3*bb**3-1)*dxx*sdfd(1)/6._dp
   end if
   if (btest(kod,2)) then
     dder=aa*sdfd(indx)+bb*sdfd(1)
   end if
 end if

end subroutine evspln
!!***
