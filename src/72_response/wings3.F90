!{\src2tex{textfont=tt}}
!!****f* ABINIT/wings3
!! NAME
!! wings3
!!
!! FUNCTION
!!  Suppress the wings of the cartesian 2DTE for which
!!  the diagonal element is not known
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  carflg(3,mpert,3,mpert)= ( 1 if the element of the cartesian
!!  2DTE matrix has been calculated correctly ; 0 otherwise )
!!  d2cart(2,3,mpert,3,mpert)=
!!   dynamical matrix, effective charges, dielectric tensor,....
!!   all in cartesian coordinates
!!  mpert =maximum number of ipert
!!
!! OUTPUT
!!  d2cart(2,3,mpert,3,mpert) without the wings
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine wings3(carflg,d2cart,mpert)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wings3'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert
!arrays
 integer,intent(inout) :: carflg(3,mpert,3,mpert)
 real(dp),intent(inout) :: d2cart(2,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: idir,idir1,ipert,ipert1

! *********************************************************************

 do ipert=1,mpert
   do idir=1,3
     if(carflg(idir,ipert,idir,ipert)==0)then
       do ipert1=1,mpert
         do idir1=1,3
           carflg(idir,ipert,idir1,ipert1)=0
           carflg(idir1,ipert1,idir,ipert)=0
           d2cart(1,idir,ipert,idir1,ipert1)=0.0_dp
           d2cart(2,idir,ipert,idir1,ipert1)=0.0_dp
           d2cart(1,idir1,ipert1,idir,ipert)=0.0_dp
           d2cart(2,idir1,ipert1,idir,ipert)=0.0_dp
         end do
       end do
     end if
   end do
 end do

end subroutine wings3
!!***
