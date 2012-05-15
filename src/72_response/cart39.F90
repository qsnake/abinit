!{\src2tex{textfont=tt}}
!!****f* ABINIT/cart39
!! NAME
!! cart39
!!
!!
!! FUNCTION
!! Transform a vector from reduced coordinates to cartesian coordinates,
!! taking into account the perturbation from which it was derived,
!! and also check the existence of the new values.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2012 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  flg1(3)=tell if information of each component of vec1
!!        is valid
!!  gprimd(3,3)=basis vector in the reciprocal space
!!  ipert=number of the perturbation
!!  natom=number of atom
!!  rprimd(3,3)=basis vector in the real space
!!  vec1(3)=input vector, in reduced coordinates
!!
!! OUTPUT
!!  flg2(3)=tell if information of each component of vec2 is valid
!!  vec2(3)=output vector, in cartesian coordinates
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!!
!! PARENTS
!!      cart29,carteig2d,gath3,nlopt
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine cart39(flg1,flg2,gprimd,ipert,natom,rprimd,vec1,vec2)

 use m_profiling

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cart39'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: ipert,natom
!arrays
 integer,intent(in) :: flg1(3)
 integer,intent(out) :: flg2(3)
 real(dp),intent(in) :: gprimd(3,3),rprimd(3,3),vec1(3)
 real(dp),intent(out) :: vec2(3)

!Local variables -------------------------
!scalars
 integer :: idir,ii

! *********************************************************************

!Treat phonon-type perturbation
 if(ipert>=1.and.ipert<=natom)then

   do idir=1,3
     vec2(idir)=0.0_dp
     flg2(idir)=1
     do ii=1,3
       if(abs(gprimd(idir,ii))>1.0d-10)then
         if(flg1(ii)==1)then
           vec2(idir)=vec2(idir)+gprimd(idir,ii)*vec1(ii)
         else
           flg2(idir)=0
         end if
       end if
     end do
     if(flg2(idir)==0)vec2(idir)=0.0_dp
   end do

!  Treat electric field perturbation
 else if(ipert==natom+2) then
!  OCL SCALAR
   do idir=1,3
     vec2(idir)=0.0_dp
     flg2(idir)=1
!    OCL SCALAR
     do ii=1,3
       if(abs(rprimd(idir,ii))>1.0d-10)then
         if(flg1(ii)==1)then
           vec2(idir)=vec2(idir)+rprimd(idir,ii)*vec1(ii)/two_pi
         else
           flg2(idir)=0
         end if
       end if
     end do
     if(flg2(idir)==0)vec2(idir)=0.0_dp
   end do

!  Treat other perturbations
 else
   do idir=1,3
     vec2(idir)=vec1(idir)
     flg2(idir)=flg1(idir)
   end do
 end if

end subroutine cart39
!!***
